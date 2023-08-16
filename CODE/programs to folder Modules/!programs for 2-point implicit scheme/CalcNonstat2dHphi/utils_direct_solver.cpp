/*
 * GENERAL REMARKS
 *  
 *  This code is freely available under the following conditions:
 *  
 *  1) The code is to be used only for non-commercial purposes.
 *  2) No changes and modifications to the code without prior permission of the developer.
 *  3) No forwarding the code to a third party without prior permission of the developer.
 *  
 *  			TDEMLineCalc
 *  This file contains code for utility functions for calculating with source grouping
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 January 16, 2023
*/

#include "stdafx.h"
#include "utils_direct_solver.h"

extern ofstream logfile;

int GetNumberOfPlaces(int &npls,int &nplsa)
{
	int i,k,p1,p2,nprof;
	ifstream inf;


	npls=0;
	inf.open("clcnplsh");
	if(!inf)
	{
		logfile<<"Error in open file "<<"clcnplsh"<<endl;
		cout<<"Error in open file "<<"clcnplsh"<<endl;
		return 1;
	}
	inf>>npls;
	inf.close();
	inf.clear();

	nplsa=0;
	inf.open("clcnplsa");
	if(!inf)
	{
		logfile<<"Error in open file "<<"clcnplsa"<<endl;
		cout<<"Error in open file "<<"clcnplsa"<<endl;
		return 1;
	}
	inf>>nplsa;
	inf.close();
	inf.clear();

	return 0;
}

int FindTimeLayer(int nt,double *times,double val)
{
	int it;

	if(nt<1)
	{
		return -1;
	}

	if(fabs(val)<1e-14 && fabs(times[0])<1e-14)
	{
		return 0;
	}

	for(it=1;it<nt;it++)
	{
		if(fabs(1.0-times[it]/val)<1e-6)
		{
			return it;
		}
	}

	return -1;
}

int read_time_mesh(int &ntime,double *(&time))
{
	int i;
	ifstream inf;
	inf.open("infite.0");
	if(!inf){return 1;}
	inf.ignore(1000, '=');
	inf >> ntime;
	time=new double[ntime];
	inf.ignore(1000, '\n');
	inf.ignore(1000, '\n');
	inf.ignore(1000, '\n');
	inf.ignore(1000, '\n');
	for (i=0; i<ntime; i++)
	{
		if(!inf.good()){return 1;}
		inf >> time[i];
		inf.ignore(1000, ';');
	}
	inf.close();
	inf.clear();
	return 0;
}

int ReadField2d(char *dpath,double *u,int n,int ind)
{
	FILE *fp;
	char buff[256];
	if(ind>=0)
	{
		sprintf(buff,"%s\\v2.%d",dpath,ind);
		fp=fopen(buff,"rb");
		if(!fp){return 1;}
		fread(u,sizeof(double),n,fp);
		fclose(fp);
		fflush(fp);
	}
	return 0;
}

int ReadField2d(char *dpath,double *u,int n,char *fname)
{
	FILE *fp;
	char buff[256];
	sprintf(buff,"%s\\%s",dpath,fname);
	fp=fopen(buff,"rb");
	if(!fp){return 1;}
	cout<<"Reading file "<<fname<<endl;
	fread(u,sizeof(double),n,fp);
	fclose(fp);
	fflush(fp);
	return 0;
}

int WriteField2d(char *dpath,double *u,int n,int ind)
{
	FILE *fp;
	char buff[256];
	if(ind>=0)
	{
		sprintf(buff,"%s\\v2.%d",dpath,ind);
		fp=fopen(buff,"wb");
		if(!fp){return 1;}
		fwrite(u,sizeof(double),n,fp);
		fclose(fp);
		fflush(fp);
	}
	return 0;
}

int WriteField2d(char *dpath,double *u,int n,char *fname)
{
	FILE *fp;
	char buff[256];
	cout<<"Writing file "<<fname<<endl;
	sprintf(buff,"%s\\%s",dpath,fname);
	fp=fopen(buff,"wb");
	if(!fp){return 1;}
	fwrite(u,sizeof(double),n,fp);
	fclose(fp);
	fflush(fp);
	return 0;
}

int ReadDecInfo(int &nDec,int &iDec,vector<int> &DecIg)
{
	ifstream inf;
	int i;

	inf.open("iDec");
	if(!inf)
	{
		logfile<<"Error in open file "<<"iDec"<<endl;
		cout<<"Error in open file "<<"iDec"<<endl;
		return 1;
	}
	inf>>iDec;
	inf.close();
	inf.clear();

	inf.open("nDec");
	if(!inf)
	{
		logfile<<"Error in open file "<<"nDec"<<endl;
		cout<<"Error in open file "<<"nDec"<<endl;
		return 1;
	}
	inf>>nDec;
	inf.close();
	inf.clear();

	cout<<"iDec= "<<iDec<<endl;
	cout<<"nDec= "<<nDec<<endl;

	DecIg.resize(nDec+1);

	inf.open("DecIg");
	if(!inf)
	{
		logfile<<"Error in open file "<<"DecIg"<<endl;
		cout<<"Error in open file "<<"DecIg"<<endl;
		return 1;
	}
	for(i=0;i<(nDec+1);i++){inf>>DecIg[i];}
	inf.close();
	inf.clear();

	return 0;
}

int FindTimeLayersForDec(int ntime,double *time,int iDec,vector<int> &DecIg,int &npre1,int &npre2,int &tbeg,int &tend)
{
	npre1=-1;
	npre2=-1;
	if(iDec)
	{
		npre1=FindTimeLayer(ntime,time,time[DecIg[iDec]]-1*(time[DecIg[iDec]+1]-time[DecIg[iDec]]));
		npre2=FindTimeLayer(ntime,time,time[DecIg[iDec]]-2*(time[DecIg[iDec]+1]-time[DecIg[iDec]]));
		if(npre1<0)
		{
			cout<<"Error in npre1"<<endl;
			logfile<<"Error in npre1"<<endl;
			return 1;
		}
		if(npre2<0)
		{
			cout<<"Error in npre2"<<endl;
			logfile<<"Error in npre2"<<endl;
			return 1;
		}
	}

	cout<<"npre1= "<<npre1<<endl;
	cout<<"npre2= "<<npre2<<endl;

	tbeg=DecIg[iDec];
	tend=DecIg[iDec+1];

	cout<<"tbeg= "<<tbeg<<endl;
	cout<<"tend= "<<tend<<endl;

	return 0;
}

void GenerateDecInfo(int ntime,double *time)
{
	int iDec,nDec,i,j;
	double htpre,htcur;
	int max_dec_size;
	ofstream ofp;
	vector<int> DecIg;

	htpre=time[1]-time[0];
	DecIg.push_back(0);
	for(i=2;i<ntime;i++)
	{
		htcur=time[i]-time[i-1];
		if(fabs(1.0-htpre/htcur)>1e-3)
		{
			DecIg.push_back(i);
			htpre=htcur;
		}
	}
	DecIg.push_back(ntime);
	
	nDec=(int)DecIg.size()-1;

	ofp.open("nDec");
	ofp<<nDec<<'\n';
	ofp.close();
	ofp.clear();

	ofp.open("DecIg");
	ofp<<0<<'\n';
	for(iDec=0;iDec<nDec;iDec++){ofp<<DecIg[iDec+1]<<'\n';}
	ofp.close();
	ofp.clear();

	max_dec_size=0;
	for(iDec=0;iDec<nDec;iDec++)
	{
		j=DecIg[iDec+1]-DecIg[iDec];
		if(j>max_dec_size)
		{
			max_dec_size=j;
		}
	}
	ofp.open("max_dec_size");
	ofp<<max_dec_size<<'\n';
	ofp.close();
	ofp.clear();
}
