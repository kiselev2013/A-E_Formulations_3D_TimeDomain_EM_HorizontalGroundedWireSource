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
 *  This file contains code for finding nodes to calculate the grounded radial source part of primary field

 *  
 *  
 *  Written by Ph.D. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 January 16, 2023
*/


#include "stdafx.h"
#include "task2d.h"

struct LineAB
{
	PointXYZ A,B;
};

ofstream logfile;

double get_angle(double px,double py,double ox,double oy)
{		
	double x[2],r;
	x[0]=px-ox;
	x[1]=py-oy;
	r=sqrt(x[0]*x[0]+x[1]*x[1]);
	if(!r)return 0;
	x[0]/=r;
	x[1]/=r;
	if(!x[1]){
		if(x[0]>0)return 0;
		else return PI;
	}
	if(!x[0]){
		if(x[1]>0)return PI/2;
		else return 3*PI/2;
	}
	if(x[0]>0 && x[1]>0)return atan2(x[1],x[0]);
	if(x[0]<0 && x[1]>0)return PI-atan2(x[1],-x[0]);
	if(x[0]<0 && x[1]<0)return PI+atan2(-x[1],-x[0]);
	return 2*PI-atan2(-x[1],x[0]);
}

void get_int(char *str, int &ipos, int &var)
{
	int l;
	bool flag;
	char ch;
	char word[32];
	flag = true;
	l = 0;
	do
	{
		ch = str[ipos];
		if (!(ch == ' ' || ch == '\t' || ch == '\n'))
		{
			word[l] = ch;
			l++;
		}
		else
		{
			if (l)
			{
				flag = false;
				word[l] = '\0';
			}
		}
		ipos++;
	} while (flag);
	var = atoi(word);
}

void get_double(char *str, int &ipos, double &var)
{
	int l;
	bool flag;
	char ch;
	char word[32];
	flag = true;
	l = 0;
	do
	{
		ch = str[ipos];
		if (!(ch == ' ' || ch == '\t' || ch == '\n'))
		{
			word[l] = ch;
			l++;
		}
		else
		{
			if (l)
			{
				flag = false;
				word[l] = '\0';
			}
		}
		ipos++;
	} while (flag);
	var = atof(word);
}

int GetNumberOfPlaces(int &npls)
{
	int i,k,p1,p2,nprof;
	ifstream inf;

	npls=0;
	inf.open("group");
	if(!inf)
	{
		logfile<<"Error in open file "<<"group"<<endl;
		cout<<"Error in open file "<<"group"<<endl;
		return 1;
	}
	inf>>nprof;
	for(i=0;i<nprof;i++)
	{
		inf>>k>>p1>>p2;
		npls+=p2-p1+1;
	}
	inf.close();
	inf.clear();

	return 0;
}

int main()
{
	int i,k,retp;
	ifstream inf;
	ofstream ofp;
	double x,y,z;
	
	int p1,p2,ipls,npls,nprof;
	vector<LineAB> GenAB;
	vector<double> zgel,zvela,zvelb;

	int nGels,nVels,nGens;

	task2d task("Hfi");





	retp=GetNumberOfPlaces(npls);
	if(retp)
	{
		logfile<<"Error in function "<<"GetNumberOfPlaces"<<endl;
		cout<<"Error in function "<<"GetNumberOfPlaces"<<endl;
		return 1;
	}

	nGels=0;
	inf.open("srsclcgsz");
	if(!inf)
	{
		logfile<<"Error in open file "<<"srsclcgsz"<<endl;
		cout<<"Error in open file "<<"srsclcgsz"<<endl;
		return 1;
	}
	for(ipls=0;ipls<npls;ipls++)
	{
		inf>>k;
		nGels+=k;
	}
	inf.close();
	inf.clear();


	nVels=0;
	inf.open("srsclcvsz");
	if(!inf)
	{
		logfile<<"Error in open file "<<"srsclcvsz"<<endl;
		cout<<"Error in open file "<<"srsclcvsz"<<endl;
		return 1;
	}
	for(ipls=0;ipls<npls;ipls++)
	{
		inf>>k;
		nVels+=k;
	}
	inf.close();
	inf.clear();

	nGens=nGels+nVels;

	cout<<"nGels= "<<nGels<<'\n';
	cout<<"nVels= "<<nVels<<'\n';
	cout<<"nGens= "<<nGens<<'\n';

	GenAB.resize(nGens);
	zgel.resize(nGels);
	zvela.resize(nVels);
	zvelb.resize(nVels);

	npls=0;
	inf.open("srsclch");
	if(!inf)
	{
		logfile<<"Error in open file "<<"srsclch"<<endl;
		cout<<"Error in open file "<<"srsclch"<<endl;
		return 1;
	}
	for(i=0;i<nGens;i++)
	{
		inf>>x>>y>>z;
		GenAB[i].A.x=x;
		GenAB[i].A.y=y;
		GenAB[i].A.z=z;
		inf>>x>>y>>z;
		GenAB[i].B.x=x;
		GenAB[i].B.y=y;
		GenAB[i].B.z=z;
	}
	inf.close();
	inf.clear();

	for(i=0;i<nGels;i++){zgel[i]=GenAB[i].A.z;}
	for(i=0;i<nVels;i++)
	{
		zvela[i]=GenAB[nGels+i].A.z;
		zvelb[i]=GenAB[nGels+i].B.z;
		if(zvela[i]>zvelb[i])
		{
			x=zvelb[i];
			zvelb[i]=zvela[i];
			zvela[i]=x;
		}
	}

	retp=task.Read();
	if(retp)
	{
		logfile<<"Function Read() for Hphi returned "<<retp<<'\n';
		return 1;
	}

	task.GenKU(nGels,nVels,zgel,zvela,zvelb);

	logfile.close();
	logfile.clear();

	return 0;
}
