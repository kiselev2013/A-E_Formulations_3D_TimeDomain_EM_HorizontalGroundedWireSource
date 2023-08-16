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
 *  This file contains the code for correcting time mesh with piecewise constant steps
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 January 16, 2023
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<vector>
#include<string>

using namespace std;

int FindIntervalInDoubleVec(const vector<double> &vec,const double &elem)
{
	int i,j,k;
	i=0;
	k=(int)vec.size()-1;
	if(!k || elem<vec[0] || elem>vec[k])return -1;
	j=(i+k)/2;
	do{		
		if(elem>vec[j])i=j;
		else k=j;
		j=(i+k)/2;
	}while(j!=i);
	return i;
}

double GetValueFromParabola(double t1,double t2,double t3,double v1,double v2,double v3,double t)
{
	double dt1,dt2,dt3;

	dt1=t3-t1;
	dt2=t2-t1;
	dt3=t3-t2;

	return ((t-t2)*(t-t3))/(dt2*dt1)*v1-((t-t1)*(t-t3))/(dt2*dt3)*v2+((t-t1)*(t-t2))/(dt1*dt3)*v3;
}

double GetValueFromLinear(double t1,double t2,double v1,double v2,double t)
{
	return ((t2-t)*v1+(t-t1)*v2)/(t2-t1);
}

double Interpolate(double t,const vector<double> &vt,const vector<double> &vv)
{
	int it,let;
	double res;

	let=(int)vt.size()-1;

	if(t<vt[0])
	{
		res=vv[0];
	}
	else if(t>vt[let])
	{
		res=vv[let];
	}
	else
	{
		it=FindIntervalInDoubleVec(vt,t);
			res=GetValueFromLinear(vt[it],vt[it+1],vv[it],vv[it+1],t);
	}

	return res;
}

void GenerateDecInfo(int ntime,double *time,int &nDec,vector<int> &DecIg)
{
	int i;
	double htpre,htcur;
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
}

void SaveTimeMesh(int ntime,vector<double> &time)
{
	int i,j;
	ifstream inf;
	ofstream ofp;
	double recrad;
	int ntold,it;
	vector<double> tim,val,value;

	inf.open("deltafunction");
	inf>>recrad;
	inf>>ntold;
	if(!inf)
	{
		cout<<"Error in open file"<<"deltafunction"<<endl;
		exit(1);
	}
	tim.resize(ntold);
	val.resize(ntold);
	for(it=0;it<ntold;it++){inf>>tim[it]>>val[it];}
	inf.close();
	inf.clear();
	
	rename("currentfunction","!!currentfunction");
	rename("infite.0","!!infite.0");
	rename("deltafunction","!!deltafunction");

	value.resize(ntime);
	for(it=0;it<ntime;it++)
	{
		value[it]=Interpolate(time[it],tim,val);
	}

	ofp.open("deltafunction");
	ofp<<scientific<<setprecision(14);
	ofp<<recrad<<'\n';
	ofp<<ntime<<'\n';
	for(it=0;it<ntime;it++)
	{
		ofp<<time[it]<<' '<<value[it]<<'\n';
	}
	ofp.close();
	ofp.clear();

	ofp.open("currentfunction");
	ofp<<scientific<<setprecision(14);
	for(it=0;it<ntime;it++)
	{
		ofp<<time[it]<<' '<<value[it]<<'\n';
	}
	ofp.close();
	ofp.clear();

	ofp.open("infite.0");
	ofp<<scientific<<setprecision(7);
	ofp<<"    ktime=  "<<ntime<<";  ntstop=  "<<ntime<<";   kiter=    1;   ntime=    1;   niter=    1;"<<'\n';
    ofp<<"   kprogm=    1;  kpropi=    1;  kpropt=   -1;  kitrel=  250;              ;"<<'\n';
    ofp<<"       u0=   0.00000e+000"<<'\n';
	ofp<<" T I M E :"<<'\n';
	for(i=0;i<ntime;i++)
	{
		for(j=0;j<5;j++)
		{
			if(!(i+j<ntime)){break;}
			ofp<<" "<<time[i+j]<<';';
		}
		ofp<<'\n';
		i+=4;
	}
	ofp.close();
	ofp.clear();
}

int main(int argc, char **argv)
{
	setlocale(LC_ALL, "English");

	ifstream inf;

	int itime,ntime;
	vector<double> time,TimeCor;

	int iDec,nDec;
	vector<int> DecIg;

	int i;
	double tcff,tcur,hpre;

	inf.open("tcff");
	if(!inf)
	{
		cout<<"Error in open file tcff"<<endl;
		exit(1);
	}
	inf>>tcff;
	inf.close();
	inf.clear();


	inf.open("infite.0");
	if(!inf)
	{
		cout<<"Error in open file infite.0"<<endl;
		exit(1);
	}
	inf.ignore(1000, '=');
	inf >> ntime;
	time.resize(ntime);
	inf.ignore(1000, '\n');
	inf.ignore(1000, '\n');
	inf.ignore(1000, '\n');
	inf.ignore(1000, '\n');
	for(itime=0;itime<ntime;itime++)
	{
		if(!inf.good())
		{
			cout<<"infite.0 not good"<<endl;
			exit(1);
		}
		inf>>time[itime];
		inf.ignore(1000, ';');
	}
	inf.close();
	inf.clear();

	GenerateDecInfo(ntime,&(time.front()),nDec,DecIg);

	TimeCor.clear();
	for(iDec=0;iDec<nDec;iDec++)
	{
		if((DecIg[iDec+1]-DecIg[iDec])>1 && (iDec==0 || (DecIg[iDec+1]-DecIg[iDec])*tcff>(DecIg[iDec]-DecIg[iDec-1])))
		{
			for(itime=DecIg[iDec];itime<DecIg[iDec+1];itime++)
			{
				TimeCor.push_back(time[itime]);
			}
		}
		else
		{
			if(iDec==0 && (DecIg[iDec+1]-DecIg[iDec])<2)
			{
				cout<<"Error in first decade"<<endl;
				return 1;
			}
			else
			{
				hpre=(time[DecIg[iDec]-1]-time[DecIg[iDec-1]])/(DecIg[iDec]-DecIg[iDec-1]-1);
				
				i=1;
				tcur=time[DecIg[iDec]-1]+hpre;
				while(tcur-0.5*hpre<time[DecIg[iDec+1]-1])
				{
					TimeCor.push_back(tcur);
					i++;
					tcur=time[DecIg[iDec]-1]+i*hpre;
				}
			}
		}
	}

	SaveTimeMesh((int)TimeCor.size(),TimeCor);

	return 0;
}
