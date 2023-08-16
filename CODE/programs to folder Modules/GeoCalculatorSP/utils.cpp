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
 *  This file contains code for utility functions
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 January 16, 2023
*/

#include "stdafx.h"
#include "Picket.h"
#include "utils.h"

extern bool CheckStop(void);

extern ofstream logfile;


int FindIntervalInDoubleVec(const vector<double> &vec,const double &elem)
{
	int i,j,k;
	i=0;
	k=(int)vec.size()-1;
	if(k<0 || (elem<vec[0] && fabs(elem-vec[0])>1e-6) || (elem>vec[k]*1.0001 && fabs(elem-vec[k])>1e-6))return -1;
	j=(i+k)/2;
	do{		
		if(elem>vec[j])i=j;
		else k=j;
		j=(i+k)/2;
	}while(j!=i);
	return i;
}

int SimpleFindInDoubleVec(const vector<double> &vec,const double &elem)
{
	int i,j,k;
	k=(int)vec.size();
	j=-1;
	for(i=0;i<k;i++)
	{
		if(fabs(elem)<1e-12)
		{
			if(fabs(vec[i])<1e-12){j=i;}
			break;
		}
		if(fabs(1.0-vec[i]/elem)<1e-6)
		{
			j=i;
			break;
		}
	}
	return j;
}

int ConCatEdsall(char *res,char *add)
{
	ifstream ifa;
	ofstream ofa;
	Picket pkt;
	char buf[1024];
	ifa.open(add);
	if(!ifa)
	{
		logfile<<"File "<<add<<" not opend"<<'\n';
		return 1;
	}
	ofa.open(res,ios::app);
	ifa.getline(buf, 1000, '\n');
	while(ifa && !ifa.eof())
	{
		pkt.Clear();
		pkt.ReadPicket(ifa,buf);
		pkt.WritePicket(ofa);
	}
	ifa.close();
	ifa.clear();
	ofa.close();
	ofa.clear();
	return 0;
}

double GetValueFromParamola(double t1,double t2,double t3,double v1,double v2,double v3,double t)
{
	double dt1,dt2,dt3;

	dt1=t3-t1;
	dt2=t2-t1;
	dt3=t3-t2;

	return v1*((t-t2)*(t-t3))/(dt2*dt1)-v2*((t-t1)*(t-t3))/(dt2*dt3)+v3*((t-t1)*(t-t2))/(dt1*dt3);
}

double GetDiffFromParabola(double arg_i2,double arg_i1,double arg_i,double val_i2,double val_i1,double val_i,double tk)
{
	double dt,dt0,dt1,mult,mult1,mult2;
	double cur_cff=1e-6;
	dt=arg_i-arg_i2;
	dt0=arg_i-arg_i1;
	dt1=arg_i1-arg_i2;
	mult=((tk-arg_i1)+(tk-arg_i2))/(dt*dt0);
	mult1=-((tk-arg_i)+(tk-arg_i2))/(dt1*dt0);
	mult2=((tk-arg_i)+(tk-arg_i1))/(dt*dt1);
	return mult*val_i+mult1*val_i1+mult2*val_i2;
}

double GetValueLin(const int i,const double& _t,const vector<double> &t,const vector<double> &v)
{
	return v[i]+(_t-t[i])*(v[i+1]-v[i])/(t[i+1]-t[i]);
}

float GetValueLin(const int i,const float& _t,const vector<float> &t,const vector<float> &v)
{
	return v[i]+(_t-t[i])*(v[i+1]-v[i])/(t[i+1]-t[i]);
}

double GetValue(const int i,const double& _t,const vector<double> &t,const vector<double> &v)
{
	int n=(int)t.size();
	double res=0.0;

	if(i-1>=0 && i+1<n)
	{
		res=GetValueFromParamola(t[i-1],t[i],t[i+1],v[i-1],v[i],v[i+1],_t);
	}
	else if(i+2<n)
	{
		res=GetValueFromParamola(t[i],t[i+1],t[i+2],v[i],v[i+1],v[i+2],_t);
	}
	else
	{
		res=GetValueLin(i,_t,t,v);
	}
	
	return res;
}

float GetValue(const int i,const float& _t,const vector<float> &t,const vector<float> &v)
{
	int n=(int)t.size();
	float res=0.0;

	if(i-1>=0 && i+1<n)
	{
		res=float(GetValueFromParamola(t[i-1],t[i],t[i+1],v[i-1],v[i],v[i+1],_t));
	}
	else if(i+2<n)
	{
		res=float(GetValueFromParamola(t[i],t[i+1],t[i+2],v[i],v[i+1],v[i+2],_t));
	}
	else
	{
		res=GetValueLin(i,_t,t,v);
	}
	
	return res;
}

void WriteMessage(char *str)
{
	cout<<str;
	logfile<<str;
}

void WriteMessage(stringstream &str)
{
	cout<<str.str();
	logfile<<str.str();
}

bool IsFileExist(char *fname)
{
	bool flag;
	ifstream inf;

	flag=false;
	inf.open(fname);
	if(inf)
	{
		flag=true;
		inf.close();
	}
	inf.clear();

	return flag;
}

void ColseProgramm(int err_code)
{
	stringbuffer<<"Closing programm with err_code "<<err_code<<endl;
	exit(err_code);
}

void StopIfErrorReturn(int err_code,char *FuncName)
{
	if(err_code)
	{
		stringbuffer<<"Function "<<FuncName<<" returned "<<err_code<<endl;
		WriteMessage(stringbuffer);
		ColseProgramm(err_code);
	}
}

int CreateProcessForEXENoWait(char *cmdline, char *workdir)
{
	int retp;
	STARTUPINFOA si;
	PROCESS_INFORMATION pi;
	ZeroMemory(&si, sizeof(si));
	si.cb = sizeof(si);
	ZeroMemory(&pi, sizeof(pi));
	cout<<"Start "<<cmdline;
	if(workdir){cout<<" in "<<workdir;}
	cout<<endl;
	retp=CreateProcessA(NULL,  // No module name (use command line). 
		(LPSTR)(const char*)cmdline,// Command line. 
		NULL,				// Process handle not inheritable. 
		NULL,				// Thread handle not inheritable. 
		FALSE,				// Set handle inheritance to FALSE. 
		CREATE_NO_WINDOW,	// No creation flags. 
		NULL,				// Use parent's environment block. 
		workdir,				// Use parent's starting directory. 
		&si,				// Pointer to STARTUPINFO structure.
		&pi);				// Pointer to PROCESS_INFORMATION structure.
	return 0;
}
