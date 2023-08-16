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
 *  This file contains the code for getting time mesh with piecewise constant steps
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 January 16, 2023
*/

#include "stdafx.h"
#include "TimeStep.h"

ofstream logfile;

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

int main()
{
	int ret;

	ifstream inf;
	ofstream ofp;

	int it,sznew,ntold;
	vector<double> tim,val,val_new;
	double tmpd,recrad;

	logfile.open("logfile.TimeStep");

	TimeStep ts;

	double tcff; 

	inf.open("tcff");
	if(!inf)
	{
		cout<<"Error in open file tcff"<<endl;
		logfile<<"Error in open file tcff"<<endl;
		return 1;
	}
	inf>>tcff;
	inf.close();
	inf.clear();

	ret=ts.Read_infite0("infite.0");
	if(ret)
	{
		cout<<"Error: Function "<<"Read_infite0"<<" returned "<<ret<<endl;
		logfile<<"Error: Function "<<"Read_infite0"<<" returned "<<ret<<endl;
		return 1;
	}
	ts.GenRegMesh(tcff);

	rename("infite.0","!infite.0");

	ts.WriteTimes("infite.0");

	sznew=(int)ts.gg_times.size();
	val_new.resize(sznew);

	inf.open("currentfunction");
	if(!inf)
	{
		cout<<"Error in open file"<<"currentfunction"<<endl;
		logfile<<"Error in open file"<<"currentfunction"<<endl;
		return 1;
	}
	tim.clear();
	val.clear();
	while(!inf.eof())
	{
		inf>>tmpd;
		if(inf.eof() || !inf.good()){break;}
		tim.push_back(tmpd);
		inf>>tmpd;
		val.push_back(tmpd);
	}
	inf.close();
	inf.clear();

	for(it=0;it<sznew;it++)
	{
		vector<double> &tim_new=ts.gg_times;
		val_new[it]=Interpolate(tim_new[it],tim,val);
	}

	rename("currentfunction","!currentfunction");

	ofp.open("currentfunction");
	ofp<<scientific<<setprecision(14);
	for(it=0;it<sznew;it++)
	{
		vector<double> &tim_new=ts.gg_times;
		ofp<<tim_new[it]<<' '<<val_new[it]<<'\n';
	}
	ofp.close();
	ofp.clear();

	inf.open("deltafunction");
	inf>>recrad;
	inf>>ntold;
	if(!inf)
	{
		cout<<"Error in open file"<<"deltafunction"<<endl;
		logfile<<"Error in open file"<<"deltafunction"<<endl;
		return 1;
	}
	tim.resize(ntold);
	val.resize(ntold);
	for(it=0;it<ntold;it++){inf>>tim[it]>>val[it];}
	inf.close();
	inf.clear();

	for(it=0;it<sznew;it++)
	{
		vector<double> &tim_new=ts.gg_times;
		val_new[it]=Interpolate(tim_new[it],tim,val);
	}

	rename("deltafunction","!deltafunction");

	ofp.open("deltafunction");
	ofp<<scientific<<setprecision(14);
	ofp<<recrad<<'\n';
	ofp<<sznew<<'\n';
	for(it=0;it<sznew;it++)
	{
		vector<double> &tim_new=ts.gg_times;
		ofp<<tim_new[it]<<' '<<val_new[it]<<'\n';
	}
	ofp.close();
	ofp.clear();

	logfile.close();
	logfile.clear();

	return 0;
}
