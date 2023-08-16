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
 *  This file contains code for read-write recevires data functions
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 January 16, 2023
*/

#include "stdafx.h"
#include "picket.h"
#include "openclose.h"
#include "utils.h"

void Picket::Clear()
{
	leg[0]='\0';
	t.clear();
	v[0].clear();
	v[1].clear();
	v[2].clear();
}

void Picket::Allocate(int n)
{
	int i;
	t.resize(n);
	v[0].resize(n);
	v[1].resize(n);
	v[2].resize(n);
	for(i=0;i<n;i++){v[0][i]=v[1][i]=v[2][i]=0.0;}
}

void Picket::AddPicket(Picket &pkt,double cff)
{
	int i,n;
	if(leg[0]=='\0')
	{
		iPx=pkt.iPx;
		iPy=pkt.iPy;
		Px=pkt.Px;
		Py=pkt.Py;
		Pz=pkt.Pz;
		strcpy(leg,pkt.leg);
	}
	n=(int)t.size();
	if(n!=(int)pkt.t.size())
	{
		CloseProgrammWithError(100);
	}
	for(i=0;i<n;i++)
	{
		v[0][i]+=cff*pkt.v[0][i];
		v[1][i]+=cff*pkt.v[1][i];
		v[2][i]+=cff*pkt.v[2][i];
	}
}

void Picket::SummValues(vector<double> &v1,vector<double> &v2)
{
	int i,n;
	n=(int)t.size();
	if(n!=(int)v1.size() || n!=(int)v2.size())
	{
		CloseProgrammWithError(100);
	}
	for(i=0;i<n;i++)
	{
		v[0][i]=v1[i];
		v[1][i]=v2[i];
		v[2][i]=v1[i]+v2[i];
	}
}

void Picket::CopyPicketInfo(Picket &pkt)
{
	iPx=pkt.iPx;
	iPy=pkt.iPy;
	Px=pkt.Px;
	Py=pkt.Py;
	Pz=pkt.Pz;
}

void Picket::SetLegend(char *pleg)
{
	strcpy(leg,pleg);
}

void ReadPickets(int nr,int nt,vector<Picket> &EdsAll,char *fname,bool fmem)
{
	int tmpi,retv;
	double tmpd;
	int ir,it;
	ifstream inf;
	retv=OpenFile(inf,fname,true);
	if(fmem)
	{
		EdsAll.clear();
		EdsAll.resize(nr);
	}
	for(ir=0;ir<nr;ir++)
	{
		Picket &pkt=EdsAll[ir];
		vector<double> &t=pkt.t;
		vector<double> &vx=pkt.v[0];
		vector<double> &vy=pkt.v[1];
		vector<double> &vz=pkt.v[2];
		if(fmem)
		{
			t.resize(nt);
			vx.resize(nt);
			vy.resize(nt);
			vz.resize(nt);
		}
		inf>>pkt.iPx>>pkt.iPy;
		inf>>pkt.Px>>pkt.Py>>pkt.Pz;
		inf>>tmpi;
		inf>>tmpi>>tmpi>>tmpi>>tmpi;
		inf>>tmpi;
		inf.ignore(1000,'\n');
		inf.getline(pkt.leg,255,'\n');
		for(it=0;it<nt;it++){inf>>t[it]>>vx[it]>>vy[it]>>vz[it];}
	}
	CloseFile(inf,fname,true);
}

void WritePickets(int nr,int nt,vector<Picket> &EdsAll,char *fname)
{
	int ir,it;
	ofstream ofp;
	ofp.open(fname);
	ofp<<scientific<<setprecision(16);
	for(ir=0;ir<nr;ir++)
	{
		Picket &pkt=EdsAll[ir];
		vector<double> &t=pkt.t;
		vector<double> &vx=pkt.v[0];
		vector<double> &vy=pkt.v[1];
		vector<double> &vz=pkt.v[2];
		ofp<<pkt.iPx<<' '<<pkt.iPy<<'\n';
		ofp<<pkt.Px<<' '<<pkt.Py<<' '<<pkt.Pz<<'\n';
		ofp<<1<<'\n';
		ofp<<1<<' '<<1<<' '<<1<<' '<<1<<'\n';
		ofp<<1<<'\n';
		ofp<<pkt.leg<<'\n';
		for(it=0;it<nt;it++){ofp<<t[it]<<' '<<vx[it]<<' '<<vy[it]<<' '<<vz[it]<<'\n';}
	}
	ofp.close();
	ofp.clear();
}

void WritePicket(int nr,int nt,vector<Picket> &EdsAll,char *fname,int ind)
{
	int ir,it;
	ofstream ofp;
	ofp.open(fname);
	ofp<<scientific<<setprecision(16);
	for(ir=0;ir<nr;ir++)
	{
		if(ind==ir)
		{
			Picket &pkt=EdsAll[ir];
			vector<double> &t=pkt.t;
			vector<double> &vx=pkt.v[0];
			vector<double> &vy=pkt.v[1];
			vector<double> &vz=pkt.v[2];
			ofp<<pkt.iPx<<' '<<pkt.iPy<<'\n';
			ofp<<pkt.Px<<' '<<pkt.Py<<' '<<pkt.Pz<<'\n';
			ofp<<1<<'\n';
			ofp<<1<<' '<<1<<' '<<1<<' '<<1<<'\n';
			ofp<<1<<'\n';
			ofp<<pkt.leg<<'\n';
			for(it=0;it<nt;it++){ofp<<t[it]<<' '<<vx[it]<<' '<<vy[it]<<' '<<vz[it]<<'\n';}
			break;
		}
	}
	ofp.close();
	ofp.clear();
}

void WritePickets(int from,int to,int nt,vector<Picket> &EdsAll,char *fname)
{
	int ir,it;
	ofstream ofp;
	ofp.open(fname);
	ofp<<scientific<<setprecision(16);
	for(ir=from;ir<to;ir++)
	{
		Picket &pkt=EdsAll[ir];
		vector<double> &t=pkt.t;
		vector<double> &vx=pkt.v[0];
		vector<double> &vy=pkt.v[1];
		vector<double> &vz=pkt.v[2];
		ofp<<pkt.iPx<<' '<<pkt.iPy<<'\n';
		ofp<<pkt.Px<<' '<<pkt.Py<<' '<<pkt.Pz<<'\n';
		ofp<<1<<'\n';
		ofp<<1<<' '<<1<<' '<<1<<' '<<1<<'\n';
		ofp<<1<<'\n';
		ofp<<pkt.leg<<'\n';
		for(it=0;it<nt;it++){ofp<<t[it]<<' '<<vx[it]<<' '<<vy[it]<<' '<<vz[it]<<'\n';}
	}
	ofp.close();
	ofp.clear();
}
