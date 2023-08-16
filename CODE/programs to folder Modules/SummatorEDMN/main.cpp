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
 *  This file contains code for E-field integration for receivers lines
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
#include<fstream>
#include<vector>
#include<string>
#include<sstream>
#include<iomanip>

using namespace std;

#include "Picket.h"

int NUMBEROFPOINTSPERLINE=20;

int SummatorE(const char *rdir,const int& i,const int& j,const int& k,const int& mn_sea,const bool& WithDDV)
{
	char bufx[1024], bufy[1024], bufz[1024], buf[1024], bufe[1024];
	ifstream infx, infy, infz, inf, infe;
	ofstream outf;

	Point3D pA, pB, pE, pS;
	Vector vAB;

	vector<Picket> Ex;
	vector<Picket> Ey;
	vector<Picket> Ez;

	vector<Picket> Edmn;
	vector<Picket> Edmoon;

	int q, m, l, ii, jj;

	if (k==-1)
	{
		sprintf(bufx, "%s\\exall%d_%d", rdir, i, j);
		sprintf(bufy, "%s\\eyall%d_%d", rdir, i, j);
		sprintf(bufz, "%s\\ezall%d_%d", rdir, i, j);
	}
	else
	{
		sprintf(bufx, "%s\\exall%d_%d_%d", rdir, i, j, k+1);
		sprintf(bufy, "%s\\eyall%d_%d_%d", rdir, i, j, k+1);
		sprintf(bufz, "%s\\ezall%d_%d_%d", rdir, i, j, k+1);
	}

	sprintf(buf, "xyzmn%d_%d",j+1,i);
	sprintf(bufe, "xyzVectorElin%d_%d",j+1,i);

	infx.open(bufx);
	if(!infx)
	{
		cout<<"Error - open file: "<<bufx<<endl;
		return 1;
	}
	infy.open(bufy);
	if(!infy)
	{
		cout<<"Error - open file: "<<bufy<<endl;
		return 1;
	}
	infz.open(bufz);
	if(!infz)
	{
		cout<<"Error - open file: "<<bufz<<endl;
		return 1;
	}
	inf.open(buf);
	if(!inf)
	{
		cout<<"Error - open file: "<<buf<<endl;
		return 1;
	}
	infe.open(bufe);
	if(!infe)
	{
		cout<<"Error - open file: "<<bufe<<endl;
		return 1;
	}

	inf>>q;
	infe>>ii;

	NUMBEROFPOINTSPERLINE=ii/q;

	Ex.resize(q*NUMBEROFPOINTSPERLINE);
	Ey.resize(q*NUMBEROFPOINTSPERLINE);
	Ez.resize(q*NUMBEROFPOINTSPERLINE);
	Edmn.resize(q);
	Edmoon.resize(q);

	infx.getline(bufx, 1000, '\n');
	for (ii=0; ii<(int)Ex.size(); ii++)
	{
		Ex[ii].ReadPicket(infx, bufx);
	}
	infy.getline(bufy, 1000, '\n');
	for (ii=0; ii<(int)Ey.size(); ii++)
	{
		Ey[ii].ReadPicket(infy, bufy);
	}
	infz.getline(bufz, 1000, '\n');
	for (ii=0; ii<(int)Ez.size(); ii++)
	{
		Ez[ii].ReadPicket(infz, bufz);
	}

	infx.close();
	infy.close();
	infz.close();


	m=0;
	for (ii=0; ii<q; ii++)
	{
		Picket& Edmn_ii=Edmn[ii];
		Picket& Edmoon_ii=Edmoon[ii];
		inf>>pA.x()>>pA.y()>>pA.z();
		inf>>pB.x()>>pB.y()>>pB.z();
		Edmn_ii.GetPnt().x()=0.5*(pA.getx()+pB.getx());
		Edmn_ii.GetPnt().y()=0.5*(pA.gety()+pB.gety());
		Edmn_ii.GetPnt().z()=0.5*(pA.getz()+pB.getz());
		strcpy(Edmn_ii.leg,"         t, ms                         dVp, mV                         dVs, mV                         dV, mV");
		if(WithDDV)
		{
			Edmoon_ii.GetPnt().x()=0.5*(pA.getx()+pB.getx());
			Edmoon_ii.GetPnt().y()=0.5*(pA.gety()+pB.gety());
			Edmoon_ii.GetPnt().z()=0.5*(pA.getz()+pB.getz());
		}
		strcpy(Edmoon_ii.leg,"         t, ms                         ddVp, mV                         ddVs, mV                         ddV, mV");
		for (jj=0; jj<NUMBEROFPOINTSPERLINE; jj++)
		{
			infe>>pE.x()>>pE.y()>>pE.z();

			if(jj)
			{
				vAB=pE-pS;
			}
			else
			{
				vAB=pE-pA;
				vAB.x()*=2.0;
				vAB.y()*=2.0;
				vAB.z()*=2.0;
			}

			pS=pE;

			const Picket& Ex_m=Ex[m];
			const Picket& Ey_m=Ey[m];
			const Picket& Ez_m=Ez[m];
			for (l=0; l<(int)Ex_m.GetTimes().size(); l++)
			{
				const Picket::EdsVal &vEx=Ex_m.GetEds()[l];
				const Picket::EdsVal &vEy=Ey_m.GetEds()[l];
				const Picket::EdsVal &vEz=Ez_m.GetEds()[l];
				Picket::EdsVal E(vEx.n*vAB.getx()+vEy.n*vAB.gety()+vEz.n*vAB.getz(),
								 vEx.a*vAB.getx()+vEy.a*vAB.gety()+vEz.a*vAB.getz(),
								 vEx.s*vAB.getx()+vEy.s*vAB.gety()+vEz.s*vAB.getz());
				if (jj==0)
				{
					Edmn_ii.GetTimes().push_back(Ex_m.GetTimes()[l]);
					Edmn_ii.GetEds().push_back(E);
					if(WithDDV)
					{
						Edmoon_ii.GetTimes().push_back(Ex_m.GetTimes()[l]);
						Edmoon_ii.GetEds().push_back(E);
					}
				}
				else
				{
					Edmn_ii.GetEds()[l]+=E;
					if (jj>=(NUMBEROFPOINTSPERLINE/2))
						E*=-1.;
					if(WithDDV)
					{
						Edmoon_ii.GetEds()[l]+=E;
					}
				}
			}
			m++;
		}
	}

	inf.close();
	inf.clear();

	infe.close();
	infe.clear();


	if (k==-1)
		sprintf(buf, "%s\\edmnall%d_%d", rdir, i, j);
	else
		sprintf(buf, "%s\\edmnall%d_%d_%d", rdir, i, j, k+1);

	outf.open(buf);
	outf<<scientific<<setprecision(16);
	for (ii=0; ii<q; ii++)
	{
		Picket& Edmn_ii=Edmn[ii];
		Edmn_ii.WritePicket(outf);
	}
	outf.close();
	outf.clear();



	if(WithDDV)
	{
		if (k==-1)
			sprintf(buf, "%s\\edmoonall%d_%d", rdir, i, j);
		else
			sprintf(buf, "%s\\edmoonall%d_%d_%d", rdir, i, j, k+1);

		outf.open(buf);
		outf<<scientific<<setprecision(16);
		for (ii=0; ii<q; ii++)
		{
			Picket& Edmoon_ii=Edmoon[ii];
			Edmoon_ii.WritePicket(outf);
		}
		outf.close();
		outf.clear();


	}

	return 0;
}

int main(int argc,char **argv)
{
	ifstream inf;
	int ipr,jpr,ret;
	bool WithDDV;
	int mn_sea;
	int nLin,nRec;

	if(argc<4)
	{
		cout<<"Error: need 3 parameters "<<endl;
		cout<<"1 - flag WithDDV"<<endl;
		cout<<"2 - number of profile"<<endl;
		cout<<"3 - number of place by profile"<<endl;
		return 1;
	}

	mn_sea=0;
	WithDDV=atoi(argv[1]);
	jpr=atoi(argv[2]);
	ipr=atoi(argv[3]);

	inf.open("xyzmn");
	if(!inf)
	{
		cout<<"Error - open file: "<<"xyzmn"<<endl;
		return 1;
	}
	inf>>nLin;
	inf.close();
	inf.clear();

	if(nLin)
	{
		ret=SummatorE(".",jpr,ipr,-1,mn_sea,WithDDV);
		if(ret)
		{
			return ret;
		}
	}

	return 0;
}
