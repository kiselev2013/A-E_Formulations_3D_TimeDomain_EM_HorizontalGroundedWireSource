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
#include "GeoPrepDocSettings.h"

extern GeoPrepDocSettings GPDocSettings;

extern bool CheckStop(void);

extern ofstream logfile;

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


double FrontGetValue(const double& _t,const vector<double> &t,const vector<double> &v,const double& _E0)
{
	int i;
	const double EPS=1e-10;
	double res=0.0;
	if(_t<0.0)res=_E0;
	else if(_t<t[0])res=_E0+(v[0]-_E0)*_t/t[0];
	else{
		i=FindIntervalInDoubleVec(t,_t);
		if(i!=-1){
			res=v[i]+(_t-t[i])*(v[i+1]-v[i])/(t[i+1]-t[i]);
		}
		else{
			res=0.0;
		}
	}
	return res;
}


void CalcFront(double &E0,								//  
			  vector<double> &t, vector<double> &v,		// 
			  vector<double> &Arg, vector<double> &Val,	// 
			  vector<double> &)						// 
{
	int i, j, n;
	n=(int)Arg.size();
	.resize(n);
	[0]=0.0;
	for(i=1;i<n;i++){
		double sum;
		sum=0;
		for(j=i-1;j>0;j--){
			sum+=0.5*(Val[j]+Val[j-1])*(FrontGetValue(Arg[i]-Arg[j],t,v,E0)-FrontGetValue(Arg[i]-Arg[j-1],t,v,E0));
		}
		[i]=sum+Val[i]*(E0-FrontGetValue(Arg[i]-Arg[i-1],t,v,E0));
	}
}

int SummatorE(
	const char *rdir,
	const int& i, // profile number
	const int& j, // profile point number
	const int& k, // parameter number
	const int& mn_sea,
	const bool& WithDDV)
{
	char bufx[1024], bufy[1024], bufz[1024], buf[1024];
	ifstream infx, infy, infz, inf;

	Point3D pA, pB;
	Vector vAB;

	vector<Picket> Ex;
	vector<Picket> Ey;
	vector<Picket> Ez;

	vector<Picket> Edmn;
	vector<Picket> Edmoon;

	vector<Picket> EdmnImp;
	vector<Picket> EdmoonImp;

	vector<Picket> EdmnFrn;
	vector<Picket> EdmoonFrn;

	int q, m, l, ii, jj;

	double Pimp, Ppause;
	int NT2, ImpType;

	inf.open("imp");
	if(!inf)return 1;
	inf>>Pimp>>Ppause>>NT2>>ImpType;
	inf.close();
	inf.clear();

	Pimp*=1000.; Ppause*=1000.;

	if (k==-1)
	{
		sprintf(bufx, "%s\\exall", rdir);
		sprintf(bufy, "%s\\eyall", rdir);
		sprintf(bufz, "%s\\ezall", rdir);
	}
	else
	{
		sprintf(bufx, "%s\\exall%d_%d_%d", rdir, i, j, k+1);
		sprintf(bufy, "%s\\eyall%d_%d_%d", rdir, i, j, k+1);
		sprintf(bufz, "%s\\ezall%d_%d_%d", rdir, i, j, k+1);
	}

	sprintf(buf, "xyzmn");

	infx.open(bufx);
	infy.open(bufy);
	infz.open(bufz);
	inf.open(buf);

	if (!infx||!infy||!infz||!inf)
	{
		if (infx)
			infx.close();
		if (infy)
			infy.close();
		if (infz)
			infz.close();
		if (inf)
			inf.close();
		return 1;
	}

	inf>>q;

	Ex.resize(q*GPDocSettings.line_mrad_coeff);
	Ey.resize(q*GPDocSettings.line_mrad_coeff);
	Ez.resize(q*GPDocSettings.line_mrad_coeff);
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
		if(CheckStop())return 0;

		Picket& Edmn_ii=Edmn[ii];
		Picket& Edmoon_ii=Edmoon[ii];
		inf>>pA.x()>>pA.y()>>pA.z();
		inf>>pB.x()>>pB.y()>>pB.z();
		Edmn_ii.GetPnt().x()=0.5*(pA.getx()+pB.getx());
		Edmn_ii.GetPnt().y()=0.5*(pA.gety()+pB.gety());
		Edmn_ii.GetPnt().z()=0.5*(pA.getz()+pB.getz());
		strcpy(Edmn_ii.leg, Ex[m].leg);
		if(WithDDV)
		{
			Edmoon_ii.GetPnt().x()=0.5*(pA.getx()+pB.getx());
			Edmoon_ii.GetPnt().y()=0.5*(pA.gety()+pB.gety());
			Edmoon_ii.GetPnt().z()=0.5*(pA.getz()+pB.getz());
		}
		strcpy(Edmoon_ii.leg, Ex[m].leg);
		vAB=pB-pA;
		vAB/=(double)GPDocSettings.line_mrad_coeff;
		for (jj=0; jj<GPDocSettings.line_mrad_coeff; jj++)
		{
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
					if (jj>=(GPDocSettings.line_mrad_coeff/2))
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

	double sea;
	vector<double> vsea;
	vsea.resize(q);

	if (mn_sea)
	{
		if (k==-1)
			sprintf(buf, "%s\\dV.000", rdir);
		else
			sprintf(buf, "%s\\dV.000_%d_%d_%d", rdir, i, j, k+1);

		inf.open(buf);
		if(!inf)return 1;
		for (ii=0; ii<q; ii++)
		{
			inf>>vsea[ii];
		}
		inf.close();
		inf.clear();
	}

	bool front;
	int itime,ktime,NtimeEd;
	vector<double> ImpT,ImpV;
	vector<double> VecVal[3],VecRes[3];
	double E0;
	inf.open("impulse");
	front=false;
	if(inf){
		front=true;
		inf>>ktime;
		ImpT.resize(ktime);
		ImpV.resize(ktime);
		for(itime=0;itime<ktime;itime++)
		{
			inf>>ImpT[itime]>>ImpV[itime];
			ImpT[itime]*=1000.0;
		}
		inf.close();
		inf.clear();
		
		if (k==-1)
			sprintf(buf, "%s\\dV.000", rdir);
		else
			sprintf(buf, "%s\\dV.000_%d_%d_%d", rdir, i, j, k+1);
		inf.open(buf);
		EdmnFrn.resize(q);
		for (ii=0; ii<q; ii++)
		{
			if(CheckStop())return 0;

			inf>>E0;
			NtimeEd=(int)Edmn[ii].GetTimes().size();
			EdmnFrn[ii].GetTimes().resize(NtimeEd);
			EdmnFrn[ii].GetEds().resize(NtimeEd);
			for(jj=0;jj<NtimeEd;jj++)
			{
				EdmnFrn[ii].GetTimes()[jj]=Edmn[ii].GetTimes()[jj];
				EdmnFrn[ii].GetEds()[jj]=Edmn[ii].GetEds()[jj];
			}
			vector<double> &lt=EdmnFrn[ii].GetTimes();
			vector<Picket::EdsVal> &lp=EdmnFrn[ii].GetEds();
			VecVal[0].resize(NtimeEd);
			VecVal[1].resize(NtimeEd);
			VecVal[2].resize(NtimeEd);
			VecRes[0].resize(NtimeEd);
			VecRes[1].resize(NtimeEd);
			VecRes[2].resize(NtimeEd);
			for(jj=0;jj<NtimeEd;jj++)
			{
				VecVal[0][jj]=lp[jj].n;
				VecVal[1][jj]=lp[jj].a;
				VecVal[2][jj]=lp[jj].s;
			}
			CalcFront(E0,lt,VecVal[0],ImpT,ImpV,VecRes[0]);
			CalcFront(E0,lt,VecVal[1],ImpT,ImpV,VecRes[1]);
			CalcFront(E0,lt,VecVal[2],ImpT,ImpV,VecRes[2]);
			for(jj=0;jj<NtimeEd;jj++)
			{
				lp[jj].n=FrontGetValue(lt[jj],ImpT,VecRes[0],E0);//VecRes[0][jj];
				lp[jj].a=FrontGetValue(lt[jj],ImpT,VecRes[1],E0);//VecRes[1][jj];
				lp[jj].s=FrontGetValue(lt[jj],ImpT,VecRes[2],E0);//VecRes[2][jj];
			}
		}
		inf.close();
		inf.clear();

		if(WithDDV)
		{
			if (k==-1)
				sprintf(buf, "%s\\ddV.000");
			else
				sprintf(buf, "%s\\ddV.000_%d_%d_%d", rdir, i, j, k+1);
			inf.open(buf);
			EdmoonFrn.resize(q);
			for (ii=0; ii<q; ii++)
			{
				inf>>E0;
				NtimeEd=(int)Edmoon[ii].GetTimes().size();
				EdmoonFrn[ii].GetTimes().resize(NtimeEd);
				EdmoonFrn[ii].GetEds().resize(NtimeEd);
				for(jj=0;jj<NtimeEd;jj++)
				{
					EdmoonFrn[ii].GetTimes()[jj]=Edmoon[ii].GetTimes()[jj];
					EdmoonFrn[ii].GetEds()[jj]=Edmoon[ii].GetEds()[jj];
				}
				vector<double> &lt=EdmoonFrn[ii].GetTimes();
				vector<Picket::EdsVal> &lp=EdmoonFrn[ii].GetEds();
				VecVal[0].resize(NtimeEd);
				VecVal[1].resize(NtimeEd);
				VecVal[2].resize(NtimeEd);
				VecRes[0].resize(NtimeEd);
				VecRes[1].resize(NtimeEd);
				VecRes[2].resize(NtimeEd);
				for (jj=0;jj<NtimeEd;jj++)
				{
					VecVal[0][jj]=lp[jj].n;
					VecVal[1][jj]=lp[jj].a;
					VecVal[2][jj]=lp[jj].s;
				}
				CalcFront(E0,lt,VecVal[0],ImpT,ImpV,VecRes[0]);
				CalcFront(E0,lt,VecVal[1],ImpT,ImpV,VecRes[1]);
				CalcFront(E0,lt,VecVal[2],ImpT,ImpV,VecRes[2]);
				for(jj=0;jj<NtimeEd;jj++)
				{
					lp[jj].n=FrontGetValue(lt[jj],ImpT,VecRes[0],E0);//VecRes[0][jj];
					lp[jj].a=FrontGetValue(lt[jj],ImpT,VecRes[1],E0);//VecRes[1][jj];
					lp[jj].s=FrontGetValue(lt[jj],ImpT,VecRes[2],E0);//VecRes[2][jj];
				}
			}
			inf.close();
			inf.clear();
		}
	}
	inf.clear();

	EdmnImp.resize(Edmn.size());
	EdmoonImp.resize(Edmoon.size());

	ofstream outf;

	if (k==-1)
		sprintf(buf, "%s\\edsall", rdir);
	else
		sprintf(buf, "%s\\edsall%d_%d_%d", rdir, i, j, k+1);

	outf.open(buf);
	for (ii=0; ii<q; ii++)
	{
		if(CheckStop())return 0;
		Picket& Edmn_ii=Edmn[ii];
		EdmnImp[ii]=Edmn_ii;
		if(front)
		{
			Picket& EdmnFrn_ii=EdmnFrn[ii];
			for (jj=0; jj<(int)Edmn_ii.GetTimes().size(); jj++)
				Edmn_ii.GetEds()[jj]+=EdmnFrn_ii.GetEds()[jj];
		}
		if (mn_sea)
		{
			sea=(mn_sea==1)? vsea[0] : vsea[ii];
			for (jj=0; jj<(int)Edmn_ii.GetTimes().size(); jj++)
			{
				Edmn_ii.GetEds()[jj]*=1.0/sea;
			}
		}
		Edmn_ii.WritePicket(outf);
	}
	outf.close();

	if (k==-1)
		sprintf(buf, "%s\\edsall_imp", rdir);
	else
		sprintf(buf, "%s\\edsall%d_%d_%d_imp", rdir, i, j, k+1);

	outf.open(buf);
	for (ii=0; ii<q; ii++)
	{
		if(CheckStop())return 0;
		Picket& Edmn_ii=EdmnImp[ii];
		Edmn_ii.HevToImp2(Pimp, Ppause, NT2, ImpType);
		if(front)
		{
			Picket& EdmnFrn_ii=EdmnFrn[ii];
			for (jj=0; jj<(int)Edmn_ii.GetTimes().size(); jj++)
				Edmn_ii.GetEds()[jj]+=EdmnFrn_ii.GetEds()[jj];
		}
		if (mn_sea)
		{
			sea=(mn_sea==1)? vsea[0] : vsea[ii];
			for (jj=0; jj<(int)Edmn_ii.GetTimes().size(); jj++)
			{
				Edmn_ii.GetEds()[jj]*=1.0/sea;
			}
		}
		Edmn_ii.WritePicket(outf);
	}
	outf.close();

	if(WithDDV){

		if (k==-1)
			sprintf(buf, "%s\\edsallddv", rdir);
		else
			sprintf(buf, "%s\\edsallddv%d_%d_%d", rdir, i, j, k+1);

		outf.open(buf);
		for (ii=0; ii<q; ii++)
		{
			if(CheckStop())return 0;
			Picket& Edmoon_ii=Edmoon[ii];
			EdmoonImp[ii]=Edmoon_ii;
			if(front)
			{
				Picket& EdmoonFrn_ii=EdmoonFrn[ii];
				for (jj=0; jj<(int)Edmoon_ii.GetTimes().size(); jj++)
					Edmoon_ii.GetEds()[jj]+=EdmoonFrn_ii.GetEds()[jj];
			}
			if (mn_sea)
			{
				sea=(mn_sea==1)? vsea[0] : vsea[ii];
				for (jj=0; jj<(int)Edmoon_ii.GetTimes().size(); jj++)
				{
					Edmoon_ii.GetEds()[jj]*=1.0/sea;
				}
			}
			Edmoon_ii.WritePicket(outf);
		}
		outf.close();

		if (k==-1)
			sprintf(buf, "%s\\edsallddv_imp", rdir);
		else
			sprintf(buf, "%s\\edsallddv%d_%d_%d_imp", rdir, i, j, k+1);

		outf.open(buf);
		for (ii=0; ii<q; ii++)
		{
			if(CheckStop())return 0;
			Picket& Edmoon_ii=EdmoonImp[ii];
			Edmoon_ii.HevToImp2(Pimp, Ppause, NT2, ImpType);
			if(front)
			{
				Picket& EdmoonFrn_ii=EdmoonFrn[ii];
				for (jj=0; jj<(int)Edmoon_ii.GetTimes().size(); jj++)
					Edmoon_ii.GetEds()[jj]+=EdmoonFrn_ii.GetEds()[jj];
			}
			if (mn_sea)
			{
				sea=(mn_sea==1)? vsea[0] : vsea[ii];
				for (jj=0; jj<(int)Edmoon_ii.GetTimes().size(); jj++)
				{
					Edmoon_ii.GetEds()[jj]*=1.0/sea;
				}
			}
			Edmoon_ii.WritePicket(outf);
		}
		outf.close();

	}

	return 0;
}

int SumEdsall(char *s1,char *s2,char *res)
{
	ifstream ifa1,ifa2;
	ofstream ofa;
	Picket pkt;
	char buf1[1024],buf2[1024];

	ifa1.clear();
	ifa1.open(s1);
	if(!ifa1)
	{
		logfile<<"File "<<s1<<" not opend"<<'\n';
		return 1;
	}
	ifa2.clear();
	ifa2.open(s2);
	if(!ifa2)
	{
		logfile<<"File "<<s2<<" not opend"<<'\n';
		return 1;
	}
	ofa.open(res);
	ifa1.getline(buf1, 1000, '\n');
	ifa2.getline(buf2, 1000, '\n');
	while(ifa1 && !ifa1.eof() && ifa2 && !ifa2.eof())
	{
		pkt.Clear();
		pkt.ReadPicket(ifa1,buf1);
		pkt.ReadPicket(ifa2,buf2,false);
		pkt.WritePicket(ofa);
	}
	ifa1.close();
	ifa1.clear();
	ifa2.close();
	ifa2.clear();
	ofa.close();
	ofa.clear();

	return 0;
}

int AddFloatFile(char *FileName0,char *FileName1,char *FileName2)
{
	int fh;
	long posb,pose;
	int i,n;
	const int size_f=sizeof(float);
	float f1,f2,fs;
	FILE *fpi1,*fpi2,*fpo;

	_sopen_s( &fh, FileName0, _O_RDONLY, _SH_DENYNO, 0 );

	/* Seek the beginning of the file: */
	posb = _lseek( fh, 0L, SEEK_SET );
	if( posb == -1L )return 1;

	/* Set the end of the file: */
	pose = _lseek( fh, 0L, SEEK_END );
	if( pose == -1L )return 1;

	_close( fh );

	n=(pose-posb)/size_f;

	fpi1=fopen(FileName0,"rb");
	fpi2=fopen(FileName1,"rb");
	fpo=fopen(FileName2,"wb");

	if(!fpi1 || !fpi2)return 1;

	for(i=0;i<n;i++)
	{
		fread(&f1,size_f,1,fpi1);
		fread(&f2,size_f,1,fpi2);
		fs=f1+f2;
		fwrite(&fs,size_f,1,fpo);
	}

	fclose(fpi1);
	fclose(fpi2);
	fclose(fpo);

	return 0;
}

int UnionEdsall(char *res,int n,char *s,...)
{
	int i;
	ifstream ifa;
	ofstream ofa;
	Picket pkt;
	char buf[1024];
	char **ps;

	ps=&s;

	ofa.open(res);

	for(i=0;i<n;i++)
	{
		ifa.open(*ps);
		if(!ifa)
		{
			logfile<<"File "<<*ps<<" not opend"<<'\n';
			ofa.close();
			ofa.clear();
			return 1;
		}
		ifa.getline(buf, 1000, '\n');
		while(ifa && !ifa.eof())
		{
			pkt.Clear();
			pkt.ReadPicket(ifa,buf);
			pkt.WritePicket(ofa);
		}
		ifa.close();
		ifa.clear();

		ps++;
	}

	ofa.close();
	ofa.clear();

	return 0;
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

int DivEdsall(char *s1,char *s2,char *res)
{
	int i,n;
	ifstream ifa1,ifa2;
	ofstream ofa;
	Picket pkt1,pkt2,pktr;
	char buf1[1024],buf2[1024];
	ifa1.clear();
	ifa1.open(s1);
	if(!ifa1)
	{
		logfile<<"File "<<s1<<" not opend"<<'\n';
		return 1;
	}
	ifa2.clear();
	ifa2.open(s2);
	if(!ifa2)
	{
		logfile<<"File "<<s2<<" not opend"<<'\n';
		return 1;
	}
	ofa.open(res);
	ifa1.getline(buf1, 1000, '\n');
	ifa2.getline(buf2, 1000, '\n');
	while(ifa1 && !ifa1.eof() && ifa2 && !ifa2.eof())
	{
		pkt1.Clear();
		pkt1.ReadPicket(ifa1,buf1);
		pktr=pkt1;
		pkt2.Clear();
		pkt2.ReadPicket(ifa2,buf2);
		vector<Picket::EdsVal> &v2=pkt2.GetEds();
		vector<Picket::EdsVal> &vr=pktr.GetEds();
		pktr=pkt1;
		n=(int)vr.size();
		for(i=0;i<n;i++)
		{
			if(fabs(v2[i].n)>1e-20)
			{
				vr[i].n/=v2[i].n;
			}
			if(fabs(v2[i].a)>1e-20)
			{
				vr[i].a/=v2[i].a;
			}
			if(fabs(v2[i].s)>1e-20)
			{
				vr[i].s/=v2[i].s;
			}
		}
		pktr.WritePicket(ofa);
	}
	ifa1.close();
	ifa1.clear();
	ifa2.close();
	ifa2.clear();
	ofa.close();
	ofa.clear();

	return 0;
}

void ClearTextFile(char *s1)
{
	ofstream ofp;
	ofp.clear();
	ofp.open(s1);
	ofp.close();
	ofp.clear();
}
