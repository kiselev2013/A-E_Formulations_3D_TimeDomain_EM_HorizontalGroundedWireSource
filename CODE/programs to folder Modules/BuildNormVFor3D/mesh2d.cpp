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
 *  This file contains the code for working with a 2D mesh and solution
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 January 16, 2023
*/

#include"stdafx.h"
#include"mesh2d.h"
#include"retcode.h"

#define CHECKSTRING0(s, n) if (int(sizeof(s)/sizeof(s[0])-1-n)<0) return RETCODE_OUTOFRANGE;

extern ofstream logfile;

ifstream& operator>(ifstream& inf, PointRZ& p)
{
	inf > p.r > p.z;
	return inf;
}

ifstream& operator>(ifstream& inf, Rect& r)
{
	inf > r.nodes[2] > r.nodes[3] > r.nodes[0] > r.nodes[1] > r.nodes[4] > r.rtype;
	return inf;
}

MeshForVP::MeshForVP(int mOptions):taskpath("."),meshoptions(mOptions)
{
	pnt=NULL;
	q_all=NULL;
	qv=NULL;
	rect=NULL;
	l1=NULL;
	ResCalcPointArray=NULL;
	ResCalcRectArray=NULL;
	ResCalcPointNum=0;
	ResCalcRectNum=0;
	ntimes=0;
	Rbak=0.0;
	reg=NULL;
	rm=NULL;
	zm=NULL;
}
MeshForVP::~MeshForVP()
{
	if (pnt) delete [] pnt;
	if (q_all) delete [] q_all;
	if (qv) delete [] qv;
	if (rect) delete [] rect;
	if (l1) delete [] l1;
	if (ResCalcPointArray) delete [] ResCalcPointArray;
	if (ResCalcRectArray) delete [] ResCalcRectArray;
	if (reg) delete [] reg;
	if (rm) delete [] rm;
	if (zm) delete [] zm;
}

int MeshForVP::Read(const char* path, bool read_times)
{
	int i, retc, qrptotal=0;
	string buf;
	ifstream inf, rzf, nvtrf, nvkatf, l1f, curf, abf;

	taskpath=path;

	strcpy__m(buf, path);
	inf.open(strcat__m(buf, "\\inf2tr.dat"));
	if (!inf) return RETCODE_NOFILE;
	inf.ignore(1000, '\n');
	inf.ignore(1000, '='); inf>>kpnt;
	inf.ignore(1000, '='); inf>>krect;
	inf.ignore(1000, '='); inf>>kt1;
	inf.close();
	inf.clear();

	strcpy__m(buf, path);
	rzf.open(strcat__m(buf, "\\rz.dat"), ios::binary);
	if (!rzf) return RETCODE_NOFILE;
	pnt=new PointRZ[kpnt];
	if (!pnt) return RETCODE_NOMEM;
	for (i=0; i<kpnt; i++)
		rzf > pnt[i];
	rzf.close();
	rzf.clear();

	q_all=new double[kpnt];
	if (q_all==NULL) return RETCODE_NOMEM;

	strcpy__m(buf, path);
	nvtrf.open(strcat__m(buf, "\\nvtr.dat"), ios::binary);
	if (!nvtrf) return RETCODE_NOFILE;
	rect=new Rect[krect];
	if (!rect) return RETCODE_NOMEM;
	for (i=0; i<krect; i++)
		nvtrf > rect[i];
	nvtrf.close();
	nvtrf.clear();

	strcpy__m(buf, path);
	nvkatf.open(strcat__m(buf, "\\nvkat2d.dat"), ios::binary);
	if (!nvkatf) return RETCODE_NOFILE;
	for (i=0; i<krect; i++)
		nvkatf > rect[i].mtr;
	nvkatf.close();
	nvkatf.clear();

	strcpy__m(buf, path);
	l1f.open(strcat__m(buf, "\\l1.dat"), ios::binary);
	if (!l1f) return RETCODE_NOFILE;
	l1=new int[kt1];
	if (!l1) return RETCODE_NOMEM;
	for (i=0; i<kt1; i++)
		l1f > l1[i];
	l1f.close();
	l1f.clear();

	nc=kpnt;

	qv=new double[nc];
	if (qv==NULL) return RETCODE_NOMEM;

	strcpy__m(buf, path);
	if ((retc=ReadMtrCatalog(strcat__m(buf, "\\sigma"), sigma))!=0)
		return retc;


	ifstream infreg;

	strcpy__m(buf, path);
	infreg.open(strcat__m(buf, "\\rz.txt"));
	if(!infreg) return RETCODE_NOFILE;
	infreg>>qr>>qz;
	infreg.close();
	infreg.clear();

	nreg=krect;

	if(!(reg=new int[nreg])) return RETCODE_NOMEM;
	if(!(rm=new double[qr])) return RETCODE_NOMEM;
	if(!(zm=new double[qz])) return RETCODE_NOMEM;

	strcpy__m(buf, path);
	infreg.open(strcat__m(buf, "\\r.dat"), ios::binary);
	if(!infreg) return RETCODE_NOFILE;
	for(i=0;i<qr;i++){infreg>rm[i];}
	infreg.close();
	infreg.clear();

	strcpy__m(buf, path);
	infreg.open(strcat__m(buf, "\\z.dat"), ios::binary);
	if(!infreg) return RETCODE_NOFILE;
	for(i=0;i<qz;i++){infreg>zm[i];}
	infreg.close();
	infreg.clear();

	for(i=0;i<nreg;i++){reg[i]=i+1;}

	return RETCODE_OK;
}

int MeshForVP::ReadRecievers(const char* path, int pntVsize, const PointXYZ* pntV)
{
	int i, j, retc, qrptotal=0;
	string buf;
	ifstream inf;

#ifdef USEREGINFO
	ResCalcElInfo *ResCalcElArray;
#endif


	ResCalcPointNum=0;
	ResCalcRectNum=0;

	if(ResCalcPointArray){delete [] ResCalcPointArray;ResCalcPointArray=NULL;}
	if(ResCalcRectArray){delete [] ResCalcRectArray; ResCalcRectArray=NULL;}

	ListOfValues<PointForResultCalculation> PointList;


	boundsVA.beginP=PointList.GetListLength();
	for (i=0; i<pntVsize; i++)
	{
		const PointXYZ& p=pntV[i];
		PointForResultCalculation Pnti;
		Pnti.p=p; 
			Pnti.p.x-=A.x;
			Pnti.p.y-=A.y;
		Pnti.r=sqrt(Pnti.p.x*Pnti.p.x+Pnti.p.y*Pnti.p.y);
		if (Pnti.r<1e-3) Pnti.r=1e-3;
		Pnti.sinfi=Pnti.p.y/Pnti.r;
		Pnti.cosfi=Pnti.p.x/Pnti.r;
		if (!PointList.AddToList(Pnti)) return RETCODE_NOMEM;
	}
	boundsVA.endP=PointList.GetListLength();

	if(meshoptions&moForPED)
	{
		boundsVB.beginP=PointList.GetListLength();
		for(j=0;j<NvB;j++)
		{
			PointXYZ &pEL=vB[j];
			for (i=0; i<pntVsize; i++)
			{
				const PointXYZ& p=pntV[i];
				PointForResultCalculation Pnti;
				Pnti.p=p; 
					Pnti.p.x-=pEL.x;
					Pnti.p.y-=pEL.y;

				Pnti.r=sqrt(Pnti.p.x*Pnti.p.x+Pnti.p.y*Pnti.p.y);
				if (Pnti.r<1e-3) Pnti.r=1e-3;
				Pnti.sinfi=Pnti.p.y/Pnti.r;
				Pnti.cosfi=Pnti.p.x/Pnti.r;
				if (!PointList.AddToList(Pnti)) return RETCODE_NOMEM;
			}
		}
		boundsVB.endP=PointList.GetListLength();
	}
	else
	{
		boundsVB.beginP=PointList.GetListLength();
		for (i=0; i<pntVsize; i++)
		{
			const PointXYZ& p=pntV[i];
			PointForResultCalculation Pnti;
			Pnti.p=p; 
				Pnti.p.x-=B.x;
				Pnti.p.y-=B.y;

			Pnti.r=sqrt(Pnti.p.x*Pnti.p.x+Pnti.p.y*Pnti.p.y);
			if (Pnti.r<1e-3) Pnti.r=1e-3;
			Pnti.sinfi=Pnti.p.y/Pnti.r;
			Pnti.cosfi=Pnti.p.x/Pnti.r;
			if (!PointList.AddToList(Pnti)) return RETCODE_NOMEM;
		}
		boundsVB.endP=PointList.GetListLength();
	}
	
	qMN=0;
		
	if (ResCalcPointNum=PointList.GetListLength())
	{
		ResCalcPointArray=new PointForResultCalculation[ResCalcPointNum];
		if (!ResCalcPointArray) return RETCODE_NOMEM;
		PointList.LoadList(ResCalcPointArray);
	}
	else return RETCODE_ERROR;

	for (i=0; i<ResCalcPointNum; i++)
		ResCalcPointArray[i].res=new PointXYZ[1];
	
	ListOfValues<RectForResultCalculation> RectList;

	if ((ResCalcElArray=new ResCalcElInfo[krect])==NULL) return RETCODE_NOMEM;

	logfile<<"vp: searching for 3d points in 2d mesh"<<endl;

	for (j=0; j<ResCalcPointNum; j++)
	{
		PointForResultCalculation& pj=ResCalcPointArray[j];
		PointRZ pjrz=PointRZ(pj.r, pj.p.z);


		i=GetRegElem2DNum(qr, rm, qz, zm, pjrz, 1e-3); //   

		if (i==-1)
		{
			pj.rnum=-1;
			pj.isfound=true;
			continue;
		}

		if (i<1||i>nreg)
			return RETCODE_ERROR;

		const int& elnum=reg[i-1];						//   

		if (elnum<1||elnum>krect)
			return RETCODE_ERROR;

		const Rect& r=rect[elnum-1];
		ResCalcElInfo& ElInfo=ResCalcElArray[elnum-1];

		pj.isfound=true;
		r.calcH(pjrz, pnt, pj.ra, pj.rb, pj.za, pj.zb);

		ElInfo.isResCalc=true;
		if (!ElInfo.ResCalcPoints.AddToList(j))
			return RETCODE_NOMEM;
	}

	logfile<<"vp: searching finished"<<endl;

	for (i=0; i<krect; i++)
	{
		const ResCalcElInfo& ElInfo=ResCalcElArray[i];
		if (ElInfo.isResCalc==false) continue;
		const Rect &r=rect[i];
		const PointRZ& p0=pnt[r.nodes[0]-1];
		const PointRZ& p3=pnt[r.nodes[3]-1];
		RectForResultCalculation RectForCalc;
		RectForCalc.rnum=i; //    
		RectForCalc.hr=p3.r-p0.r;
		RectForCalc.hz=p3.z-p0.z;
		j=RectList.GetListLength();
		if (!RectList.AddToList(RectForCalc))
			return RETCODE_NOMEM;
		ElInfo.ResCalcPoints.ResetScan();
		while (ElInfo.ResCalcPoints.GetCurrentElPtr())
		{
			ResCalcPointArray[ElInfo.ResCalcPoints.GetCurrentElPtr()->v].rnum=j;
			ElInfo.ResCalcPoints.MoveScan();
		}
	}

	delete [] ResCalcElArray;

	for (j=0; j<ResCalcPointNum; j++)
		if (!ResCalcPointArray[j].isfound) return RETCODE_ERROR;
	
	if (ResCalcRectNum=RectList.GetListLength())
	{
		ResCalcRectArray=new RectForResultCalculation[ResCalcRectNum];
		if (ResCalcRectArray==NULL) return RETCODE_NOMEM;
		RectList.LoadList(ResCalcRectArray);
	}
	else return RETCODE_ERROR;

	return 0;
}

/*!      */
int MeshForVP::GetWeights(const double* q, bool storeq)
{
	int i, j;
	double sum;
	if (storeq)
	{
		for (i=0; i<nc; i++)
			qv[i]=q[i];
	}
	for (i=0; i<nc; i++)
		q_all[i]=q[i];
	return RETCODE_OK;
}

/*!     */
int MeshForVP::ReadMtrCatalog(const char* fname, double* cat)
{
	ifstream inf(fname);
	if (!inf) return RETCODE_NOFILE;
	int i;
	double v;
	i=0;
	while (inf&&!inf.eof())
	{
		if (!inf.good()) return RETCODE_BADFILE;
		inf>>i;
		if (!inf.eof()) {
			inf>>v;
			cat[i]=v;
		}
	}
	inf.close();
	return RETCODE_OK;
}

/*!   xyzMN */
int MeshForVP::ReadResMNCalcFile(const char* fname, ListOfValues<PointForResultCalculation>& _listP, ResCalcPointBounds& _boundsP,
#ifndef USEREGINFO
	PointXYZ& gp1, PointXYZ& gp2, // gabarit points
	const bool& first, 
#endif
	const PointXYZ& pEL)
{
	int i, sz;
	ifstream xyzf(fname);
	_boundsP.beginP=_listP.GetListLength();
	if (xyzf)
	{
		xyzf>>sz;
		qMN=sz;
		for (i=0; i<sz; i++)
		{
			PointForResultCalculation PntM, PntN, PntO;
			if (!xyzf.good()) 
				return RETCODE_BADFILE;
			xyzf>>PntM.p>>PntN.p;
			PntO.p=(PntM.p+PntN.p)*0.5;
				PntM.p.x-=pEL.x;
				PntM.p.y-=pEL.y;

				PntN.p.x-=pEL.x;
				PntN.p.y-=pEL.y;

				PntO.p.x-=pEL.x;
				PntO.p.y-=pEL.y;


			PntM.r=sqrt(PntM.p.x*PntM.p.x+PntM.p.y*PntM.p.y);
			if (PntM.r<1e-3) PntM.r=1e-3;
			PntM.sinfi=PntM.p.y/PntM.r;
			PntM.cosfi=PntM.p.x/PntM.r;

			PntN.r=sqrt(PntN.p.x*PntN.p.x+PntN.p.y*PntN.p.y);
			if (PntN.r<1e-3) PntN.r=1e-3;
			PntN.sinfi=PntN.p.y/PntN.r;
			PntN.cosfi=PntN.p.x/PntN.r;

			PntO.r=sqrt(PntO.p.x*PntO.p.x+PntO.p.y*PntO.p.y);
			if (PntO.r<1e-3) PntO.r=1e-3;
			PntO.sinfi=PntO.p.y/PntO.r;
			PntO.cosfi=PntO.p.x/PntO.r;

			if (!_listP.AddToList(PntM)) 
				return RETCODE_NOMEM;
			if (!_listP.AddToList(PntN)) 
				return RETCODE_NOMEM;
			if (!_listP.AddToList(PntO)) 
				return RETCODE_NOMEM;
		}
		xyzf.close();
	}
	_boundsP.endP=_listP.GetListLength();
	return RETCODE_OK;
}

/*!   xyzVector* */
int MeshForVP::ReadResCalcFile(const char* fname, ListOfValues<PointForResultCalculation>& _listP, ResCalcPointBounds& _boundsP,
#ifndef USEREGINFO
	PointXYZ& gp1, PointXYZ& gp2, // gabarit points
	const bool& first, 
#endif
	const PointXYZ& pEL)
{
	int i, sz;
	ifstream xyzf(fname);
	_boundsP.beginP=_listP.GetListLength();
	if (xyzf)
	{
		xyzf>>sz;
		for (i=0; i<sz; i++)
		{
			PointForResultCalculation Pnt;
			if (!xyzf.good()) 
				return RETCODE_BADFILE;
			xyzf>>Pnt.p;
				Pnt.p.x-=pEL.x;
				Pnt.p.y-=pEL.y;

			Pnt.r=sqrt(Pnt.p.x*Pnt.p.x+Pnt.p.y*Pnt.p.y);
			if (Pnt.r<1e-3) Pnt.r=1e-3;
			Pnt.sinfi=Pnt.p.y/Pnt.r;
			Pnt.cosfi=Pnt.p.x/Pnt.r;

#ifndef USEREGINFO
			if (first&&i==0)
			{
				gp1.x=gp2.x=Pnt.p.x;
				gp1.y=gp2.y=Pnt.p.y;
				gp1.z=gp2.z=Pnt.p.z;
			}
			else
			{
				gp1.x=min(gp1.x, Pnt.p.x);
				gp2.x=max(gp2.x, Pnt.p.x);
				gp1.y=min(gp1.y, Pnt.p.y);
				gp2.y=max(gp2.y, Pnt.p.y);
				gp1.z=min(gp1.z, Pnt.p.z);
				gp2.z=max(gp2.z, Pnt.p.z);
			}
#endif
			if (!_listP.AddToList(Pnt)) 
				return RETCODE_NOMEM;
		}
		xyzf.close();
	}
	_boundsP.endP=_listP.GetListLength();
	return RETCODE_OK;
}

/*!    */
int MeshForVP::ReadTimesFile(const char* fname)
{
	ifstream inf(fname);
	int i, sz;
	if (!inf) return RETCODE_NOFILE;
	inf >> sz; inf.ignore(1000, '\n');
	ntimes=sz;
	for (i=0; i<sz; i++)
	{
		if (!inf.good()) return RETCODE_BADFILE;
		inf >> times[i];
	}
	inf.close();
	return RETCODE_OK;
}

#ifdef USEREGINFO
int GetRegElem1DNum(const int& _qc0, const int& _qc, const double* _cm, const double& c, const double& eps)
{
	int _qcc=(_qc0+_qc)>>1;
	if (c<(_cm[_qcc]+eps))
	{
		if ((_qc0+1)==_qcc)
		{
			if ((_cm[_qc0]-eps)<c)
				return _qc0;
			return -1;
		}
		return GetRegElem1DNum(_qc0, _qcc, _cm, c, eps);
	}
	if (_qcc==(_qc-1))
	{
		if (c<(_cm[_qc]+eps))
			return _qcc;
		return -1;
	}
	return GetRegElem1DNum(_qcc, _qc, _cm, c, eps);
}
int GetRegElem2DNum(const int& _qr, const double* _rm,
					const int& _qz, const double* _zm,
					const PointRZ& p, const double& eps)
{
	int i,j;
	for(i=0;i<_qr;i++)
	{
		if(p.r<_rm[i])
		{
			i-=(i>0);
			break;
		}
	}
	for(j=0;j<_qz;j++)
	{
		if(p.z<_zm[j]+eps)
		{
			j-=1;
			break;
		}
	}
	if (i==_qr||j==_qz||j<0)
	{
		logfile<<"Error: point is not in mesh."<<endl;
		return -1;
	}
	return ((_qr-1)*j+i+1);
}
#endif
