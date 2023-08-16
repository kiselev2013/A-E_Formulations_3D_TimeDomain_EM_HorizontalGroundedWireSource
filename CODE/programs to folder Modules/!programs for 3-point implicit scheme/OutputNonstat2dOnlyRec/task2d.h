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
 *  This file contains the headers for working with a 2D mesh and solution
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 January 16, 2023
*/

#pragma once
#include "utils_direct_solver.h"
#include "SimpleOutput2D.h"
#include "SmoothOutput2D.h"

struct SqLoop
{
	PointXYZ A,B;
};

bool inTriangle(PointRZ *trg,const PointRZ &ptmp);
bool CalcBaseFunctions(PointRZ *trg,double a[3][3]);
double scal(const PointRZ &p1,const PointRZ &p2);
double dist(const PointRZ &p1,const PointRZ &p2);

// Base class for working with the results of calculating two-dimensional problems
struct task2d
{
	char path[256];

	int kpnt, krect, nc;
	PointRZ* pnt;
	Rect* rect;

	int ntimes;
	double *times;

	int nmat3d;
	int *mtr3d2d;

	int nmat2d;
	double *sigma,*sigmaZ;
	double *mu;

	int nreg, qr, qz;
	int *reg;
	double *rm, *zm;

	task2d();
	~task2d();
	void SetPath(char *_path);
	int Read(int fSigmaZ);
	int ReadTimesFile();
	int ReadMtr3d2d(char *path3d);
	int	GetField(int it,double *A,char *str);

	int iDec,nDec,npre1,npre2,tbeg,tend,ntimesdec;
	vector<int> DecIg;

	ifstream inf;
	ofstream ofp,ofpstat;
	char buf[256];
	int ipls,npls;
};

struct task2d_Ax : public task2d
{
	double **A,**Bx,**By,**Bz,**E,*H,**Fx,**Fy,**Fz,*E0;
	double U[3],ux,uy,uz;
	double ftime, ltime;
	int time_gaps_for_loop;
	PointXYZ *xyzVectorE,*xyzVectorB,*xyzVectorE0;
	double x,y,z,r,len;
	double cosphi,sinphi;
	int npntEor,npntBor,npntE0or;
	int npntE,npntB,npntE0;
	int *matrec,*matrec0;
	PointRZ *recE,*recB,*recE0;
	double *sinfiE,*cosfiE;
	double *sinfiB,*cosfiB;
	double *sinfiE0,*cosfiE0;
	double DminSrc;
	PointXYZ Tp;
	int p1,p2,nprof,npls_main,ipm,npc;
	vector<SqLoop> GenSq;
	vector<int> RecvPlsIgB, RecvPlsIgE, RecvPlsIgE0, RecToSourceB, RecToSourceE, RecToSourceE0;
	vector<double> cosphiB,sinphiB,cosphiE,sinphiE,cosphiE0,sinphiE0,lenB,lenE,lenE0;
	vector<int> nGenByPls,GenToPls,NrecB,NrecE;
	vector<PointXYZ> TgCompE0;
	char PathTo2d[1024];
	int nthreads;
	PointXYZ pA,pB,pT,vAB,Ic;
	vector<PointXYZ> IcE,IcE0,IcB;
	SmoothOutput2D_Ax SobjB,SobjE;
	double rc0;
	RectOutData *vRecOutData;
	double rmin,zmin,rmax,zmax;
	int i,j,idp,k,l,m,it,nt,retp,kk;
	
	int init();
	void output(int itime,vector<float> &enor);
	void output_rec();
	void clear();
	void GetSource(PointXYZ &GenMin,PointXYZ &GenMax,PointXYZ &pA,PointXYZ &pB,PointXYZ &vAB,PointXYZ &Ic,double &len,double &cosphi,double &sinphi);
};

struct task2d_Ar : public task2d
{
	double **B,**E,**A,**F,*H,*E0;
	double ftime, ltime;
	int time_gaps_for_loop;
	ifstream inf;
	ofstream ofp;
	PointXYZ *xyzVectorE,*xyzVectorB,*xyzVectorE0;
	double x,y,z;
	char str[256];
	PointXYZ Av,Bv,Tp;			//    
	int npntEor,npntBor,npntE0or;
	int npntE,npntB,npntE0;
	int *matrec,*matrec0;
	PointRZ *recE,*recB,*recE0;
	double *sinfiE,*cosfiE;
	double *sinfiB,*cosfiB;
	double *sinfiE0,*cosfiE0;
	double DminSrc;
	char buf[256];
	int p1,p2,nprof,npls_main,ipm,npc;
	vector<SqLoop> GenSq;
	vector<int> RecvPlsIgB, RecvPlsIgE, RecvPlsIgE0, RecToSourceB, RecToSourceE, RecToSourceE0;
	vector<int> nGenByPls,NrecB,NrecE;
	vector<PointXYZ> TgCompE0;
	vector<int> TgCompFlagX,TgCompFlagY,TgCompFlagZ,TgCompFlagXAll,TgCompFlagYAll,TgCompFlagZAll,TgCompFlagRAll;
	char PathTo2d[1024];
	SmoothOutput2D_Ar SobjB,SobjE/*,SobjE0*/;
	double rc0;
	int i,j,k,l,m,it,nt,retp;
	RectOutData *vRecOutData;

	int init();
	void output(int itime,vector<float> &enor);
	void output_rec();
	void clear();
};

struct task2d_Hfi : public task2d
{
	double **B,**Er,**Ez,**Hphi,*v2,**F,*THphi,*Er0,*Ez0;
	double ftime, ltime;
	int time_gaps_for_loop;
	ifstream inf;
	ofstream ofp;
	PointXYZ *xyzVectorE,*xyzVectorB,*xyzVectorE0;
	double x,y,z;
	char str[256];
	float fv;
	PointXYZ Av,Bv,Tp;
	int npntEor,npntBor,npntE0or;
	int npntE,npntB,npntE0;
	int *matrec,*matrec0;
	PointRZ *recE,*recB,*recE0;
	double *sinfiE,*cosfiE;
	double *sinfiB,*cosfiB;
	double *sinfiE0,*cosfiE0;
	double DminSrc;
	int *DsFlgNode;
	int p1,p2,nprof,npls_main,ipm,npc;
	vector<SqLoop> GenSq;
	vector<int> RecvPlsIgB, RecvPlsIgE, RecvPlsIgE0, RecToSourceB, RecToSourceE, RecToSourceE0;
	vector<int> nGenByPls,nGels,nVels,GenToPls,GenType,GenTypeTrue,NrecB,NrecE;
	int nGelsAll;
	vector<PointXYZ> TgCompE0;
	vector<int> TgCompFlagX,TgCompFlagY,TgCompFlagZ,TgCompFlagXAll,TgCompFlagYAll,TgCompFlagZAll,TgCompFlagRAll;
	char PathTo2d[1024];
	SmoothOutput2D_Hfi SobjB,SobjE/*,SobjE0*/;
	double rc0;
	int i,j,k,l,m,it,nt,retp,kk;
	RectOutData *vRecOutData;

	vector<int> kt1dc;
	vector<vector<int>> l1dc;

	int init();
	void output(int itime,vector<float> &enor);
	void output_rec();
	void clear();

	int ReadDiscontinues(int npls);
	void RepairGeneratorOrder(vector<int> &GenToPls,vector<int> &nGels,vector<int> &nVels,vector<SqLoop> &GenSq,int n,int npls);
	void RepairFieldOrder(vector<int> &GenToPls,vector<int> &nGels,vector<int> &nVels,double *(&Hphi),double *(&THphi),int krect,int n,int npls);
	int GetField(int it,double *v2,double *(&Hphi),char *str,int *DsFlgNode,
			 double *(&THphi),vector<int> &nGels,vector<int> &nVels,vector<int> &GenToPls,int npls_main);
};
