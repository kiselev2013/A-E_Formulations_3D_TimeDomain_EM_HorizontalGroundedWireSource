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
 *  This file contains headers for smooth outputing desired primary field components functions
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 January 16, 2023
*/

#pragma once


struct Output2D_PAR
{
	int fl;
	int lem[4];
	float coef[4][4];
	float ves;
	float v[2][3];
};

struct Output2D
{
	Output2D_PAR or;
	Output2D_PAR oz;
	Output2D_PAR orfz[4],ozfr[4];
};

// The struct is used to smooth the output of the primary field
struct SmoothOutput2D
{
	int nrec,kpnt,krect;
	PointRZ *rec;
	PointRZ *pnt;
	Rect *rect;
	vector<int> RecToElem;
	double *sigma,*sigmaZ;
	double *mu;
	int nreg;
	int qr,qz;
	int *reg;
	double *rm,*zm;

	Output2D *OutRecList;

	SmoothOutput2D();
	~SmoothOutput2D();
	int Init(int p_nrec,int p_kpnt,int p_krect,PointRZ *p_rec,PointRZ *p_pnt,Rect *p_rect,int MatFlag,int *matrec,double *p_sigma,double *p_sigmaZ,
		double *p_mu,int p_nreg,int p_qr,int p_qz,int *p_reg,double *p_rm,double *p_zm, int stype);
	int Init(int p_nrec,int p_kpnt,int p_krect,PointRZ *p_rec,PointRZ *p_pnt,Rect *p_rect,int MatFlag,int *matrec,double *p_sigma,double *p_sigmaZ,
		double *p_mu,int p_nreg,int p_qr,int p_qz,int *p_reg,double *p_rm,double *p_zm, int stype,vector<int> &OutCompFlag);
	void GetSmoothForR(int i,int j,int k,int ef,int l,int MatFlag,PointRZ &reck,int rlem[4],float rcoef[4][4],float rv[2][3],float &ves,int &fl,int stype);
	void GetSmoothForZ(int i,int j,int k,int ef,int l,int MatFlag,PointRZ &reck,int zlem[4],float zcoef[4][4],float zv[2][3],float &ves,int &fl,int stype);
};

struct SmoothOutput2D_Ax : public SmoothOutput2D
{
	void OutputE(double *v2,double *result,int OutType,vector<int> &RecToSourceE);
	void OutputB(double *v2,double *resultx,double *resulty,double *resultz,vector<PointXYZ> &Ic,double *sinfi,double *cosfi,int OutType,vector<int> &RecToSourceB);
};

struct SmoothOutput2D_Ar : public SmoothOutput2D
{
	void Output_dArdz(double *v2,double *result,int OutType,vector<int> &RecToSourceB);
	void Output_Er(double *v2a,double *result,int OutType,vector<int> &RecToSourceE);
	void Output_Er_Pr(double *v2a,double *result,int OutType,vector<int> &RecToSourceE,vector<int> &OutCompFlag);
};

struct SmoothOutput2D_Hfi : public SmoothOutput2D
{
	void Output_Er(double *v2,double *result,int OutType,vector<int> &RecToSourceE);
	void Output_Ez(double *v2,double *result,int OutType,vector<int> &RecToSourceE);
	void Output_Er_Pr(double *v2,double *result,int OutType,vector<int> &RecToSourceE,vector<int> &OutCompFlag);
	void Output_Ez_Pr(double *v2,double *result,int OutType,vector<int> &RecToSourceE,vector<int> &OutCompFlag);
	void Output_Bphi(double *v2a,double *result,int OutType,vector<int> &RecToSourceB);
	double GetSmoothValForZ(int k,int jj,double *v2,int OutType);
	double GetSmoothValForR(int k,int jj,double *v2,int OutType);
};
