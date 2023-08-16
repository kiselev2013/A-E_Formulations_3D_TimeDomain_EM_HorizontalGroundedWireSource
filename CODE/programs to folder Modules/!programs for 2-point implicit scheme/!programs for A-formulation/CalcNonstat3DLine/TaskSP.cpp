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
 *  This file contains helper code for calculating nonstationary 3D VFEM tasks
 *  
 *  
 *  Written by Prof. Yuri G. Soloveichik, Prof. Marina G. Persova, Ph.D. Denis V. Vagin
 *  Novosibirsk State Technical University,                                            
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                  
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                   
 *  Version 2           
*/

#include "stdafx.h"

#include "TaskSP.h"

#include "time_approx.h"

#include"TimeWatcher.h"

using namespace std;

ofstream logfile;

extern double Scal(double *a, double *b, int n);

extern int GetStiffnessMatrixForRect3D(
	const double &x1, const double &y1, const double &z1, // coordinates of local node 0
	const double &x2, const double &y2, const double &z2, // coordinates of local node 7
	double m[8][8]);


extern bool CheckStop(void);

void CopyV(const int &n, const double *from, double *to)
{
	int i;
	for (i=0; i<n; i++)
		to[i]=from[i];
}

double Scal(const int &n, const double *v1, const double *v2)
{
	double s=0;
	int i;
	for (i=0; i<n; i++)
		s+=v1[i]*v2[i];
	return s;
}

void  mult_symmetr(int *ig, int *jg, double *gg, double *di, double* x, double *y, int n)
{
	int i, j, k;

		for(i=0; i<n; i++)
		{
			y[i] = di[i]*x[i];
			for(j=ig[i]; j<=ig[i+1]-1; j++)
			{
				k = jg[j];
				y[i] += gg[j]*x[k];
				y[k] += gg[j]*x[i];
			}
		}
}

struct ListEl 
{
	int v;
	ListEl *next;
	ListEl() { next=NULL; };
};

inline int indX(const int& i)
{
	return div(i, 2).rem;
}

inline int indY(const int& i)
{
	return div(div(i, 2).quot, 2).rem;
}

inline int indZ(const int& i)
{
	return div(i, 4).quot;
}

extern double G1d[2][2];
extern double M1d[2][2];


void Calc_J_optimized(double *x, double *y, double *z, 
					  double (*J_1_T)[3], double *det_J_abs,  int n_of_point)
{
	int i, j;
	double J[3][3];
	double det_J;

	for(i=0; i<3; i++)
		for(j=0; j<3; j++)
			J[i][j] = 0.0;

	for(i=0; i<8; i++)
	{
		J[0][0] += x[i]*gauss_3_d_phi[n_of_point][i][0]; //  d_xi[i];
		J[0][1] += x[i]*gauss_3_d_phi[n_of_point][i][1]; //  d_eta[i];
		J[0][2] += x[i]*gauss_3_d_phi[n_of_point][i][2]; //  d_zeta[i];

		J[1][0] += y[i]*gauss_3_d_phi[n_of_point][i][0]; //  d_xi[i];
		J[1][1] += y[i]*gauss_3_d_phi[n_of_point][i][1]; //  d_eta[i];
		J[1][2] += y[i]*gauss_3_d_phi[n_of_point][i][2]; //  d_zeta[i];

		J[2][0] += z[i]*gauss_3_d_phi[n_of_point][i][0]; //  d_xi[i];
		J[2][1] += z[i]*gauss_3_d_phi[n_of_point][i][1]; //  d_eta[i];
		J[2][2] += z[i]*gauss_3_d_phi[n_of_point][i][2]; //  d_zeta[i];
	}
	det_J = J[0][0]*J[1][1]*J[2][2] - J[0][0]*J[1][2]*J[2][1] + J[1][0]*J[2][1]*J[0][2]
	- J[1][0]*J[0][1]*J[2][2] + J[2][0]*J[0][1]*J[1][2] - J[2][0]*J[1][1]*J[0][2];

	*det_J_abs = fabs(det_J);

	J_1_T[0][0] = (J[1][1]*J[2][2]-J[2][1]*J[1][2])/det_J;
	J_1_T[1][0] = (J[2][1]*J[0][2]-J[0][1]*J[2][2])/det_J;
	J_1_T[2][0] = (J[0][1]*J[1][2]-J[1][1]*J[0][2])/det_J;
	J_1_T[0][1] = (-J[1][0]*J[2][2]+J[2][0]*J[1][2])/det_J;
	J_1_T[1][1] = (J[0][0]*J[2][2]-J[2][0]*J[0][2])/det_J;
	J_1_T[2][1] = (-J[0][0]*J[1][2]+J[1][0]*J[0][2])/det_J;
	J_1_T[0][2] = (J[1][0]*J[2][1]-J[2][0]*J[1][1])/det_J;
	J_1_T[1][2] = (-J[0][0]*J[2][1]+J[2][0]*J[0][1])/det_J;
	J_1_T[2][2] = (J[0][0]*J[1][1]-J[1][0]*J[0][1])/det_J;
}

void Calc_Hex_Local_Matrix_Dx(double *x, double *y, double *z, double (*dx)[8])
{
	int i, j, i1, j1;
	double phi1;
	double gauss_3_mult;
	double J_1_T[3][3];
	double det_J_abs;
	double grad_all[8]; //   .    

	for(i=0; i<8; i++) //  
		for(j=0; j<8; j++)
		{
			dx[i][j] = 0.0;
		}

		for(i=0; i<27; i++) //    
		{
			Calc_J_optimized(x, y, z, J_1_T, &det_J_abs, i); //  
			gauss_3_mult = gauss_3_A_all[i]*det_J_abs; // A_i*A_j*A_k*|J|

			for(j=0; j<8; j++) //    -   ...
				grad_all[j]    = J_1_T[0][0]*gauss_3_d_phi[i][j][0]+  
								 J_1_T[0][1]*gauss_3_d_phi[i][j][1]+  
								 J_1_T[0][2]*gauss_3_d_phi[i][j][2];

			for(i1=0; i1<8; i1++)
			{
				phi1 = gauss_3_phi[i][i1]*gauss_3_mult;	
				for(j1=0; j1<8; j1++)
				{
					dx[i1][j1] += grad_all[j1]*phi1;
				}//j1
			}//i1
		}// i

	for(i=0; i<8; i++) // transpose
		for(j=0; j<i; j++)
			swap(dx[i][j], dx[j][i]);			
}
void Calc_Hex_Local_Matrix_Dy(double *x, double *y, double *z, double (*dy)[8])
{
	int i, j, i1, j1;
	double phi1;
	double gauss_3_mult;
	double J_1_T[3][3];
	double det_J_abs;
	double grad_all[8]; //   .    

	for(i=0; i<8; i++) //  
		for(j=0; j<8; j++)
		{
			dy[i][j] = 0.0;
		}

		for(i=0; i<27; i++) //    
		{
			Calc_J_optimized(x, y, z, J_1_T, &det_J_abs, i); //  
			gauss_3_mult = gauss_3_A_all[i]*det_J_abs; // A_i*A_j*A_k*|J|

			for(j=0; j<8; j++) //    -   ...
				grad_all[j]   =		J_1_T[1][0]*gauss_3_d_phi[i][j][0]+  
									J_1_T[1][1]*gauss_3_d_phi[i][j][1]+  
									J_1_T[1][2]*gauss_3_d_phi[i][j][2];

			for(i1=0; i1<8; i1++)
			{
				phi1 = gauss_3_phi[i][i1]*gauss_3_mult;	
				for(j1=0; j1<8; j1++)
				{
					dy[i1][j1] += grad_all[j1]*phi1;
				}//j1
			}//i1
		}// i

	for(i=0; i<8; i++) // transpose
		for(j=0; j<i; j++)
			swap(dy[i][j], dy[j][i]);
}
void Calc_Hex_Local_Matrix_Dz(double *x, double *y, double *z, double (*dz)[8])
{
	int i, j, i1, j1;
	double phi1;
	double gauss_3_mult;
	double J_1_T[3][3];
	double det_J_abs;
	double grad_all[8]; //   .    

	for(i=0; i<8; i++) //  
		for(j=0; j<8; j++)
		{
			dz[i][j] = 0.0;
		}

		for(i=0; i<27; i++) //    
		{
			Calc_J_optimized(x, y, z, J_1_T, &det_J_abs, i); //  
			gauss_3_mult = gauss_3_A_all[i]*det_J_abs; // A_i*A_j*A_k*|J|

			for(j=0; j<8; j++) //    -   ...
				grad_all[j] =		J_1_T[2][0]*gauss_3_d_phi[i][j][0]+  
									J_1_T[2][1]*gauss_3_d_phi[i][j][1]+  
									J_1_T[2][2]*gauss_3_d_phi[i][j][2];

			for(i1=0; i1<8; i1++)
			{
				phi1 = gauss_3_phi[i][i1]*gauss_3_mult;	
				for(j1=0; j1<8; j1++)
				{
					dz[i1][j1] += grad_all[j1]*phi1;
				}//j1
			}//i1
		}// i

	for(i=0; i<8; i++) // transpose
		for(j=0; j<i; j++)
			swap(dz[i][j], dz[j][i]);
}

void Mult_Plot(double *a, double *x, double *y, int n)
{
	int i, j, temp;
	double sum;

	for(i=0;i<n;i++)
	{
		sum = 0.0;
		temp = i*n;
		for(j=0;j<n;j++)
			sum += a[temp+j]*x[j];

		y[i] = sum;
	}
}

void Calc_Hex_Local_Matrix_C(double *x, double *y, double *z, double (*c)[8]) //  
							 
{
	int i, j, i1, j1;
	double gauss_3_mult;
	double J_1_T[3][3];
	double det_J_abs;
	double grad_all[8][3]; //   .    

	for(i=0; i<8; i++) //  
		for(j=0; j<8; j++)
		{
			c[i][j] = 0.0;
		}

	for(i=0; i<27; i++) //    
	{
		Calc_J_optimized(x, y, z, J_1_T, &det_J_abs, i); //  
		gauss_3_mult = gauss_3_A_all[i]*det_J_abs; // A_i*A_j*A_k*|J|

		for(j=0; j<8; j++) //    -   ...
			Mult_Plot((double*)J_1_T, (double*)gauss_3_d_phi[i][j],	(double*)grad_all[j], 3);

		for(i1=0; i1<8; i1++)
		{
			for(j1=0; j1<8; j1++)
			{
				c[i1][j1] += Scal((double*)grad_all[i1],(double*)grad_all[j1], 3)*gauss_3_mult;
			}//j1
		}//i1
	}// i
}

void AddV(const int &n, const double *from, double *to, const double& mlt=1)
{
	for (int i=0; i<n; i++)
		to[i]+=mlt*from[i];
}

void LocalMatrixOnVector(const double m[4][4], const double v[4], double res[4])
{
	int i, j;
	for (i=0; i<4; i++)
	{
		res[i]=0;
		for (j=0; j<4; j++)
			res[i]+=m[i][j]*v[j];
	}
}

int FindTimeLayer(int nt,double *times,double val)
{
	int it;

	if(nt<1)
	{
		return -1;
	}

	if(fabs(val)<1e-14 && fabs(times[0])<1e-14)
	{
		return 0;
	}

	for(it=1;it<nt;it++)
	{
		if(fabs(1.0-times[it]/val)<1e-6)
		{
			return it;
		}
	}

	return -1;
}

void GenerateDecInfo(int ntime,double *time)
{
	int iDec,nDec,i,j;
	double htpre,htcur;
	int max_dec_size;
	ofstream ofp;
	vector<int> DecIg;

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

	ofp.open("nDec");
	ofp<<nDec<<'\n';
	ofp.close();
	ofp.clear();

	ofp.open("DecIg");
	ofp<<0<<'\n';
	for(iDec=0;iDec<nDec;iDec++){ofp<<DecIg[iDec+1]<<'\n';}
	ofp.close();
	ofp.clear();

	max_dec_size=0;
	for(iDec=0;iDec<nDec;iDec++)
	{
		j=DecIg[iDec+1]-DecIg[iDec];
		if(j>max_dec_size)
		{
			max_dec_size=j;
		}
	}
	ofp.open("max_dec_size");
	ofp<<max_dec_size<<'\n';
	ofp.close();
	ofp.clear();
}

int ReadDecInfo(int &nDec,int &max_dec_size,vector<int> &DecIg)
{
	ifstream inf;
	int i;

	inf.open("max_dec_size");
	if(!inf)
	{
		logfile<<"Error in open file "<<"max_dec_size"<<endl;
		cout<<"Error in open file "<<"max_dec_size"<<endl;
		return 1;
	}
	inf>>max_dec_size;
	inf.close();
	inf.clear();

	inf.open("nDec");
	if(!inf)
	{
		logfile<<"Error in open file "<<"nDec"<<endl;
		cout<<"Error in open file "<<"nDec"<<endl;
		return 1;
	}
	inf>>nDec;
	inf.close();
	inf.clear();

	cout<<"nDec= "<<nDec<<endl;

	DecIg.resize(nDec+1);

	inf.open("DecIg");
	if(!inf)
	{
		logfile<<"Error in open file "<<"DecIg"<<endl;
		cout<<"Error in open file "<<"DecIg"<<endl;
		return 1;
	}
	for(i=0;i<(nDec+1);i++){inf>>DecIg[i];}
	inf.close();
	inf.clear();
	
	return 0;
}

int CalcSP()
{
	char buf[256];

	int i,j;
	bool iVP;
	ifstream inf;
	Time_approx_for_vfem *ta=NULL;
	int retv;
	vector<int> RecvPlsIgB,RecvPlsIgE,RecToSourceB,RecToSourceE;
	bool fstop;

	int tbeg_pre;

	sprintf(buf,"Start CalcSP");
	TimeWatcher::getInstance().AddTimeSpot(clock(),buf);

	if(CheckStop()){return 1;}

	iVP=false;
	fstop=false;

	int npls,ngrp,p1,p2;

	npls=0;
	inf.open("group");
	if(!inf)
	{
		cout<<"Error in open file "<<"group"<<endl;
		return 1;
	}
	inf>>ngrp;
	for(i=0;i<ngrp;i++)
	{
		inf>>j>>p1>>p2;
		npls+=p2-p1+1;
	}
	inf.close();
	inf.clear();

	RecvPlsIgB.resize(npls+1);
	RecvPlsIgE.resize(npls+1);

	inf.open("recvsb");
	if(!inf)
	{
		logfile<<"Error in open file "<<"recvsb"<<endl;
		cout<<"Error in open file "<<"recvsb"<<endl;
		return 1;
	}
	RecvPlsIgB[0]=0;
	for(i=0;i<npls;i++){inf>>RecvPlsIgB[i+1];}
	inf.close();
	inf.clear();

	inf.open("recvse");
	if(!inf)
	{
		logfile<<"Error in open file "<<"recvse"<<endl;
		cout<<"Error in open file "<<"recvse"<<endl;
		return 1;
	}
	RecvPlsIgE[0]=0;
	for(i=0;i<npls;i++){inf>>RecvPlsIgE[i+1];}
	inf.close();
	inf.clear();

	for(i=0;i<npls;i++)
	{
		RecvPlsIgB[i+1]+=RecvPlsIgB[i];
		RecvPlsIgE[i+1]+=RecvPlsIgE[i];
	}

	RecToSourceB.resize(RecvPlsIgB[npls]);
	RecToSourceE.resize(RecvPlsIgE[npls]);

	for(i=0;i<npls;i++)
	{
		for(j=RecvPlsIgB[i];j<RecvPlsIgB[i+1];j++)
		{
			RecToSourceB[j]=i;
		}
	}

	for(i=0;i<npls;i++)
	{
		for(j=RecvPlsIgE[i];j<RecvPlsIgE[i+1];j++)
		{
			RecToSourceE[j]=i;
		}
	}

	ta=new Time_approx_for_vfem(npls);
	ta->iVP=iVP;

	ta->AnomalType=0;
	inf.open("AnomalType");
	if(inf)
	{
		inf>>ta->AnomalType;
		inf.close();
	}
	inf.clear();

	ta->ForInv=0;
	inf.open("sigmafit");
	if(inf)
	{
		inf>>ta->ForInv;
		inf.close();
	}
	inf.clear();
	
	ta->nthreads=1;
	if(!ta->AnomalType)
	{
		inf.open("nthreads.txt");
		if(inf)
		{
			inf>>ta->nthreads;
			inf.close();
		}
		inf.clear();

		inf.open("nthreads");
		if(inf)
		{
			inf>>ta->nthreads;
			inf.close();
		}
		inf.clear();
	}
	else
	{
		inf.open("nProcAnom.txt");
		if(inf)
		{
			inf>>ta->nthreads;
			inf.close();
		}
		inf.clear();

		inf.open("ia");
		if(!inf)
		{
			cout<<"Error in open file "<<"ia"<<endl;
			logfile<<"Error in open file "<<"ia"<<endl;
			exit(1);
		}
		inf>>ta->ia;
		inf.close();
		inf.clear();
	}

	ta->Read_data(RecvPlsIgB,RecvPlsIgE);

	int iDec,nDec,npre1,npre2,nlst;
	vector<int> DecIg;


	int max_dec_size;


	GenerateDecInfo(ta->ntime,ta->time);

	retv=ReadDecInfo(nDec,max_dec_size,DecIg);
	if(retv)
	{
		logfile<<"Error in function "<<"ReadDecInfo"<<endl;
		cout<<"Error in function "<<"ReadDecInfo"<<endl;
		return retv;
	}

	ta->max_dec_size=max_dec_size;

	sprintf(buf,"Start Prepare_To_Schema");
	TimeWatcher::getInstance().AddTimeSpot(clock(),buf);

	ta->Prepare_To_Schema();

	sprintf(buf,"Finish Prepare_To_Schema");
	TimeWatcher::getInstance().AddTimeSpot(clock(),buf);

	for(iDec=0;iDec<nDec;iDec++)
	{
		sprintf(buf,"Start decade iDec= %d",iDec);
		TimeWatcher::getInstance().AddTimeSpot(clock(),buf);

		Subdomain::fdfct=!iDec;

		tbeg_pre=(iDec>0)? DecIg[iDec-1] : 0;
		
		npre1=-1;
		npre2=-1;
		if(iDec)
		{
			npre1=FindTimeLayer(ta->ntime,ta->time,ta->time[DecIg[iDec]]-1*(ta->time[DecIg[iDec]+1]-ta->time[DecIg[iDec]]));
			npre2=FindTimeLayer(ta->ntime,ta->time,ta->time[DecIg[iDec]]-2*(ta->time[DecIg[iDec]+1]-ta->time[DecIg[iDec]]));
			if(npre1<0 || npre2<0)
			{
				cout<<"Error in npre1 or npre2"<<endl;
				logfile<<"Error in npre1 or npre2"<<endl;
				return 1;
			}
		}

		nlst=DecIg[iDec+1];
		if(iDec<nDec-1)
		{
			nlst=FindTimeLayer(ta->ntime,ta->time,ta->time[DecIg[iDec+1]]-2*(ta->time[DecIg[iDec+1]+1]-ta->time[DecIg[iDec+1]]));
		}

		sprintf(buf,"Start schema iDec= %d",iDec);
		TimeWatcher::getInstance().AddTimeSpot(clock(),buf);

		ta->Three_Layers_Schema(DecIg[iDec],DecIg[iDec+1],npre1,npre2,nlst,tbeg_pre,(iDec==nDec-1));

		sprintf(buf,"Processing decade results iDec= %d",iDec);
		TimeWatcher::getInstance().AddTimeSpot(clock(),buf);

		if(iDec==nDec-1)
		{
			for(i=0;i<nDec-1;i++)
			{
				npre1=-1;
				npre2=-1;
				if(i)
				{
					npre1=FindTimeLayer(ta->ntime,ta->time,ta->time[DecIg[i]]-1*(ta->time[DecIg[i]+1]-ta->time[DecIg[i]]));
					npre2=FindTimeLayer(ta->ntime,ta->time,ta->time[DecIg[i]]-2*(ta->time[DecIg[i]+1]-ta->time[DecIg[i]]));
					if(npre1<0 || npre2<0)
					{
						cout<<"Error in npre1 or npre2 in pre out"<<endl;
						logfile<<"Error in npre1 or npre2 in pre out"<<endl;
						return 1;
					}
				}

				ta->InitvPre(DecIg[i],DecIg[i+1],npre1,npre2);

				retv=ta->LoadPreviousResults(i,DecIg[i]-1,DecIg[i+1]-1,RecToSourceB,RecToSourceE);
				if(retv)
				{
					cout<<"Function LoadPreviousResults returned "<<retv<<endl;
					logfile<<"Function LoadPreviousResults returned "<<retv<<endl;
					return 1;
				}
			}
			ta->FinishCalculation();
		}
		else
		{
			ta->SaveCurrentResults(iDec,DecIg[iDec]-1,DecIg[iDec+1]-1,RecToSourceB,RecToSourceE);
		}

		sprintf(buf,"Finish decade iDec= %d",iDec);
		TimeWatcher::getInstance().AddTimeSpot(clock(),buf);
	}

	if(ta)
	{
		delete ta;
		ta=NULL;
	}

	return 0;
}
