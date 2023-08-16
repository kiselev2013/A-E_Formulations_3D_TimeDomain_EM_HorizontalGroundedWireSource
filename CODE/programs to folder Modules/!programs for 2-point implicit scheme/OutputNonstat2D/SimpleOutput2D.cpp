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
 *  This file contains code for direct outputing desired primary field components functions
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 January 16, 2023
*/

#include "stdafx.h"
#include "SimpleOutput2D.h"
#include "task2d.h"

extern ofstream logfile;

void Output_Field(int nrec, double *v2, double *result, RectOutData *vRecOutData,vector<int> &RecToSourceE,int kpnt)
{
	int k, node0, node1, node2, node3;
	double psi, eta, psi2, eta2;
	int ipls;
	double *v2a;

	for(k=0;k<nrec;k++)
	{
		RectOutData &rod = vRecOutData[k];

		ipls=RecToSourceE[k/NUMBEROFABDIPOLES];
		v2a=v2+kpnt*ipls;

		node0 = rod.node0;
		node1 = rod.node1;
		node2 = rod.node2;
		node3 = rod.node3;

		psi = rod.psi;
		eta = rod.eta;

		psi2 = rod.psi2;
		eta2 = rod.eta2;

		result[k] = psi2*eta2*v2a[node0] + psi*eta2*v2a[node1] + psi2*eta*v2a[node2] + psi*eta*v2a[node3];
	}
}

void Output_Field_Pr(int nrec, double *v2, double *result, RectOutData *vRecOutData,vector<int> &RecToSourceE,int kpnt)
{
	int k, node0, node1, node2, node3;
	double psi, eta, psi2, eta2;
	int ipls;
	double *v2a;

	#pragma omp parallel private(k,node0,node1,node2,node3,psi,eta,psi2,eta2,ipls,v2a)
	{
		#pragma omp for
		for(k=0;k<nrec;k++)
		{
			RectOutData &rod = vRecOutData[k];

			ipls=RecToSourceE[k/NUMBEROFABDIPOLES];
			v2a=v2+kpnt*ipls;

			node0 = rod.node0;
			node1 = rod.node1;
			node2 = rod.node2;
			node3 = rod.node3;

			psi = rod.psi;
			eta = rod.eta;

			psi2 = rod.psi2;
			eta2 = rod.eta2;

			result[k] = psi2*eta2*v2a[node0] + psi*eta2*v2a[node1] + psi2*eta*v2a[node2] + psi*eta*v2a[node3];
		}
	}
}

void Output_FieldAr_Pr(int nrec, double *v2, double *result, RectOutData *vRecOutData,vector<int> &RecToSourceE,int kpnt)
{
	int k, node0, node1, node2, node3;
	double psi, eta, psi2, eta2;
	int ipls;
	double *v2a;

	#pragma omp parallel private(k,node0,node1,node2,node3,psi,eta,psi2,eta2,ipls,v2a)
	{
		#pragma omp for
		for(k=0;k<nrec;k++)
		{
			RectOutData &rod = vRecOutData[k];

			ipls=RecToSourceE[k/2];
			v2a=v2+kpnt*ipls;

			node0 = rod.node0;
			node1 = rod.node1;
			node2 = rod.node2;
			node3 = rod.node3;

			psi = rod.psi;
			eta = rod.eta;

			psi2 = rod.psi2;
			eta2 = rod.eta2;

			result[k] = psi2*eta2*v2a[node0] + psi*eta2*v2a[node1] + psi2*eta*v2a[node2] + psi*eta*v2a[node3];
		}
	}
}

void Output_ErForHphi_Pr(int nrec, double *v2, double *result, RectOutData *vRecOutData,vector<int> &RecToSourceE,int krect)
{
	int k, i;
	double psi;
	int ipls;
	double *v2a;

	#pragma omp parallel private(k,i,psi,ipls,v2a)
	{
		#pragma omp for
		for(k=0;k<nrec;k++)
		{
			RectOutData &rod = vRecOutData[k];

			ipls=RecToSourceE[k/2];
			v2a=v2+4*krect*ipls;

			psi = rod.psi;

			i=rod.elem;

			result[k]=-(-(1.0-psi)*v2a[4*i+0]-psi*v2a[4*i+1]+(1.0-psi)*v2a[4*i+2]+psi*v2a[4*i+3])/rod.hz;

			result[k]/=rod.sigma;
		}
	}
}

void Output_EzForHphi_Pr(int nrec, double *v2, double *result, RectOutData *vRecOutData,vector<int> &RecToSourceE,int krect,PointRZ *rec)
{
	int k, i;
	double psi, eta, psi2, eta2;
	int ipls;
	double *v2a;

	#pragma omp parallel private(k,i,psi,eta,psi2,eta2,ipls,v2a)
	{
		#pragma omp for
		for(k=0;k<nrec;k++)
		{
			RectOutData &rod = vRecOutData[k];

			ipls=RecToSourceE[k/2];
			v2a=v2+4*krect*ipls;

			psi = rod.psi;
			eta = rod.eta;

			psi2 = rod.psi2;
			eta2 = rod.eta2;

			i=rod.elem;

			result[k]=(psi2*eta2*v2a[4*i+0] + psi*eta2*v2a[4*i+1] + psi2*eta*v2a[4*i+2] + psi*eta*v2a[4*i+3])/rec[k].r+
				(-(1.0-eta)*v2a[4*i+0]+(1.0-eta)*v2a[4*i+1]-eta*v2a[4*i+2]+eta*v2a[4*i+3])/rod.hr;
		
			result[k]/=rod.sigmaZ;
		}
	}
}

int FindElementsForReceivers(int nrec, int qr, int qz, double *rm, double *zm, int *reg, PointRZ *rec, RectOutData *vRecOutData, 
	PointRZ *pnt, Rect *rect, double *sigma, double *sigmaZ,double rmin, double  zmin, double  rmax, double  zmax)
{
	bool fstop;
	int i, j, k, l, m;
	int imin, imax, jmin, jmax;

	if(!nrec)
	{
		return 0;
	}

	imin = FindIntervalInDoubleMas(rm, qr, rmin);
	imax = FindIntervalInDoubleMas(rm, qr, rmax);
	jmin = FindIntervalInDoubleMas(zm, qz, zmin);
	jmax = FindIntervalInDoubleMas(zm, qz, zmax);

	if (imin == -1)
	{
		logfile << "No element for rmin= " << rmin << '\n';
		return 1;
	}

	if (imax == -1)
	{
		logfile << "No element for rmax= " << rmax << '\n';
		return 1;
	}

	if (jmin == -1)
	{
		logfile << "No element for zmin= " << zmin << '\n';
		return 1;
	}

	if (jmax == -1)
	{
		logfile << "No element for zmax= " << zmax << '\n';
		return 1;
	}

	fstop = false;
	for (k = 0; (k < nrec && !fstop); k++)
	{
		RectOutData &rod = vRecOutData[k];
		i = FindIntervalInDoubleMas(rm, imin, imax+1, rec[k].r);
		j = FindIntervalInDoubleMas(zm, jmin, jmax+1, rec[k].z);
		if (i != -1 && j != -1)
		{
			m = j*(qr - 1) + i;
			l = reg[m] - 1;

			rod.elem=l;
			rod.sigma=sigma[rect[l].mtr-1];
			rod.sigmaZ=sigmaZ[rect[l].mtr-1];

			rod.node0 = rect[l].nodes[0] - 1;
			rod.node1 = rect[l].nodes[1] - 1;
			rod.node2 = rect[l].nodes[2] - 1;
			rod.node3 = rect[l].nodes[3] - 1;

			rod.hr = pnt[rod.node3].r - pnt[rod.node0].r;
			rod.hz = pnt[rod.node3].z - pnt[rod.node0].z;
			rod.psi = (rec[k].r - pnt[rod.node0].r) / rod.hr;
			rod.eta = (rec[k].z - pnt[rod.node0].z) / rod.hz;
			rod.psi2 = 1.0 - rod.psi;
			rod.eta2 = 1.0 - rod.eta;
		}
		else
		{
			logfile << "No element for reciver " << k + 1 << " r= " << rec[k].r << " z= " << rec[k].z << '\n';
			fstop = true;
		}
	}
	return fstop;
}
