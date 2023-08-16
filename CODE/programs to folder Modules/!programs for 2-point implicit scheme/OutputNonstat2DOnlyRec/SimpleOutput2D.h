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
 *  This file contains headers for direct outputing desired primary field components functions
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 January 16, 2023
*/

#pragma once

// The struct is used for direct output of the primary field
struct RectOutData
{
	float hr, hz, psi, eta, psi2, eta2, sigma, sigmaZ;
	int node0, node1, node2, node3, elem;
};

void Output_Field(int nrec, double *v2, double *result, RectOutData *vRecOutData,vector<int> &RecToSourceE,int kpnt);
void Output_Field_Pr(int nrec, double *v2, double *result, RectOutData *vRecOutData,vector<int> &RecToSourceE,int kpnt);
int FindElementsForReceivers(int nrec, int qr, int qz, double *rm, double *zm, int *reg, PointRZ *rec, RectOutData *vRecOutData,
	PointRZ *pnt, Rect *rect,double *sigma,double *sigmaZ,double rmin,double zmin,double rmax,double zmax);
void Output_FieldAr_Pr(int nrec, double *v2, double *result, RectOutData *vRecOutData,vector<int> &RecToSourceE,int kpnt);
void Output_ErForHphi_Pr(int nrec, double *v2, double *result, RectOutData *vRecOutData,vector<int> &RecToSourceE,int krect);
void Output_EzForHphi_Pr(int nrec, double *v2, double *result, RectOutData *vRecOutData,vector<int> &RecToSourceE,int krect,PointRZ *rec);
