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
 *  This file contains the code for working with edge mesh in 3D VFEM
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin                         
 *  Novosibirsk State Technical University,                 
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia       
 *  vdv_wk@mail.ru                                          
 *  Version 2.0 January 16, 2023                            
*/


#include "stdafx.h"
#include "T_Mapping.h"
double Scal(double *a, double *b, int n)
{
	long i;
	double sum = 0.0;

	for(i=0; i<n; i++)
		sum += a[i]*b[i];

	return sum;
}
double Norm_Euclid(double *a, int n)
{
	double value;

	value = sqrt(Scal(a,a,n));

	return value;	
}
double Projection_On_Axis(double *v,double *o) //  v   o
{
	double value;

	value = Scal(v,o,3)/sqrt(Scal(o,o,3));

	return value;
}
T_Mapping_Vec::T_Mapping_Vec()
{
	this->nver = 0;
	this->xyz = 0;
	this->kuzlov = 0;
	this->kpar = 0;
	ed = NULL;
	edges = NULL;
}
T_Mapping_Vec::T_Mapping_Vec(int (*nver)[14], double (*xyz)[3], int kuzlov, int kpar)
{
	this->nver = nver;
	this->xyz = xyz;
	this->kuzlov = kuzlov;
	this->kpar = kpar;
	ed = NULL;
	edges = NULL;
}
T_Mapping_Vec::~T_Mapping_Vec()
{
	if(ed) {delete [] ed; ed = NULL;}
	if(edges) {delete [] edges; edges = NULL;}
}
