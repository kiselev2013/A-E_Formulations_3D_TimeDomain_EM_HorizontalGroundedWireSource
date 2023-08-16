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
 *  This file contains the headers for using MKL PARDISO solver
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin                         
 *  Novosibirsk State Technical University,                 
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia       
 *  vdv_wk@mail.ru                                          
 *  Version 2.0 January 16, 2023                            
*/

#pragma once
#include "FormatConverter.h"

// Class - adapter to use MKL PARDISO solver
struct pardiso_solver
{
	FormatConverter f;

	MKL_INT64 *ia;
	MKL_INT64 *ja;
	double *a;

	MKL_INT64 n;
	MKL_INT64 mtype; // real and symmetric positive definite
	MKL_INT64 nrhs;
	void *pt[64];
	MKL_INT64 maxfct;
	MKL_INT64 mnum;
	MKL_INT64 msglvl;
	MKL_INT64 phase;
	MKL_INT64 *perm;
	MKL_INT64 iparm[64];
	MKL_INT64 info;

	pardiso_solver();
	~pardiso_solver();

	void factorize(int nb,int *ig,int *jg,double *ggl,double *di,int nthreads);
	void solve_nrhs(int nrhs,double *pr,double *x);
	void stop_solver();
	void clear();
	void init();
};
