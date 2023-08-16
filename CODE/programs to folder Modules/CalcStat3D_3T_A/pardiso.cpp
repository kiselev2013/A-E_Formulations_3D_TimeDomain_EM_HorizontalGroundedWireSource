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
 *  This file contains the code for using MKL PARDISO solver
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin                         
 *  Novosibirsk State Technical University,                 
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia       
 *  vdv_wk@mail.ru                                          
 *  Version 2.0 January 16, 2023                            
*/

#include "stdafx.h"
#include "pardiso.h"

extern ofstream logfile;

pardiso_solver::pardiso_solver()
{
	ia=NULL;
	ja=NULL;
	a=NULL;
	perm=NULL;

	n = 0;

	init();
}

pardiso_solver::~pardiso_solver()
{
	clear();
}

void pardiso_solver::init()
{
	int i;

	mtype = 2;
	nrhs = 1;
	maxfct = 1;
	mnum = 1;
	msglvl = 1;
	phase = 13;

	for(i=0;i<64;i++){pt[i]=0;}
	for(i=0;i<64;i++){iparm[i]=0;}


	msglvl=0;
}

void pardiso_solver::clear()
{
	if(a){delete [] a;  a=NULL;}
	if(ia){delete [] ia; ia=NULL;}
	if(ja){delete [] ja; ja=NULL;}
	if(perm){delete [] perm; perm=NULL;}
}

void pardiso_solver::factorize(int nb,int *ig,int *jg,double *ggl,double *di,int nthreads)
{
	int i;
	int ig_n_1=0;
	int sz_iptr=0;
	int sz_jptr=0;

	init();

	if(n!=nb)
	{
		clear();
	}

	mkl_set_num_threads(nthreads);

	f.FromRSFToCSR_Real_1_Sym(nb, ig, &sz_iptr, &sz_jptr);

	if(n!=nb)
	{
		ia = new MKL_INT64[sz_iptr];
		ja = new MKL_INT64[sz_jptr];
		a = new double[sz_jptr];
	}

	f.FromRSFToCSR_Real_2_Sym(nb, ig, jg, di, ggl, ia, ja, a);

	for(i=0;i<sz_iptr;i++){ia[i]++;}
	for(i=0;i<sz_jptr;i++){ja[i]++;}

	if(n!=nb)
	{
		perm = new MKL_INT64[nb];
	}

	phase = 12;
	n = nb;

	pardiso_64(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, NULL, NULL, &info);

	cout<<"info= "<<info<<'\n';

	if(info)
	{
		ofstream ofp;
		ofp.open("fitting", ios::app);
		ofp<<"Error spline info= "<<info<<'\n';
		ofp.close();
		ofp.clear();
		exit(1);
	}
}

void pardiso_solver::solve_nrhs(int _nrhs,double *pr,double *x)
{
	int i;
	for(i=0;i<n;i++){x[i]=0.0;}
	phase = 33;
	nrhs = _nrhs;
	pardiso_64(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, pr, x, &info);
}

void pardiso_solver::stop_solver()
{
	phase = -1;
	pardiso_64(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, NULL, NULL, &info);
}
