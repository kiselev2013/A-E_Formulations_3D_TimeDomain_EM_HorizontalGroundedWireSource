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
 *  This file contains the headers for assembling SLAE in 3D VFEM
 *  
 *  
 *  Written by Prof. Yuri G. Soloveichik, Prof. Marina G. Persova, Ph.D. Denis V. Vagin   
 *  Novosibirsk State Technical University,                                               
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                     
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                      
 *  Version 2                                                                             
*/

#pragma once
#include "vec_prep_data.h"
#include "T_Mapping.h"
#include "T_Portrait.h"
#include "T_Brick.h"

// The class contains procedures for calculating local matrices and 
// assembling a global matrix for a non-stationary 3D VFEM problem
class T_Global_SLAE_Vec
{
public:
	int *ig;    //    
	int *jg;    //    
	double *di;  //   
	double *ggl; //    ( )
	double *ggu; //    ( )
	double *pr;  //    

	Vec_Prep_Data *d;

	int *ig_t; 
	int *jg_t;
	double *gg_t;
	
	int n_elem;    //    
	int n_edges;   //   ()
	int n_edges_c; //   ()
	int n;         //  
	int ig_n_1;    // ig(n+1)-1    

	int (*nver)[14];
	int (*ed)[25];
	double (*xyz)[3];
	int (*edges)[2];

	int n_of_materials;
	int *nvkat;
	double *sigma3d;
	double *sigma0;
	double *mu3d;

	double *En; //    

	T_Global_SLAE_Vec(Vec_Prep_Data *d, T_Mapping_Vec *T, T_Portrait *P,int npls);

	~T_Global_SLAE_Vec();

	void Assembling(const int what_compute, double *En);

	void Assembling_pr(double *En,int ipls);

	void AsmPrMuEpsSigma(double *An, double *En, double *d2An, int ipls);
};
