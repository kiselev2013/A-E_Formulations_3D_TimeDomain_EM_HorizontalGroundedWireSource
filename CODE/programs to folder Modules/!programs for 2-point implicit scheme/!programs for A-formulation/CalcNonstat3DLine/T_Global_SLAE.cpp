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
 *  This file contains the code for assembling SLAE in 3D VFEM
 *  
 *  
 *  Written by Prof. Yuri G. Soloveichik, Prof. Marina G. Persova, Ph.D. Denis V. Vagin 
 *  Novosibirsk State Technical University,                                             
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                   
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                    
 *  Version 2                                                                           
*/

#include "stdafx.h" 
#include "T_Global_SLAE.h"
#include "in_out.h"
extern ofstream logfile;
extern void Memory_allocation_error(const char *var, const char *func);
T_Global_SLAE_Vec::T_Global_SLAE_Vec(Vec_Prep_Data *d, T_Mapping_Vec *T, T_Portrait *P,int npls)
{
	this->d = d;

	this->ig = P->ig;
	this->jg = P->jg;

	this->n_elem = d->kpar;
	this->n_edges = T->n; 
	this->n_edges_c = T->n_c; 

	this->n = this->n_edges_c; //    
	this->ig_n_1 = this->ig[n];

	this->nver = d->nver;
	this->ed = T->ed;
	this->edges = T->edges;
	this->xyz = d->xyz;

	this->n_of_materials = d->n_materials;

	this->sigma3d = d->sigma3d;
	this->sigma0 = d->sigma0;
	this->mu3d = d->mu3d;
	this->nvkat = d->nvkat;

	di = NULL;
	ggl= NULL;
	ggu= NULL;
	pr = NULL;

	di = new double[n];
	if(di == 0)
		Memory_allocation_error("di", "T_Global_SLAE_Vec::T_Global_SLAE_Vec");

	pr = new double[n*npls];
	if(pr == 0)
		Memory_allocation_error("pr", "T_Global_SLAE_Vec::T_Global_SLAE_Vec");

	ggl = new double[ig_n_1];
	if(ggl == NULL)
		Memory_allocation_error("ggl", "T_Global_SLAE_Vec::T_Global_SLAE_Vec");
}
T_Global_SLAE_Vec::~T_Global_SLAE_Vec()
{
	if(di) {delete [] di; di=NULL;}
	if(pr) {delete [] pr; pr=NULL;}
	if(ggl) {delete [] ggl; ggl=NULL;}
	if(ggu) {delete [] ggu; ggu=NULL;}
}
void T_Global_SLAE_Vec::Assembling(const int what_compute, double *En)
{
	/*
	what_compute
	0 - B and C and pr
	1 - C and pr
	2 - add B
	3 - C_eps   
	*/
	int i, j, k, m, it, jt, i_mu, j_nu;
	int ii, jj; //  

	logfile << "Assembling global SLAE using T-mapping...\n";
	
	switch(what_compute)
	{
	case 0: // B and C and pr
	case 1: // C and pr
	case 3: // C_eps   
	case 4: // C_sig and pr
	case 5: // B_mu and pr
		for(i=0; i<n; i++)	pr[i] = di[i] = 0.0;
		for(i=0; i<ig_n_1; i++) ggl[i] = 0.0;
		break;
	}

	for(i=0; i<this->n_elem; i++)
	{
		T_Brick L(i, nver, ed, edges, xyz, nvkat, d->sigma3d, sigma0, mu3d, En);

		L.Set_dpr(d->dpr3d[d->nvkat[i]]);
		L.Set_dpr0(d->dpr0[d->nvkat[i]]);

		L.Compute_Local_Matrix_And_Vector(what_compute);

		for(j=0; j<12; j++)
		{
			ii = ed[i][j];

			if(ii >= n_edges_c)
			{//   
				for(it = ig_t[ii]; it<=ig_t[ii+1]-1; it++)
				{
					i_mu = jg_t[it];
					pr[i_mu] += L.g[j]*gg_t[it]; 
				}
			}
			else
			{   //   
				pr[ii] += L.g[j];
				di[ii] += L.a[j][j];
			}

			for(k=0; k<12; k++)
			{
				jj = ed[i][k];

				if(ii < n_edges_c && jj < n_edges_c) 
				{
					if(jj < ii) 
						for(m=ig[ii]; m<=ig[ii+1]-1; m++)
							if(jg[m]==jj)
							{
								ggl[m] += L.a[j][k];
								break;
							}
				}

				else if(ii >= n_edges_c && jj < n_edges_c)
				{
					for(it = ig_t[ii]; it<=ig_t[ii+1]-1; it++)
					{
						i_mu = jg_t[it];

						if(jj < i_mu) //   
						{
							for(m=ig[i_mu]; m<=ig[i_mu+1]-1; m++)
								if(jg[m]==jj)
								{
									ggl[m] += L.a[j][k]*gg_t[it];
									break;
								}
						}
						else if(jj == i_mu)
						{
							di[i_mu] += L.a[j][k]*gg_t[it];
						}
					}
				}
				else if(ii < n_edges_c && jj >= n_edges_c)
				{
					for(it = ig_t[jj]; it<=ig_t[jj+1]-1; it++)
					{
						j_nu = jg_t[it];

						if(j_nu < ii) //   
						{
							for(m=ig[ii]; m<=ig[ii+1]-1; m++)
								if(jg[m]==j_nu)
								{
									ggl[m] += L.a[j][k]*gg_t[it];
									break;
								}
						}
						else if(j_nu == ii)
						{
							di[j_nu] += L.a[j][k]*gg_t[it];
						}
					}
				}
				else if(ii >= n_edges_c && jj >= n_edges_c)
				{
					for(it = ig_t[ii]; it<=ig_t[ii+1]-1; it++)
					{
						i_mu = jg_t[it];
						
						for(jt = ig_t[jj]; jt<=ig_t[jj+1]-1; jt++)
						{
							j_nu = jg_t[jt];

							if(j_nu < i_mu) //   
							{
								for(m=ig[i_mu]; m<=ig[i_mu+1]-1; m++)
									if(jg[m]==j_nu)
									{
										ggl[m] += L.a[j][k]*gg_t[it]*gg_t[jt];
										break;
									}
							}
							else if(j_nu == i_mu)
							{
								di[i_mu] += L.a[j][k]*gg_t[it]*gg_t[jt];
							}
						}//jt
					}//it
				}//else
			}// k
		}// j
	}// i
}
void T_Global_SLAE_Vec::Assembling_pr(double *En,int ipls)
{
	int i, j, it, i_mu;
	int ii; //  

	logfile << "Assembling global vector using T-mapping...\n";

	for(i=0;i<n;i++){pr[ipls*n+i]=0.0;}

	for(i=0; i<this->n_elem; i++)
	{
		T_Brick L(i, nver, ed, edges, xyz, nvkat, sigma3d, sigma0, mu3d, En);

		L.Calc_local_matrix_c_for_hexahedron();
		L.Compute_Local_Vector_For_Anomal_Problem();

		for(j=0; j<12; j++)
		{
			ii = ed[i][j];

			if(ii >= n_edges_c)
			{//   
				for(it = ig_t[ii]; it<=ig_t[ii+1]-1; it++)
				{
					i_mu = jg_t[it];
					pr[ipls*n_edges_c+i_mu] += L.g[j]*gg_t[it]; 
				}
			}
			else
			{   //   
				pr[ipls*n_edges_c+ii] += L.g[j];
			}
		}
	}// i
}
void T_Global_SLAE_Vec::AsmPrMuEpsSigma(double *An, double *En, double *d2An, int ipls)
{
	int i, j, it, i_mu;
	int ii; //  

	logfile << "Assembling global vector using T-mapping (T_Global_SLAE_Vec::AsmPrMuEpsSigma) ...\n";

	for(i=0; i<n; i++){pr[ipls*n+i]=0.0;}

	for(i=0; i<this->n_elem; i++)
	{
		T_Brick L(i,nver,ed,edges,xyz,nvkat,
			sigma3d,sigma0,mu3d, En);

		L.Set_dpr(d->dpr3d[d->nvkat[i]]);
		L.Set_dpr0(d->dpr0[d->nvkat[i]]);
		L.Set_mu0(d->mu0[d->nvkat[i]]);

		L.Calc_local_matrix_c_for_hexahedron();
		L.Calc_local_matrix_b_for_hexahedron();
		L.ComputeLocalVectorMuEpsSigma(An, d2An);

		for(j=0; j<12; j++)
		{
			ii = ed[i][j];

			if(ii >= n_edges_c)
			{//   
				for(it = ig_t[ii]; it<=ig_t[ii+1]-1; it++)
				{
					i_mu = jg_t[it];
					pr[ipls*n_edges_c+i_mu] += L.g[j]*gg_t[it]; 
				}
			}
			else
			{   //   
				pr[ipls*n_edges_c+ii] += L.g[j];
			}
		}
	}// i
}
