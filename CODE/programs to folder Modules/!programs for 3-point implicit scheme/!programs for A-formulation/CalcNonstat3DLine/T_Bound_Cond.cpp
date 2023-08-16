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
 *  This file contains the code for application the main homogeneous boundary conditions in 3D VFEM
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin                         
 *  Novosibirsk State Technical University,                 
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia       
 *  vdv_wk@mail.ru                                          
 *  Version 2.0 January 16, 2023                            
*/

#include "stdafx.h"
#include "T_Bound_Cond.h"
extern ofstream logfile;
extern void Memory_allocation_error(const char *var, const char *func);
T_Bound_Cond::T_Bound_Cond(int n_nodes, int n_edges, int n_bound_nodes,
						   int *bound_nodes,
							double (*xyz)[3], int (*edges)[2],
							int n_elem, int (*ed)[25], int (*nver)[14])
{
	this->n_nodes = n_nodes;
	this->n_edges = n_edges;
	this->n_bound_nodes = n_bound_nodes;
	this->bound_nodes = bound_nodes;
	this->xyz = xyz;
	this->edges = edges;
	this->n_elem = n_elem;
	this->ed = ed;
	this->nver = nver;

	dirichlet = NULL;
	is_edge_bound = NULL;
	is_node_bound = NULL;

	Make_List_Of_Bound_Edges();

	dirichlet = new double[n_edges];
	if(dirichlet == 0)
		Memory_allocation_error("dirichlet", "T_Bound_Cond::T_Bound_Cond");

	Make_Dirichlet();
}
T_Bound_Cond::~T_Bound_Cond()
{
	if(dirichlet) {delete [] dirichlet; dirichlet=NULL;}
	if(is_edge_bound) {delete [] is_edge_bound; is_edge_bound=NULL;}
	if(is_node_bound) {delete [] is_node_bound; is_node_bound=NULL;}
}
void T_Bound_Cond::Make_List_Of_Bound_Edges()
{
	int i, tmp;

	is_node_bound = new int[n_nodes];
	if(is_node_bound == 0)
		Memory_allocation_error("is_node_bound", "T_Bound_Cond::Make_List_Of_Bound_Edges");


	for(i=0; i<n_nodes; i++)
		is_node_bound[i] = 0;

	for(i=0; i<n_bound_nodes; i++)
		is_node_bound[bound_nodes[i]] = 1;

	is_edge_bound = new int[n_edges];
	if(is_edge_bound == 0)
		Memory_allocation_error("is_edge_bound", "T_Bound_Cond::Make_List_Of_Bound_Edges");

	for(i=0; i<n_edges; i++)
		is_edge_bound[i] = 0;

	tmp = 0;
	for(i=0; i<n_edges; i++)
		if(is_node_bound[edges[i][0]]==1 && is_node_bound[edges[i][1]]==1)		
		{
			is_edge_bound[i] = 1;
			tmp++;
		}

		if(is_node_bound) {delete [] is_node_bound; is_node_bound=NULL;}
}
void T_Bound_Cond::Set_Dirichlet_Cond(int *ig, int *jg, double *di, double *ggl, double *pr,double time, int ntest)
{
	int i, j, k;
	double u_g;

	this->time = time;

	logfile << "Set_Dirichlet_Cond...\n";


	for(i=0; i<n_edges; i++)   
		if(is_edge_bound[i]==1)
		{
			u_g = dirichlet[i];

			di[i] = 1.0; //  -   

			for(k=ig[i]; k<=ig[i+1]-1; k++)
			{
				ggl[k] = 0.0;
			}
		}	

	for(i=0; i<n_edges; i++)
		for(j=ig[i]; j<=ig[i+1]-1; j++)
		{
			k = jg[j];
			if(is_edge_bound[k]==1) 
			{
				ggl[j] = 0.0;
			}
		}	
}
void T_Bound_Cond::Set_Dirichlet_Cond_Pr(double *pr, double time, int ntest,int ipls)
{
	int i, j, k;
	double u_g;

	this->time = time;

	for(i=0; i<n_edges; i++)
	{
		if(is_edge_bound[i]==1)
		{
			u_g = dirichlet[i];
			pr[ipls*n_edges+i] = u_g;
		}
	}
}
void T_Bound_Cond::Make_Dirichlet()
{
	for (int i=0; i<n_edges; i++)
	{
		dirichlet[i] = 0.0;
	}  
}
