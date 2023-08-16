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
 *  This file contains the headers for extracting smoothing subdomains and obtaining a solution in the elements belonging to them
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin                         
 *  Novosibirsk State Technical University,                 
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia       
 *  vdv_wk@mail.ru                                          
 *  Version 2.0 January 16, 2023                            
*/

#pragma once

#include "AbstractFEM3D.h"
#include "Portret.h"
#include "T_Brick.h"

#include "pardiso.h"

struct Local_Vector
{
	double g[8];
	Local_Vector()
	{
		for(int i=0;i<8;i++)g[i]=0.0;
	}
};

// Class - subdomain for smoothing the solution of a three-dimensional problem in subarea
class Subdomain
{
public:

	int my_ind,ncmp;
	char dir_name[32];

	static int subind;
	static int npls;
	static int ntimes;
	static int ftout;
	static int fdfct;

	int material; //     nvkat ( == -1,     )

	vector<int> renumElemFromOldToNew;
	vector<int> renumElemFromNewToOld;
	vector<int> renumNodeFromNewToOld;
	vector< vector<int> > renumEdgeFromNewToOld;
	int n_elem;
	int n_nodes;
	int n_nodes_c;
	
	vector<double> ValueInCenter;	

	int (*nver)[14];
	double (*xyz)[3];
	double (*xyz_r)[3];
	Portret *p;

	double *di;
	double *gg;
	double *pr;
	double *x;

	double *d;
	double *sg;

	bool llt;

	int Init(int material, int levelNeighbors, vector< vector<int> > &PointresForElem, AbstractFEM3D *TaskCalcMesh);

	void AsmGlobalMatrix();
	void AsmGlobalVector(int ipls,int itime,int icmp);
	void CalcRightPartVect(T_Brick &L);

	Subdomain();
	~Subdomain();

	pardiso_solver prds;

	vector<Local_Vector> vLV;

	int nthreads;
};

