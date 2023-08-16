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
 *  This file contains the headers for extracting anomal edges for 3D VFEM task
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin                         
 *  Novosibirsk State Technical University,                 
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia       
 *  vdv_wk@mail.ru                                          
 *  Version 2.0 January 16, 2023                            
*/

#pragma once
#include "T_Mapping.h"
#include "vec_prep_data.h"
class Time_approx_for_vfem
{
public:

	int n;         //   

	T_Mapping_Vec *Tmap;

	int ntime;     //   
	double *time;   //  
	double finSP;

	bool *is_field_in_node;
	int n_nodes_f;
	int *nodes_f;

	int nmat;
	vector<int> ElemsSizeForNodes,ElemsForNodes,ElemsIgNodes;
	vector<int> MatsForNodes,MatsIgNodes,MatsSizeForNodes;
	vector<int> ElemsForMats,ElemsIgMats,ElemsSizeForMats;
	vector<double> RegularLocalCoord;

	int *anomal_edges; //   
	int n_anomal_edges; //   

	Vec_Prep_Data *d;

	Time_approx_for_vfem(bool forPolygon, const void* _NormBFromPolygon, 
		bool flag_vfem_for_line, const void* _NormResFromLine);
	~Time_approx_for_vfem();

	int Read_data();

	int Unload_anomal_nodes(char *fname);
	int Unload_anomal_nodesA0(char *fname);

	vector<int> elemRenum;
	vector<bool> isElemAnomal;
	int nAnomalElem;

	vector<int> isEdgesAnomal;
	vector<double> EdgesCrd[3];

	int nzl;
	vector<double> zlay;
};//------------------------------------------------------------------------
