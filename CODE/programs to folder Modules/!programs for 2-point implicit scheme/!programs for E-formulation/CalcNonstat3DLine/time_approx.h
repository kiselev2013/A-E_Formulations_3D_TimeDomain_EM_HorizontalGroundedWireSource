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
 *  This file contains the headers for calculating nonstationary 3D VFEM task
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
#include "T_Portrait.h"
#include "T_Mapping.h"
#include "T_Global_SLAE.h"
#include "T_Bound_Cond.h"
#include "Give_out_vec_loop.h"
#include "vec_prep_data.h"
#include "pardiso.h"

// The class contains procedures for solving a non-stationary 3D VFEM problem using a three-layer implicit scheme in time
class Time_approx_for_vfem  : public AbstractFEM3D
{
public:
	int n;         //   

	double *u_j;    //   j- () 
	double *u_j1;   //   (j-1)-  (  )
	double *u_j2;   //   (j-2)- 
	double *u_r;    //   
	double *cuj;    // C*U^{j-1}  C*U^{j-2}
	double *cujsum;
	double *y_omp;

	Give_out_vec_loop *give_out;

	T_Mapping_Vec *Tmap;

	T_Portrait *Pt;

	T_Global_SLAE_Vec *Slae;

	double *b_gg_all, *c_gg_all, *b_di_all, *c_di_all, *c_gg_all_dsig, *c_di_all_dsig;

	T_Bound_Cond *Bc;

	int n_max_iter_per_layer; //    
	double eps; // epsilon  

	double *help;
	double *help1;
	double *help2;
	double *help3;
	double *help4;
	double *help5;
	double *help6;
	double *help7;
	double *help8;

	int ntime;     //   
	double *time;   //  
	double dt;      //  2-: t_{j}-t_{j-1};  3-: t_{j}- t_{j-2}.
	double dt0;     // t_{j}   -  t_{j-1}
	double dt1;     // t_{j-1} -  t_{j-2}
	double dt2;
	double dt3;
	double dt4;
	double finSP;

	bool *is_field_in_node;
	int n_nodes_f;
	int *nodes_f;

	double (*En_nodes)[3]; // dA/dt   (Ez=0   ) 
	double *En_edges; // dA/dt  
	double *En_edges_0; // dA/dt  
	double *En_edges_1; // dA/dt  
	double *En_edges_2; // dA/dt  
	double *An_edges;
	double *d2An_edges; // d2A/dt2  

	int nmat;
	vector<int> ElemntsForRegularEdges,RegEdgesCoords;
	vector<double> RegularLocalCoord;

	int *anomal_edges; //   
	int n_anomal_edges; //   
	vector<bool> isEdgesAnomal;
	int npntE0;

	int AnomalType,ia;

	Vec_Prep_Data *d;

	Time_approx_for_vfem(int _npls);
	~Time_approx_for_vfem();

	int Read_data(vector<int> &RecvPlsIgB,vector<int> &RecvPlsIgE);

	int Prepare_To_Schema();

    int Two_Layers_Schema(int nlst);

	int Three_Layers_Schema(int tbeg,int tend,int i_u_1,int i_u_2,int nlst,int tbeg_pre,bool f_last_dec);

	void FinishCalculation();

	int Projection_from_nodes_to_edges(int n_edges_c, 
		double (*En_nodes)[3], double *En_edges, double (*xyz)[3],
		int (*edges)[2], bool *is_field_in_node, int n_elem, int (*nver)[14],
		int (*ed)[25]);

	double Calc_dof(double *J, double *func, int n_local_edge);

	void GetE0FromRZ(int tnum, double (*En_nodes)[3]);

	int GetNumberOfNodes();
	int GetNumberOfElements();
	int GetElementNodesNumber();
	const pv::Point3D GetNode(const int& i_node);
	const pv::Point3D GetNodeTrue(const int& i_node);
	int GetNodeNumberOnElement(const int& i_element, const int& i_node);
	int GetElementMaterial(const int& i_element);
	int GetTypeOfElement(const int& i_element);
	double GetValueInElementCenter(const int& i_element, const Res3DValueType& r_type);
	int GetNumberOfResPoints(const Res3DValueType& r_type);
	pv::Point3D GetResPoint(const Res3DValueType& r_type, const int& i_point);
	int * GetPointerToRegular();
	int GetXSize();
	int GetYSize();
	int GetZSize();
	double *GetPointerToX();
	double *GetPointerToY();
	double *GetPointerToZ();
	void SaveResult(const Res3DValueType& r_type, const double& r_value, const int& j, const int& tnum, int ipls);


	bool iVP;
	bool forLine;
	bool forCED;

	int ReadV3(int nlayer, double *ex, double *ey, double *ez);
	vector<int> elemRenum;
	vector<bool> isElemAnomal;
	int nAnomalElem;

	int npls,ipls_cur,nthreads;

	void SaveCurrentResults(int iDec,int tbeg,int tend,vector<int> &RecToSourceB,vector<int> &RecToSourceE);
	int LoadPreviousResults(int iDec,int tbeg,int tend,vector<int> &RecToSourceB,vector<int> &RecToSourceE);
	void WriteTemporaryAllFile(int nrec,vector<int> &RecToSource,int tbeg,int tend,double *FieldNorm,double *FieldAnom,char *fname);
	int ReadTemporaryAllFile(int nrec,vector<int> &RecToSource,int tbeg,int tend,double *FieldNorm,double *FieldAnom,char *fname);

	int InitvPre(int tbeg,int tend,int i_u_1,int i_u_2);

	int cur_dec_strt,max_dec_size;

	vector<int> RecvPlsIgB,RecvPlsIgE;

	pardiso_solver prds;

	int LoadVectorE0();

	void DiffAnorm(vector<int> &RecToSourceE);

	char *Buff,*Buff_A;

	double eps_reg;

	double dA_dt(double t,
		double u_j, double u_j1, double u_j2, 
		double dt, double dt0, double dt1, 
		double t_j, double t_j1, double t_j2);

	void _dA_dt_vec(double *du_j,double *u_j, double *u_j1, double *u_j2, 
		double t_j,double t_j1,double t_j2,int size,double t);

	int *anomal_elem_mat; //   

	int *normal_elem_for_anomal_edges[3][8]; // 
	int *mat_elem_for_anomal_edges[8]; // 

	int OpenEnorm(int nlayer,ifstream &fp);

	int ForInv;

	int LoadE0FromRZ(int tnum,double *En_edges);
	int LoadE0DiffFromRZ(int tnum,int tnum1,int tnum2);
};//------------------------------------------------------------------------
