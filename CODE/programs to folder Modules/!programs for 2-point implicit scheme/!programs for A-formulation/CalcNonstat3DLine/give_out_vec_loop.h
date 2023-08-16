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
 *  This file contains the headers for outputting E and EMF fields from solution in 3D VFEM
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin                         
 *  Novosibirsk State Technical University,                 
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia       
 *  vdv_wk@mail.ru                                          
 *  Version 2.0 January 16, 2023                            
*/

#pragma once
#include "vec_prep_data.h"
#include "T_Mapping.h"
#include "TaskSP.h"
#include "OutputResultant3d.h"
class Give_out_vec_loop 
{
public:

	Give_out_vec_loop(Vec_Prep_Data *d, T_Mapping_Vec *tmap,
		vector<int> &_RecvPlsIgB,vector<int> &_RecvPlsIgE, int _npls);

	~Give_out_vec_loop();

	OutputResultant3d *resultantA, *resultantB;

	bool forLine;
	
	Vec_Prep_Data *d;					//  (     )
	
	T_Mapping_Vec *tmap;				// T-    
	
	const double *v3_j, *v3_j1, *v3_j2; //     
	double dt, dt0, dt1;				//     3- 
	double t_j, t_j1, t_j2;				//   
	
	int tnum;							//     ( -  )
	int tnum0, tnum1, tnum2;			//       


	double time_interval_for_print_left; 
	double time_interval_for_print_right; 


	int time_layer_first;
	int time_layer_last;



	double *Edsz_all, *Edsx_all, *Edsy_all;
	double *Bz_all, *Bx_all, *By_all;

	double *Edszn_all, *Edsxn_all, *Edsyn_all;
	double *Bzn_all, *Bxn_all, *Byn_all;

	double *Ex_all, *Ey_all, *Ez_all; 

	double *Exn_all, *Eyn_all, *Ezn_all; 



	void Work2(double *u_j2, double *u_j1, double *u_j,
		int t_2, int t_1, int t_0, int t);

	void Work_on_current_time_layer();

	void Gather_edsall(bool forCED);

	double dA_dt(double t,
		double u_j, double u_j1, double u_j2, 
		double dt, double dt0, double dt1, 
		double t_j, double t_j1, double t_j2);


	void Clear(double *u, const int& usize);

	vector<int> elemForPointB; //     -,   
	vector<int> elemForPointE; //     -,   
	int PrepareForParallelepipedOutput();

	int npls,ipls_cur;
	vector<int> vPre1,vPre2;
	bool fsdiff;

	int cur_dec_strt_res;

	vector<int> RecvPlsIgB,RecvPlsIgE;

	vector<Res3DValueType> vvta,vvtb;
};
