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
 *  This file contains the headers for calculating the values of basic functions and their derivatives on an element
 *  
 *  
 *  Written by Prof. Yuri G. Soloveichik, Prof. Marina G. Persova                              
 *  Novosibirsk State Technical University,                                                                               
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                                                     
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)        
 *  Version 2.0 January 16, 2023                                                                                          
*/

#pragma once
#include "Vec_Prep_Data.h"

// The class contains procedures for calculating the value of local basis functions for a 3D VFEM problem
class T_Brick
{
public:
	double hx, hy, hz; //  -
	double xk, xk1, yk, yk1, zk, zk1; //  
    int num; //     

	double b[12][12]; //   
	double c[12][12]; //   
	double f_re[12];  //      

	double a[12][12]; //   -  
	double g[12];     //   -
	double g8[8]; //  

	double mu;
	double mu0;
	double sigma;
	double sigma0;
	double dpr;
	double dpr0;
	int n_mat; //   

	int (*nver)[14]; //     (13- )
	int (*ed)[25];   // -    +   +  -
	int *nvkat;      //    -
	int (*edges)[2]; // ,  2- 
	double (*xyz)[3]; //  

	double *En; //   ( )

	int tasktype;

	double mrot[3][3];

	T_Brick(double *x_coords, double *y_coords, double *z_coords);

	T_Brick(int num, int (*nver)[14], int (*ed)[25], int (*edges)[2], double (*xyz)[3], int *nvkat,
		double *sigma3d, double *sigma0, double *mu3d, double *En);

	T_Brick(int num, int (*nver)[14], double (*xyz)[3]);

	T_Brick();

	~T_Brick();

	void Compute_Local_Matrix_And_Vector(const int what_compute); 

	void Compute_Local_Matrix_B(); 

	void Compute_Local_Matrix_C(); 

	void Compute_Local_Vector_For_Anomal_Problem();

	void ComputeLocalVectorMuEpsSigma(double *An, double *d2An);

	void Calc_J_Node_2(double x, double y, double z);
	double dPhi_node_2(int i, int j, double x, double y, double z);
	void Calc_V_Node_2(double *q,double x, double y, double z);
	double V[3];


	double x[8], y[8], z[8]; //   	
	double J[3][3];          //  
	double J_1[3][3];		 // J^{-1}
	double J_1_T[3][3];		 // J^{-T}
	double det_J;            //  (  )
	double det_J_abs;        //  
	double phi_all[12][3];    //     
	double rot_all[12][3];    //   .    

	void Calc_J(int n_of_point); //       
	void Calc_J_on_face(int n_of_point); //    ,   
	void Calc_J_in_parallelepiped(); //      

	void Calc_local_matrix_b_for_hexahedron(); //   ()
	void Calc_local_matrix_c_for_hexahedron(); //   ()

	void Calc_value_inside_hex(double *ves, double *in, double *out); 

	void Basis_func_on_vec_par(int i, double *in, double *out);

	void Basis_func_on_vec_par(int i, double ves, double *in, double *out);

	void Basis_func_on_reference_vec_par(int i, double *in, double *out);	

	void Basis_func_on_vec_hex(int i, double ves, double *in, double *out);

	void Calc_rotor_inside_hex(double *ves, double *in, double *out); 

	void Calc_rotor_inside_hex_crd(double *ves, double *in, double *out,int ind,
		/*double p_J[3][3],double p_J_1_T[3][3],double p_det_J_abs,*/ double *cfbx, double *cfby, double *cfbz); 

	void Calc_rotor_coefficients_inside_hex_crd(double *in, int ind, double *cfbx, double *cfby, double *cfbz);

	void Rot_of_basis_func_on_reference_vec_par(int i, double *in, double *out);

	void Rotx_of_basis_func_on_reference_vec_par(int i, double x, double *out);

	void Roty_of_basis_func_on_reference_vec_par(int i, double y, double *out);

	void Rotz_of_basis_func_on_reference_vec_par(int i, double z, double *out);

	void Rot_of_basis_func_on_vec_par(int i, double ves, double *in, double *out);

	void Rot_of_basis_func_on_vec_par_x(int i, double ves, double *in, double &out);
	void Rot_of_basis_func_on_vec_par_y(int i, double ves, double *in, double &out);
	void Rot_of_basis_func_on_vec_par_z(int i, double ves, double *in, double &out);

	void Rot_of_basis_func_on_vec_par_x(int i, double *in, double &out);
	void Rot_of_basis_func_on_vec_par_y(int i, double *in, double &out);
	void Rot_of_basis_func_on_vec_par_z(int i, double *in, double &out);

	void Get_rotz_on_face(double *ves1, double *ves2, double *ves3,
		double *out1, double *out2, double *out3);

	void Get_rot_on_face(double *ves1, double *ves2, double *ves3,
						  double *out1, double *out2, double *out3,
						  double *outx1, double *outx2, double *outx3,
						  double *outy1, double *outy2, double *outy3);

	double l0(double x);
	double l1(double x);

	double Phi_node(int i, double x, double y, double z);

	double dPhi_node(int i, int j, double x, double y, double z);

	void VectorFieldOnPar(double x, double y, double z, double *ves,
		double *x_out, double *y_out, double *z_out);
	double ScalarFieldOnPar(double x, double y, double z, double *ves);
	double DxOfScalarFieldOnPar(double x, double y, double z, double *ves);
	double DyOfScalarFieldOnPar(double x, double y, double z, double *ves);
	double DzOfScalarFieldOnPar(double x, double y, double z, double *ves);

	void ScalarFieldOnParCff(double x, double y, double z, double *cff);

	double GetValueInHexCenter(double *q);
	void GetGradInHexCenter(double *q, double *out, int cj);

	void VectorFieldXOnPar3(double y, double z, double *ves_j2, double *ves_j1, double *ves_j,
		double *out_j2, double *out_j1, double *out_j);

	void VectorFieldYOnPar3(double x, double z, double *ves_j2, double *ves_j1, double *ves_j,
		double *out_j2, double *out_j1, double *out_j);

	void RotXOnPar(double x, double *ves, double *out, bool loc_c=false);
	void RotYOnPar(double y, double *ves, double *out, bool loc_c=false);
	void RotZOnPar(double z, double *ves, double *out);
	void RotZOnPar3(double z,
		double *ves1, double *ves2, double *ves3,
		double *out1, double *out2, double *out3);

	void Transformation_of_variables(const double *in, double *out);
	void Transformation_of_variables(double *x, double *y, double *z);
	double Xi(double x);
	double Eta(double y);
	double Zeta(double z);



	void GetAonFace(double *ves1, double *ves2, double *ves3,
		double *out1, double *out2, double *out3);

	void Set_dpr(double dpr);
	void Set_dpr0(double dpr0);
	void Set_mu0(double mu0);
	Vec_Prep_Data *d;
	

	double asin0n[8][3], acos0n[8][3]; //    
	double asin0c[3], acos0c[3]; //     

	void LocalVectHarmConst();

	void GetVectorFieldNodes(double *ves, double *ax, double *ay, double *az);

};
const double MIDDLE_OF_LOCAL_EDGE[12][3] = {
	 0.0, -1.0, -1.0,
	 0.0,  1.0, -1.0,
	 0.0, -1.0,  1.0,
	 0.0,  1.0,  1.0,

	-1.0,  0.0, -1.0,
	-1.0,  0.0,  1.0,
	 1.0,  0.0, -1.0,
	 1.0,  0.0,  1.0,

	-1.0, -1.0,  0.0,
	 1.0, -1.0,  0.0,
	-1.0,  1.0,  0.0,
	 1.0,  1.0,  0.0
};
const double LOCAL_COORDS_OF_NODES[8][3] = 
{
	-1.0, -1.0, -1.0,
	 1.0, -1.0, -1.0,
	-1.0,  1.0, -1.0,
	 1.0,  1.0, -1.0,

	-1.0, -1.0,  1.0,
	 1.0, -1.0,  1.0,
	-1.0,  1.0,  1.0,
	 1.0,  1.0,  1.0
};
const double TANGENT_VECTORS_ON_REFERENCE_CUBE[12][3] = {
	1.0, 0.0, 0.0,
	1.0, 0.0, 0.0,
	1.0, 0.0, 0.0,
	1.0, 0.0, 0.0,

	0.0, 1.0, 0.0,
	0.0, 1.0, 0.0,
	0.0, 1.0, 0.0,
	0.0, 1.0, 0.0,

	0.0, 0.0, 1.0,
	0.0, 0.0, 1.0,
	0.0, 0.0, 1.0,
	0.0, 0.0, 1.0
};
const int REG_EDGES[12][2]={ 	//      
	0,1, 2,3, 4,5, 6,7, 0,2, 4,6, 1,3, 5,7, 0,4, 1,5, 2,6, 3,7 };
