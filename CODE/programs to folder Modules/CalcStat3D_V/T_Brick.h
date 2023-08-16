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

// The class contains procedures for calculating the value of local basis functions for a 3D VFEM problem
class T_Brick
{
public:
	double hx, hy, hz; //  -
	double xk, xk1, yk, yk1, zk, zk1; //  
    int num; //     

	double g8[8]; //  

	double mu;
	double mu0;
	double sigma;
	double sigma0;
	double dpr;
	double dpr0;
	int n_mat; //   

	int (*nver)[14]; //     (13- )
	int *nvkat;      //    -
	double (*xyz)[3]; //  


	T_Brick(double *x_coords, double *y_coords, double *z_coords);
	T_Brick(int num, int (*nver)[14], double (*xyz)[3]);

	~T_Brick();

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

	void Mapping(double *in, double *out);

	void Calc_J(int n_of_point); //       
	void Calc_J(double x, double y, double z); //       
	void Calc_J_on_face(int n_of_point); //    ,   
	void Calc_J_in_parallelepiped(); //      

	void Calc_value_inside_hex(double *ves, double *in, double *out); 

	double l0(double x);
	double l1(double x);

	double Phi_node(int i, double x, double y, double z);

	double dPhi_node(int i, int j, double x, double y, double z);

	double ScalarFieldOnPar(double x, double y, double z, double *ves);
	double DxOfScalarFieldOnPar(double x, double y, double z, double *ves);
	double DyOfScalarFieldOnPar(double x, double y, double z, double *ves);
	double DzOfScalarFieldOnPar(double x, double y, double z, double *ves);

	double GetValueInHexCenter(double *q);
	void GetGradInHexCenter(double *q, double *out, int cj);

	void Transformation_of_variables(const double *in, double *out);
	void Transformation_of_variables(double *x, double *y, double *z);
	double Xi(double x);
	double Eta(double y);
	double Zeta(double z);
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
