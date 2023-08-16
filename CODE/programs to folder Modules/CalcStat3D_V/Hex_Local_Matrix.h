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
 *  This file contains headers for calculating local matrices on parallelepipeds
 *  
 *  
 *  Written by Prof. Yuri G. Soloveichik, Prof. Marina G. Persova, Ph.D. Denis V. Vagin 
 *  Novosibirsk State Technical University,                                             
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                   
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                    
 *  Version 2                                                                         
*/

#pragma once

// The class contains the data of a 3D grid element
class Hex_Local_Matrix
{
public:
	int number_of_element;
	double x[8], y[8], z[8]; //   	
	double J[3][3];          //  
	double J_1_T[3][3];		 // J^{-T}
	double det_J;            //  (  )
	double det_J_abs;        //  
	double grad_all[8][3]; //   .    
	double b[8][8]; //  
	double hx, hy, hz; //  
	bool JforParCalc; //       
	void Calc_J(int n_of_point);
	void Calc_local_matrix_b_for_parallelepiped();
	void CalcMassMatrix();
	Hex_Local_Matrix(int i, int (*nver)[14], double (*xyz)[3]);
	~Hex_Local_Matrix();
};
