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
 *  This file contains the headers for working with edge mesh in 3D VFEM
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin                         
 *  Novosibirsk State Technical University,                 
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia       
 *  vdv_wk@mail.ru                                          
 *  Version 2.0 January 16, 2023                            
*/

#pragma once

// The class contains a 3D edge mesh
class T_Mapping_Vec
{
public: 
	T_Mapping_Vec();
	T_Mapping_Vec(int (*nver)[14], double (*xyz)[3], int kuzlov, int kpar);
	~T_Mapping_Vec();

	int n_c;  //    (continuous)
	int n_dc; //    (discontinuous)
	int n;    //   n = n_c + n_dc

	int kuzlov; // -   
	int kpar;   // - 
	int (*nver)[14]; //     +   +  -
	double (*xyz)[3]; //   ( )

	int (*edges)[2]; // ,  2- 
	int (*ed)[25];   //     +   +  -  !!!!!   -   
};
