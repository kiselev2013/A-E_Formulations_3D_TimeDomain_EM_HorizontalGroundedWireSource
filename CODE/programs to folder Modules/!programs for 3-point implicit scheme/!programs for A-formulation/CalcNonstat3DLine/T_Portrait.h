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
 *  This file contains the headers for building portrait of SLAE matrix for 3D VFEM
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin                         
 *  Novosibirsk State Technical University,                 
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia       
 *  vdv_wk@mail.ru                                          
 *  Version 2.0 January 16, 2023                            
*/

#pragma once

// The class contains procedures for constructing a portrait of the SLAE matrix in a 3D VFEM
class T_Portrait
{
public:
	int (*ed)[25]; // ,   
	int n_elem; //    
	int n; //   
	int n_c; //   

	int *ig; 
	int *jg;
	int size_jg;

	int *idi; //    
	int *ijg; //    

	T_Portrait(int *ed, int n, int n_c, int n_elem);
	~T_Portrait();

	void Gen_Portrait();


	void Gen_idi_ijg(int *nvkat, int (*nver)[14]);

	void Set_type_of_block(int *target_array, int adr, int type);
};
const int FILTER_MASS_MATRIX_VEC[12][12] = { 
   2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1,
   2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1,
   2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1,
   2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1,
   1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 
   1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 
   1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 
   1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2,  
   1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2,  
   1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2,  
   1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2
};
