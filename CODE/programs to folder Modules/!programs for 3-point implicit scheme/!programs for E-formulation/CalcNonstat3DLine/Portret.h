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
 *  This file contains the headers for generating a portrait of SLAE matrix for 3D scalar tasks
 *  
 *  
 *  Written by Prof. Yuri G. Soloveichik, Prof. Marina G. Persova, Ph.D. Denis V. Vagin  
 *  Novosibirsk State Technical University,                                              
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                    
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                     
 *  Version 2.0 January 16, 2023                                                                          
*/

#pragma once

struct _list
{
	int number;
	_list *next;
};

// The class contains procedures for constructing a portrait of the SLAE matrix for 3D scalar problem
class Portret
{
public:

	int *ig;
	int *jg;
	_list *s; // 
	bool is_mem_ig_jg_allocated;
	int n; //   
	int n_of_elements;
	int size_jg;
	int (*nver)[14];

	int n_c; //   

	Portret(int (*nver)[14], int n_of_elements, int n_of_nodes, int n_c);
	
	Portret();
	~Portret();

	void Gen_T_Portrait2(); //  3D

	void Add_In_Ordered_List(_list *s, int x); //  -   
	void Clear_List(_list *s);
	int Size_Of_List(_list *s);
};
