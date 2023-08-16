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
 *  This file contains the code for integration using Gauss-3 quadratures
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin                          
 *  Novosibirsk State Technical University,                  
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia        
 *  vdv_wk@mail.ru                                           
 *  Version 2.0 January 16, 2023                             
*/

const double gauss_3_t[] = {-0.774596669241483, 0.0, 0.774596669241483};
const double gauss_3_A[] = {5.0/9.0, 8.0/9.0, 5.0/9.0};
const double gauss_3_t_2d[9][2] = {	gauss_3_t[0], gauss_3_t[0],
									gauss_3_t[1], gauss_3_t[0],
									gauss_3_t[2], gauss_3_t[0],
									gauss_3_t[0], gauss_3_t[1],
									gauss_3_t[1], gauss_3_t[1],
									gauss_3_t[2], gauss_3_t[1],
									gauss_3_t[0], gauss_3_t[2],
									gauss_3_t[1], gauss_3_t[2],
									gauss_3_t[2], gauss_3_t[2]};
const double gauss_3_t_3d[27][3] = {
	gauss_3_t[0], gauss_3_t[0], gauss_3_t[0], 
	gauss_3_t[1], gauss_3_t[0], gauss_3_t[0], 
	gauss_3_t[2], gauss_3_t[0], gauss_3_t[0], 
	gauss_3_t[0], gauss_3_t[1], gauss_3_t[0], 
	gauss_3_t[1], gauss_3_t[1], gauss_3_t[0], 
	gauss_3_t[2], gauss_3_t[1], gauss_3_t[0], 
	gauss_3_t[0], gauss_3_t[2], gauss_3_t[0], 
	gauss_3_t[1], gauss_3_t[2], gauss_3_t[0], 
	gauss_3_t[2], gauss_3_t[2], gauss_3_t[0], 
	gauss_3_t[0], gauss_3_t[0], gauss_3_t[1], 
	gauss_3_t[1], gauss_3_t[0], gauss_3_t[1], 
	gauss_3_t[2], gauss_3_t[0], gauss_3_t[1], 
	gauss_3_t[0], gauss_3_t[1], gauss_3_t[1], 
	gauss_3_t[1], gauss_3_t[1], gauss_3_t[1], 
	gauss_3_t[2], gauss_3_t[1], gauss_3_t[1], 
	gauss_3_t[0], gauss_3_t[2], gauss_3_t[1], 
	gauss_3_t[1], gauss_3_t[2], gauss_3_t[1], 
	gauss_3_t[2], gauss_3_t[2], gauss_3_t[1], 
	gauss_3_t[0], gauss_3_t[0], gauss_3_t[2], 
	gauss_3_t[1], gauss_3_t[0], gauss_3_t[2], 
	gauss_3_t[2], gauss_3_t[0], gauss_3_t[2], 
	gauss_3_t[0], gauss_3_t[1], gauss_3_t[2], 
	gauss_3_t[1], gauss_3_t[1], gauss_3_t[2], 
	gauss_3_t[2], gauss_3_t[1], gauss_3_t[2], 
	gauss_3_t[0], gauss_3_t[2], gauss_3_t[2], 
	gauss_3_t[1], gauss_3_t[2], gauss_3_t[2], 
	gauss_3_t[2], gauss_3_t[2], gauss_3_t[2] 
};
const double gauss_3_A_all[27] =
{
	1.7146776406036e-001,
	2.7434842249657e-001,
	1.7146776406036e-001,
	2.7434842249657e-001,
	4.3895747599451e-001,
	2.7434842249657e-001,
	1.7146776406036e-001,
	2.7434842249657e-001,
	1.7146776406036e-001,
	2.7434842249657e-001,
	4.3895747599451e-001,
	2.7434842249657e-001,
	4.3895747599451e-001,
	7.0233196159122e-001,
	4.3895747599451e-001,
	2.7434842249657e-001,
	4.3895747599451e-001,
	2.7434842249657e-001,
	1.7146776406036e-001,
	2.7434842249657e-001,
	1.7146776406036e-001,
	2.7434842249657e-001,
	4.3895747599451e-001,
	2.7434842249657e-001,
	1.7146776406036e-001,
	2.7434842249657e-001,
	1.7146776406036e-001
};
const double gauss_3_phi[27][8] =
{
 6.9856850115867e-001,	 8.8729833462074e-002,	 8.8729833462074e-002,	 1.1270166537926e-002,	 8.8729833462074e-002,	 1.1270166537926e-002,	 1.1270166537926e-002,	 1.4314988413325e-003,	
 3.9364916731037e-001,	 5.0000000000000e-002,	 5.0000000000000e-002,	 6.3508326896292e-003,	 3.9364916731037e-001,	 5.0000000000000e-002,	 5.0000000000000e-002,	 6.3508326896292e-003,	
 8.8729833462074e-002,	 1.1270166537926e-002,	 1.1270166537926e-002,	 1.4314988413325e-003,	 6.9856850115867e-001,	 8.8729833462074e-002,	 8.8729833462074e-002,	 1.1270166537926e-002,	
 3.9364916731037e-001,	 5.0000000000000e-002,	 3.9364916731037e-001,	 5.0000000000000e-002,	 5.0000000000000e-002,	 6.3508326896292e-003,	 5.0000000000000e-002,	 6.3508326896292e-003,	
 2.2182458365519e-001,	 2.8175416344815e-002,	 2.2182458365519e-001,	 2.8175416344815e-002,	 2.2182458365519e-001,	 2.8175416344815e-002,	 2.2182458365519e-001,	 2.8175416344815e-002,	
 5.0000000000000e-002,	 6.3508326896292e-003,	 5.0000000000000e-002,	 6.3508326896292e-003,	 3.9364916731037e-001,	 5.0000000000000e-002,	 3.9364916731037e-001,	 5.0000000000000e-002,	
 8.8729833462074e-002,	 1.1270166537926e-002,	 6.9856850115867e-001,	 8.8729833462074e-002,	 1.1270166537926e-002,	 1.4314988413325e-003,	 8.8729833462074e-002,	 1.1270166537926e-002,	
 5.0000000000000e-002,	 6.3508326896292e-003,	 3.9364916731037e-001,	 5.0000000000000e-002,	 5.0000000000000e-002,	 6.3508326896292e-003,	 3.9364916731037e-001,	 5.0000000000000e-002,	
 1.1270166537926e-002,	 1.4314988413325e-003,	 8.8729833462074e-002,	 1.1270166537926e-002,	 8.8729833462074e-002,	 1.1270166537926e-002,	 6.9856850115867e-001,	 8.8729833462074e-002,	
 3.9364916731037e-001,	 3.9364916731037e-001,	 5.0000000000000e-002,	 5.0000000000000e-002,	 5.0000000000000e-002,	 5.0000000000000e-002,	 6.3508326896292e-003,	 6.3508326896292e-003,	
 2.2182458365519e-001,	 2.2182458365519e-001,	 2.8175416344815e-002,	 2.8175416344815e-002,	 2.2182458365519e-001,	 2.2182458365519e-001,	 2.8175416344815e-002,	 2.8175416344815e-002,	
 5.0000000000000e-002,	 5.0000000000000e-002,	 6.3508326896292e-003,	 6.3508326896292e-003,	 3.9364916731037e-001,	 3.9364916731037e-001,	 5.0000000000000e-002,	 5.0000000000000e-002,	
 2.2182458365519e-001,	 2.2182458365519e-001,	 2.2182458365519e-001,	 2.2182458365519e-001,	 2.8175416344815e-002,	 2.8175416344815e-002,	 2.8175416344815e-002,	 2.8175416344815e-002,	
 1.2500000000000e-001,	 1.2500000000000e-001,	 1.2500000000000e-001,	 1.2500000000000e-001,	 1.2500000000000e-001,	 1.2500000000000e-001,	 1.2500000000000e-001,	 1.2500000000000e-001,	
 2.8175416344815e-002,	 2.8175416344815e-002,	 2.8175416344815e-002,	 2.8175416344815e-002,	 2.2182458365519e-001,	 2.2182458365519e-001,	 2.2182458365519e-001,	 2.2182458365519e-001,	
 5.0000000000000e-002,	 5.0000000000000e-002,	 3.9364916731037e-001,	 3.9364916731037e-001,	 6.3508326896292e-003,	 6.3508326896292e-003,	 5.0000000000000e-002,	 5.0000000000000e-002,	
 2.8175416344815e-002,	 2.8175416344815e-002,	 2.2182458365519e-001,	 2.2182458365519e-001,	 2.8175416344815e-002,	 2.8175416344815e-002,	 2.2182458365519e-001,	 2.2182458365519e-001,	
 6.3508326896292e-003,	 6.3508326896292e-003,	 5.0000000000000e-002,	 5.0000000000000e-002,	 5.0000000000000e-002,	 5.0000000000000e-002,	 3.9364916731037e-001,	 3.9364916731037e-001,	
 8.8729833462074e-002,	 6.9856850115867e-001,	 1.1270166537926e-002,	 8.8729833462074e-002,	 1.1270166537926e-002,	 8.8729833462074e-002,	 1.4314988413325e-003,	 1.1270166537926e-002,	
 5.0000000000000e-002,	 3.9364916731037e-001,	 6.3508326896292e-003,	 5.0000000000000e-002,	 5.0000000000000e-002,	 3.9364916731037e-001,	 6.3508326896292e-003,	 5.0000000000000e-002,	
 1.1270166537926e-002,	 8.8729833462074e-002,	 1.4314988413325e-003,	 1.1270166537926e-002,	 8.8729833462074e-002,	 6.9856850115867e-001,	 1.1270166537926e-002,	 8.8729833462074e-002,	
 5.0000000000000e-002,	 3.9364916731037e-001,	 5.0000000000000e-002,	 3.9364916731037e-001,	 6.3508326896292e-003,	 5.0000000000000e-002,	 6.3508326896292e-003,	 5.0000000000000e-002,	
 2.8175416344815e-002,	 2.2182458365519e-001,	 2.8175416344815e-002,	 2.2182458365519e-001,	 2.8175416344815e-002,	 2.2182458365519e-001,	 2.8175416344815e-002,	 2.2182458365519e-001,	
 6.3508326896292e-003,	 5.0000000000000e-002,	 6.3508326896292e-003,	 5.0000000000000e-002,	 5.0000000000000e-002,	 3.9364916731037e-001,	 5.0000000000000e-002,	 3.9364916731037e-001,	
 1.1270166537926e-002,	 8.8729833462074e-002,	 8.8729833462074e-002,	 6.9856850115867e-001,	 1.4314988413325e-003,	 1.1270166537926e-002,	 1.1270166537926e-002,	 8.8729833462074e-002,	
 6.3508326896292e-003,	 5.0000000000000e-002,	 5.0000000000000e-002,	 3.9364916731037e-001,	 6.3508326896292e-003,	 5.0000000000000e-002,	 5.0000000000000e-002,	 3.9364916731037e-001,	
 1.4314988413325e-003,	 1.1270166537926e-002,	 1.1270166537926e-002,	 8.8729833462074e-002,	 1.1270166537926e-002,	 8.8729833462074e-002,	 8.8729833462074e-002,	 6.9856850115867e-001	
};
const double gauss_3_d_phi[27][8][3] =
{
	-3.9364916731037e-001,	-3.9364916731037e-001,	-3.9364916731037e-001,
	3.9364916731037e-001,	-5.0000000000000e-002,	-5.0000000000000e-002,
	-5.0000000000000e-002,	 3.9364916731037e-001,	-5.0000000000000e-002,
	5.0000000000000e-002,	 5.0000000000000e-002,	-6.3508326896292e-003,
	-5.0000000000000e-002,	-5.0000000000000e-002,	 3.9364916731037e-001,
	5.0000000000000e-002,	-6.3508326896292e-003,	 5.0000000000000e-002,
	-6.3508326896292e-003,	 5.0000000000000e-002,	 5.0000000000000e-002,
	6.3508326896292e-003,	 6.3508326896292e-003,	 6.3508326896292e-003,
	-2.2182458365519e-001,	-2.2182458365519e-001,	-3.9364916731037e-001,
	2.2182458365519e-001,	-2.8175416344815e-002,	-5.0000000000000e-002,
	-2.8175416344815e-002,	 2.2182458365519e-001,	-5.0000000000000e-002,
	2.8175416344815e-002,	 2.8175416344815e-002,	-6.3508326896292e-003,
	-2.2182458365519e-001,	-2.2182458365519e-001,	 3.9364916731037e-001,
	2.2182458365519e-001,	-2.8175416344815e-002,	 5.0000000000000e-002,
	-2.8175416344815e-002,	 2.2182458365519e-001,	 5.0000000000000e-002,
	2.8175416344815e-002,	 2.8175416344815e-002,	 6.3508326896292e-003,
	-5.0000000000000e-002,	-5.0000000000000e-002,	-3.9364916731037e-001,
	5.0000000000000e-002,	-6.3508326896292e-003,	-5.0000000000000e-002,
	-6.3508326896292e-003,	 5.0000000000000e-002,	-5.0000000000000e-002,
	6.3508326896292e-003,	 6.3508326896292e-003,	-6.3508326896292e-003,
	-3.9364916731037e-001,	-3.9364916731037e-001,	 3.9364916731037e-001,
	3.9364916731037e-001,	-5.0000000000000e-002,	 5.0000000000000e-002,
	-5.0000000000000e-002,	 3.9364916731037e-001,	 5.0000000000000e-002,
	5.0000000000000e-002,	 5.0000000000000e-002,	 6.3508326896292e-003,
	-2.2182458365519e-001,	-3.9364916731037e-001,	-2.2182458365519e-001,
	2.2182458365519e-001,	-5.0000000000000e-002,	-2.8175416344815e-002,
	-2.2182458365519e-001,	 3.9364916731037e-001,	-2.2182458365519e-001,
	2.2182458365519e-001,	 5.0000000000000e-002,	-2.8175416344815e-002,
	-2.8175416344815e-002,	-5.0000000000000e-002,	 2.2182458365519e-001,
	2.8175416344815e-002,	-6.3508326896292e-003,	 2.8175416344815e-002,
	-2.8175416344815e-002,	 5.0000000000000e-002,	 2.2182458365519e-001,
	2.8175416344815e-002,	 6.3508326896292e-003,	 2.8175416344815e-002,
	-1.2500000000000e-001,	-2.2182458365519e-001,	-2.2182458365519e-001,
	1.2500000000000e-001,	-2.8175416344815e-002,	-2.8175416344815e-002,
	-1.2500000000000e-001,	 2.2182458365519e-001,	-2.2182458365519e-001,
	1.2500000000000e-001,	 2.8175416344815e-002,	-2.8175416344815e-002,
	-1.2500000000000e-001,	-2.2182458365519e-001,	 2.2182458365519e-001,
	1.2500000000000e-001,	-2.8175416344815e-002,	 2.8175416344815e-002,
	-1.2500000000000e-001,	 2.2182458365519e-001,	 2.2182458365519e-001,
	1.2500000000000e-001,	 2.8175416344815e-002,	 2.8175416344815e-002,
	-2.8175416344815e-002,	-5.0000000000000e-002,	-2.2182458365519e-001,
	2.8175416344815e-002,	-6.3508326896292e-003,	-2.8175416344815e-002,
	-2.8175416344815e-002,	 5.0000000000000e-002,	-2.2182458365519e-001,
	2.8175416344815e-002,	 6.3508326896292e-003,	-2.8175416344815e-002,
	-2.2182458365519e-001,	-3.9364916731037e-001,	 2.2182458365519e-001,
	2.2182458365519e-001,	-5.0000000000000e-002,	 2.8175416344815e-002,
	-2.2182458365519e-001,	 3.9364916731037e-001,	 2.2182458365519e-001,
	2.2182458365519e-001,	 5.0000000000000e-002,	 2.8175416344815e-002,
	-5.0000000000000e-002,	-3.9364916731037e-001,	-5.0000000000000e-002,
	5.0000000000000e-002,	-5.0000000000000e-002,	-6.3508326896292e-003,
	-3.9364916731037e-001,	 3.9364916731037e-001,	-3.9364916731037e-001,
	3.9364916731037e-001,	 5.0000000000000e-002,	-5.0000000000000e-002,
	-6.3508326896292e-003,	-5.0000000000000e-002,	 5.0000000000000e-002,
	6.3508326896292e-003,	-6.3508326896292e-003,	 6.3508326896292e-003,
	-5.0000000000000e-002,	 5.0000000000000e-002,	 3.9364916731037e-001,
	5.0000000000000e-002,	 6.3508326896292e-003,	 5.0000000000000e-002,
	-2.8175416344815e-002,	-2.2182458365519e-001,	-5.0000000000000e-002,
	2.8175416344815e-002,	-2.8175416344815e-002,	-6.3508326896292e-003,
	-2.2182458365519e-001,	 2.2182458365519e-001,	-3.9364916731037e-001,
	2.2182458365519e-001,	 2.8175416344815e-002,	-5.0000000000000e-002,
	-2.8175416344815e-002,	-2.2182458365519e-001,	 5.0000000000000e-002,
	2.8175416344815e-002,	-2.8175416344815e-002,	 6.3508326896292e-003,
	-2.2182458365519e-001,	 2.2182458365519e-001,	 3.9364916731037e-001,
	2.2182458365519e-001,	 2.8175416344815e-002,	 5.0000000000000e-002,
	-6.3508326896292e-003,	-5.0000000000000e-002,	-5.0000000000000e-002,
	6.3508326896292e-003,	-6.3508326896292e-003,	-6.3508326896292e-003,
	-5.0000000000000e-002,	 5.0000000000000e-002,	-3.9364916731037e-001,
	5.0000000000000e-002,	 6.3508326896292e-003,	-5.0000000000000e-002,
	-5.0000000000000e-002,	-3.9364916731037e-001,	 5.0000000000000e-002,
	5.0000000000000e-002,	-5.0000000000000e-002,	 6.3508326896292e-003,
	-3.9364916731037e-001,	 3.9364916731037e-001,	 3.9364916731037e-001,
	3.9364916731037e-001,	 5.0000000000000e-002,	 5.0000000000000e-002,
	-3.9364916731037e-001,	-2.2182458365519e-001,	-2.2182458365519e-001,
	3.9364916731037e-001,	-2.2182458365519e-001,	-2.2182458365519e-001,
	-5.0000000000000e-002,	 2.2182458365519e-001,	-2.8175416344815e-002,
	5.0000000000000e-002,	 2.2182458365519e-001,	-2.8175416344815e-002,
	-5.0000000000000e-002,	-2.8175416344815e-002,	 2.2182458365519e-001,
	5.0000000000000e-002,	-2.8175416344815e-002,	 2.2182458365519e-001,
	-6.3508326896292e-003,	 2.8175416344815e-002,	 2.8175416344815e-002,
	6.3508326896292e-003,	 2.8175416344815e-002,	 2.8175416344815e-002,
	-2.2182458365519e-001,	-1.2500000000000e-001,	-2.2182458365519e-001,
	2.2182458365519e-001,	-1.2500000000000e-001,	-2.2182458365519e-001,
	-2.8175416344815e-002,	 1.2500000000000e-001,	-2.8175416344815e-002,
	2.8175416344815e-002,	 1.2500000000000e-001,	-2.8175416344815e-002,
	-2.2182458365519e-001,	-1.2500000000000e-001,	 2.2182458365519e-001,
	2.2182458365519e-001,	-1.2500000000000e-001,	 2.2182458365519e-001,
	-2.8175416344815e-002,	 1.2500000000000e-001,	 2.8175416344815e-002,
	2.8175416344815e-002,	 1.2500000000000e-001,	 2.8175416344815e-002,
	-5.0000000000000e-002,	-2.8175416344815e-002,	-2.2182458365519e-001,
	5.0000000000000e-002,	-2.8175416344815e-002,	-2.2182458365519e-001,
	-6.3508326896292e-003,	 2.8175416344815e-002,	-2.8175416344815e-002,
	6.3508326896292e-003,	 2.8175416344815e-002,	-2.8175416344815e-002,
	-3.9364916731037e-001,	-2.2182458365519e-001,	 2.2182458365519e-001,
	3.9364916731037e-001,	-2.2182458365519e-001,	 2.2182458365519e-001,
	-5.0000000000000e-002,	 2.2182458365519e-001,	 2.8175416344815e-002,
	5.0000000000000e-002,	 2.2182458365519e-001,	 2.8175416344815e-002,
	-2.2182458365519e-001,	-2.2182458365519e-001,	-1.2500000000000e-001,
	2.2182458365519e-001,	-2.2182458365519e-001,	-1.2500000000000e-001,
	-2.2182458365519e-001,	 2.2182458365519e-001,	-1.2500000000000e-001,
	2.2182458365519e-001,	 2.2182458365519e-001,	-1.2500000000000e-001,
	-2.8175416344815e-002,	-2.8175416344815e-002,	 1.2500000000000e-001,
	2.8175416344815e-002,	-2.8175416344815e-002,	 1.2500000000000e-001,
	-2.8175416344815e-002,	 2.8175416344815e-002,	 1.2500000000000e-001,
	2.8175416344815e-002,	 2.8175416344815e-002,	 1.2500000000000e-001,
	-1.2500000000000e-001,	-1.2500000000000e-001,	-1.2500000000000e-001,
	1.2500000000000e-001,	-1.2500000000000e-001,	-1.2500000000000e-001,
	-1.2500000000000e-001,	 1.2500000000000e-001,	-1.2500000000000e-001,
	1.2500000000000e-001,	 1.2500000000000e-001,	-1.2500000000000e-001,
	-1.2500000000000e-001,	-1.2500000000000e-001,	 1.2500000000000e-001,
	1.2500000000000e-001,	-1.2500000000000e-001,	 1.2500000000000e-001,
	-1.2500000000000e-001,	 1.2500000000000e-001,	 1.2500000000000e-001,
	1.2500000000000e-001,	 1.2500000000000e-001,	 1.2500000000000e-001,
	-2.8175416344815e-002,	-2.8175416344815e-002,	-1.2500000000000e-001,
	2.8175416344815e-002,	-2.8175416344815e-002,	-1.2500000000000e-001,
	-2.8175416344815e-002,	 2.8175416344815e-002,	-1.2500000000000e-001,
	2.8175416344815e-002,	 2.8175416344815e-002,	-1.2500000000000e-001,
	-2.2182458365519e-001,	-2.2182458365519e-001,	 1.2500000000000e-001,
	2.2182458365519e-001,	-2.2182458365519e-001,	 1.2500000000000e-001,
	-2.2182458365519e-001,	 2.2182458365519e-001,	 1.2500000000000e-001,
	2.2182458365519e-001,	 2.2182458365519e-001,	 1.2500000000000e-001,
	-5.0000000000000e-002,	-2.2182458365519e-001,	-2.8175416344815e-002,
	5.0000000000000e-002,	-2.2182458365519e-001,	-2.8175416344815e-002,
	-3.9364916731037e-001,	 2.2182458365519e-001,	-2.2182458365519e-001,
	3.9364916731037e-001,	 2.2182458365519e-001,	-2.2182458365519e-001,
	-6.3508326896292e-003,	-2.8175416344815e-002,	 2.8175416344815e-002,
	6.3508326896292e-003,	-2.8175416344815e-002,	 2.8175416344815e-002,
	-5.0000000000000e-002,	 2.8175416344815e-002,	 2.2182458365519e-001,
	5.0000000000000e-002,	 2.8175416344815e-002,	 2.2182458365519e-001,
	-2.8175416344815e-002,	-1.2500000000000e-001,	-2.8175416344815e-002,
	2.8175416344815e-002,	-1.2500000000000e-001,	-2.8175416344815e-002,
	-2.2182458365519e-001,	 1.2500000000000e-001,	-2.2182458365519e-001,
	2.2182458365519e-001,	 1.2500000000000e-001,	-2.2182458365519e-001,
	-2.8175416344815e-002,	-1.2500000000000e-001,	 2.8175416344815e-002,
	2.8175416344815e-002,	-1.2500000000000e-001,	 2.8175416344815e-002,
	-2.2182458365519e-001,	 1.2500000000000e-001,	 2.2182458365519e-001,
	2.2182458365519e-001,	 1.2500000000000e-001,	 2.2182458365519e-001,
	-6.3508326896292e-003,	-2.8175416344815e-002,	-2.8175416344815e-002,
	6.3508326896292e-003,	-2.8175416344815e-002,	-2.8175416344815e-002,
	-5.0000000000000e-002,	 2.8175416344815e-002,	-2.2182458365519e-001,
	5.0000000000000e-002,	 2.8175416344815e-002,	-2.2182458365519e-001,
	-5.0000000000000e-002,	-2.2182458365519e-001,	 2.8175416344815e-002,
	5.0000000000000e-002,	-2.2182458365519e-001,	 2.8175416344815e-002,
	-3.9364916731037e-001,	 2.2182458365519e-001,	 2.2182458365519e-001,
	3.9364916731037e-001,	 2.2182458365519e-001,	 2.2182458365519e-001,
	-3.9364916731037e-001,	-5.0000000000000e-002,	-5.0000000000000e-002,
	3.9364916731037e-001,	-3.9364916731037e-001,	-3.9364916731037e-001,
	-5.0000000000000e-002,	 5.0000000000000e-002,	-6.3508326896292e-003,
	5.0000000000000e-002,	 3.9364916731037e-001,	-5.0000000000000e-002,
	-5.0000000000000e-002,	-6.3508326896292e-003,	 5.0000000000000e-002,
	5.0000000000000e-002,	-5.0000000000000e-002,	 3.9364916731037e-001,
	-6.3508326896292e-003,	 6.3508326896292e-003,	 6.3508326896292e-003,
	6.3508326896292e-003,	 5.0000000000000e-002,	 5.0000000000000e-002,
	-2.2182458365519e-001,	-2.8175416344815e-002,	-5.0000000000000e-002,
	2.2182458365519e-001,	-2.2182458365519e-001,	-3.9364916731037e-001,
	-2.8175416344815e-002,	 2.8175416344815e-002,	-6.3508326896292e-003,
	2.8175416344815e-002,	 2.2182458365519e-001,	-5.0000000000000e-002,
	-2.2182458365519e-001,	-2.8175416344815e-002,	 5.0000000000000e-002,
	2.2182458365519e-001,	-2.2182458365519e-001,	 3.9364916731037e-001,
	-2.8175416344815e-002,	 2.8175416344815e-002,	 6.3508326896292e-003,
	2.8175416344815e-002,	 2.2182458365519e-001,	 5.0000000000000e-002,
	-5.0000000000000e-002,	-6.3508326896292e-003,	-5.0000000000000e-002,
	5.0000000000000e-002,	-5.0000000000000e-002,	-3.9364916731037e-001,
	-6.3508326896292e-003,	 6.3508326896292e-003,	-6.3508326896292e-003,
	6.3508326896292e-003,	 5.0000000000000e-002,	-5.0000000000000e-002,
	-3.9364916731037e-001,	-5.0000000000000e-002,	 5.0000000000000e-002,
	3.9364916731037e-001,	-3.9364916731037e-001,	 3.9364916731037e-001,
	-5.0000000000000e-002,	 5.0000000000000e-002,	 6.3508326896292e-003,
	5.0000000000000e-002,	 3.9364916731037e-001,	 5.0000000000000e-002,
	-2.2182458365519e-001,	-5.0000000000000e-002,	-2.8175416344815e-002,
	2.2182458365519e-001,	-3.9364916731037e-001,	-2.2182458365519e-001,
	-2.2182458365519e-001,	 5.0000000000000e-002,	-2.8175416344815e-002,
	2.2182458365519e-001,	 3.9364916731037e-001,	-2.2182458365519e-001,
	-2.8175416344815e-002,	-6.3508326896292e-003,	 2.8175416344815e-002,
	2.8175416344815e-002,	-5.0000000000000e-002,	 2.2182458365519e-001,
	-2.8175416344815e-002,	 6.3508326896292e-003,	 2.8175416344815e-002,
	2.8175416344815e-002,	 5.0000000000000e-002,	 2.2182458365519e-001,
	-1.2500000000000e-001,	-2.8175416344815e-002,	-2.8175416344815e-002,
	1.2500000000000e-001,	-2.2182458365519e-001,	-2.2182458365519e-001,
	-1.2500000000000e-001,	 2.8175416344815e-002,	-2.8175416344815e-002,
	1.2500000000000e-001,	 2.2182458365519e-001,	-2.2182458365519e-001,
	-1.2500000000000e-001,	-2.8175416344815e-002,	 2.8175416344815e-002,
	1.2500000000000e-001,	-2.2182458365519e-001,	 2.2182458365519e-001,
	-1.2500000000000e-001,	 2.8175416344815e-002,	 2.8175416344815e-002,
	1.2500000000000e-001,	 2.2182458365519e-001,	 2.2182458365519e-001,
	-2.8175416344815e-002,	-6.3508326896292e-003,	-2.8175416344815e-002,
	2.8175416344815e-002,	-5.0000000000000e-002,	-2.2182458365519e-001,
	-2.8175416344815e-002,	 6.3508326896292e-003,	-2.8175416344815e-002,
	2.8175416344815e-002,	 5.0000000000000e-002,	-2.2182458365519e-001,
	-2.2182458365519e-001,	-5.0000000000000e-002,	 2.8175416344815e-002,
	2.2182458365519e-001,	-3.9364916731037e-001,	 2.2182458365519e-001,
	-2.2182458365519e-001,	 5.0000000000000e-002,	 2.8175416344815e-002,
	2.2182458365519e-001,	 3.9364916731037e-001,	 2.2182458365519e-001,
	-5.0000000000000e-002,	-5.0000000000000e-002,	-6.3508326896292e-003,
	5.0000000000000e-002,	-3.9364916731037e-001,	-5.0000000000000e-002,
	-3.9364916731037e-001,	 5.0000000000000e-002,	-5.0000000000000e-002,
	3.9364916731037e-001,	 3.9364916731037e-001,	-3.9364916731037e-001,
	-6.3508326896292e-003,	-6.3508326896292e-003,	 6.3508326896292e-003,
	6.3508326896292e-003,	-5.0000000000000e-002,	 5.0000000000000e-002,
	-5.0000000000000e-002,	 6.3508326896292e-003,	 5.0000000000000e-002,
	5.0000000000000e-002,	 5.0000000000000e-002,	 3.9364916731037e-001,
	-2.8175416344815e-002,	-2.8175416344815e-002,	-6.3508326896292e-003,
	2.8175416344815e-002,	-2.2182458365519e-001,	-5.0000000000000e-002,
	-2.2182458365519e-001,	 2.8175416344815e-002,	-5.0000000000000e-002,
	2.2182458365519e-001,	 2.2182458365519e-001,	-3.9364916731037e-001,
	-2.8175416344815e-002,	-2.8175416344815e-002,	 6.3508326896292e-003,
	2.8175416344815e-002,	-2.2182458365519e-001,	 5.0000000000000e-002,
	-2.2182458365519e-001,	 2.8175416344815e-002,	 5.0000000000000e-002,
	2.2182458365519e-001,	 2.2182458365519e-001,	 3.9364916731037e-001,
	-6.3508326896292e-003,	-6.3508326896292e-003,	-6.3508326896292e-003,
	6.3508326896292e-003,	-5.0000000000000e-002,	-5.0000000000000e-002,
	-5.0000000000000e-002,	 6.3508326896292e-003,	-5.0000000000000e-002,
	5.0000000000000e-002,	 5.0000000000000e-002,	-3.9364916731037e-001,
	-5.0000000000000e-002,	-5.0000000000000e-002,	 6.3508326896292e-003,
	5.0000000000000e-002,	-3.9364916731037e-001,	 5.0000000000000e-002,
	-3.9364916731037e-001,	 5.0000000000000e-002,	 5.0000000000000e-002,
	3.9364916731037e-001,	 3.9364916731037e-001,	 3.9364916731037e-001
};
const double gauss_3_d_phi_face[9][8][3] =
{
	0.0000000000000e+000,	 0.0000000000000e+000,	-3.9364916731037e-001,
	0.0000000000000e+000,	 0.0000000000000e+000,	-5.0000000000000e-002,
	0.0000000000000e+000,	 0.0000000000000e+000,	-5.0000000000000e-002,
	0.0000000000000e+000,	 0.0000000000000e+000,	-6.3508326896292e-003,
	-4.4364916731037e-001,	-4.4364916731037e-001,	 3.9364916731037e-001,
	4.4364916731037e-001,	-5.6350832689629e-002,	 5.0000000000000e-002,
	-5.6350832689629e-002,	 4.4364916731037e-001,	 5.0000000000000e-002,
	5.6350832689629e-002,	 5.6350832689629e-002,	 6.3508326896292e-003,
	0.0000000000000e+000,	 0.0000000000000e+000,	-2.2182458365519e-001,
	0.0000000000000e+000,	 0.0000000000000e+000,	-2.8175416344815e-002,
	0.0000000000000e+000,	 0.0000000000000e+000,	-2.2182458365519e-001,
	0.0000000000000e+000,	 0.0000000000000e+000,	-2.8175416344815e-002,
	-2.5000000000000e-001,	-4.4364916731037e-001,	 2.2182458365519e-001,
	2.5000000000000e-001,	-5.6350832689629e-002,	 2.8175416344815e-002,
	-2.5000000000000e-001,	 4.4364916731037e-001,	 2.2182458365519e-001,
	2.5000000000000e-001,	 5.6350832689629e-002,	 2.8175416344815e-002,
	0.0000000000000e+000,	 0.0000000000000e+000,	-5.0000000000000e-002,
	0.0000000000000e+000,	 0.0000000000000e+000,	-6.3508326896292e-003,
	0.0000000000000e+000,	 0.0000000000000e+000,	-3.9364916731037e-001,
	0.0000000000000e+000,	 0.0000000000000e+000,	-5.0000000000000e-002,
	-5.6350832689629e-002,	-4.4364916731037e-001,	 5.0000000000000e-002,
	5.6350832689629e-002,	-5.6350832689629e-002,	 6.3508326896292e-003,
	-4.4364916731037e-001,	 4.4364916731037e-001,	 3.9364916731037e-001,
	4.4364916731037e-001,	 5.6350832689629e-002,	 5.0000000000000e-002,
	0.0000000000000e+000,	 0.0000000000000e+000,	-2.2182458365519e-001,
	0.0000000000000e+000,	 0.0000000000000e+000,	-2.2182458365519e-001,
	0.0000000000000e+000,	 0.0000000000000e+000,	-2.8175416344815e-002,
	0.0000000000000e+000,	 0.0000000000000e+000,	-2.8175416344815e-002,
	-4.4364916731037e-001,	-2.5000000000000e-001,	 2.2182458365519e-001,
	4.4364916731037e-001,	-2.5000000000000e-001,	 2.2182458365519e-001,
	-5.6350832689629e-002,	 2.5000000000000e-001,	 2.8175416344815e-002,
	5.6350832689629e-002,	 2.5000000000000e-001,	 2.8175416344815e-002,
	0.0000000000000e+000,	 0.0000000000000e+000,	-1.2500000000000e-001,
	0.0000000000000e+000,	 0.0000000000000e+000,	-1.2500000000000e-001,
	0.0000000000000e+000,	 0.0000000000000e+000,	-1.2500000000000e-001,
	0.0000000000000e+000,	 0.0000000000000e+000,	-1.2500000000000e-001,
	-2.5000000000000e-001,	-2.5000000000000e-001,	 1.2500000000000e-001,
	2.5000000000000e-001,	-2.5000000000000e-001,	 1.2500000000000e-001,
	-2.5000000000000e-001,	 2.5000000000000e-001,	 1.2500000000000e-001,
	2.5000000000000e-001,	 2.5000000000000e-001,	 1.2500000000000e-001,
	0.0000000000000e+000,	 0.0000000000000e+000,	-2.8175416344815e-002,
	0.0000000000000e+000,	 0.0000000000000e+000,	-2.8175416344815e-002,
	0.0000000000000e+000,	 0.0000000000000e+000,	-2.2182458365519e-001,
	0.0000000000000e+000,	 0.0000000000000e+000,	-2.2182458365519e-001,
	-5.6350832689629e-002,	-2.5000000000000e-001,	 2.8175416344815e-002,
	5.6350832689629e-002,	-2.5000000000000e-001,	 2.8175416344815e-002,
	-4.4364916731037e-001,	 2.5000000000000e-001,	 2.2182458365519e-001,
	4.4364916731037e-001,	 2.5000000000000e-001,	 2.2182458365519e-001,
	0.0000000000000e+000,	 0.0000000000000e+000,	-5.0000000000000e-002,
	0.0000000000000e+000,	 0.0000000000000e+000,	-3.9364916731037e-001,
	0.0000000000000e+000,	 0.0000000000000e+000,	-6.3508326896292e-003,
	0.0000000000000e+000,	 0.0000000000000e+000,	-5.0000000000000e-002,
	-4.4364916731037e-001,	-5.6350832689629e-002,	 5.0000000000000e-002,
	4.4364916731037e-001,	-4.4364916731037e-001,	 3.9364916731037e-001,
	-5.6350832689629e-002,	 5.6350832689629e-002,	 6.3508326896292e-003,
	5.6350832689629e-002,	 4.4364916731037e-001,	 5.0000000000000e-002,
	0.0000000000000e+000,	 0.0000000000000e+000,	-2.8175416344815e-002,
	0.0000000000000e+000,	 0.0000000000000e+000,	-2.2182458365519e-001,
	0.0000000000000e+000,	 0.0000000000000e+000,	-2.8175416344815e-002,
	0.0000000000000e+000,	 0.0000000000000e+000,	-2.2182458365519e-001,
	-2.5000000000000e-001,	-5.6350832689629e-002,	 2.8175416344815e-002,
	2.5000000000000e-001,	-4.4364916731037e-001,	 2.2182458365519e-001,
	-2.5000000000000e-001,	 5.6350832689629e-002,	 2.8175416344815e-002,
	2.5000000000000e-001,	 4.4364916731037e-001,	 2.2182458365519e-001,
	0.0000000000000e+000,	 0.0000000000000e+000,	-6.3508326896292e-003,
	0.0000000000000e+000,	 0.0000000000000e+000,	-5.0000000000000e-002,
	0.0000000000000e+000,	 0.0000000000000e+000,	-5.0000000000000e-002,
	0.0000000000000e+000,	 0.0000000000000e+000,	-3.9364916731037e-001,
	-5.6350832689629e-002,	-5.6350832689629e-002,	 6.3508326896292e-003,
	5.6350832689629e-002,	-4.4364916731037e-001,	 5.0000000000000e-002,
	-4.4364916731037e-001,	 5.6350832689629e-002,	 5.0000000000000e-002,
	4.4364916731037e-001,	 4.4364916731037e-001,	 3.9364916731037e-001
};

const double gauss_3_A_all_2[9] =
{
	0.3086419753086419,
	0.4938271604938271,
	0.3086419753086419,
	0.4938271604938271,
	0.7901234567901234,
	0.4938271604938271,
	0.3086419753086419,
	0.4938271604938271,
	0.3086419753086419
};

const double gauss_3_phi_2[9][4] =
{
	0.78729833462074,	0.1,				0.1,				0.012701665379258,
	0.44364916731037,	0.44364916731037,	0.056350832689629,	0.056350832689629,
	0.1,				0.78729833462074,	0.012701665379258,	0.1,
	0.44364916731037,	0.056350832689629,	0.44364916731037,	0.056350832689629,
	0.25,				0.25,				0.25,				0.25,
	0.056350832689629,	0.44364916731037,	0.056350832689629,	0.44364916731037,
	0.1,				0.012701665379258,	0.78729833462074,	0.1,
	0.056350832689629,	0.056350832689629,	0.44364916731037,	0.44364916731037,
	0.012701665379258,	0.1,				0.1,				0.78729833462074
};

const double gauss_3_d_phi_2[9][4][3] =
{
	-0.44364916731037,	-0.44364916731037,	-0.39364916731037,
	 0.44364916731037,	-0.056350832689629,	-0.05,
	-0.056350832689629,	 0.44364916731037,	-0.05,
	 0.056350832689629,	 0.056350832689629,	-0.0063508326896291,

	-0.44364916731037,	-0.25,				-0.221824583655185,
	 0.44364916731037,	-0.25,				-0.221824583655185,
	-0.056350832689629,	 0.25,				-0.028175416344815,
	 0.056350832689629,	 0.25,				-0.028175416344815,

	-0.44364916731037,	-0.056350832689629,	-0.05,
	 0.44364916731037,	-0.44364916731037,	-0.39364916731037,
	-0.056350832689629,	 0.056350832689629,	-0.0063508326896291,
	 0.056350832689629,	 0.44364916731037,	-0.05,

	-0.25,				-0.44364916731037,	-0.221824583655185,
	 0.25,				-0.056350832689629,	-0.028175416344815,
	-0.25,				 0.44364916731037,	-0.221824583655185,
	 0.25,				 0.056350832689629,	-0.028175416344815,

 	-0.25,				-0.25,				-0.125,
	 0.25,				-0.25,				-0.125,
	-0.25,				 0.25,				-0.125,
	 0.25,				 0.25,				-0.125,

 	-0.25,				-0.056350832689629,	-0.028175416344815,
	 0.25,				-0.44364916731037,	-0.221824583655185,
	-0.25,				 0.056350832689629,	-0.028175416344815,
	 0.25,				 0.44364916731037,	-0.221824583655185,

	-0.056350832689629,	-0.44364916731037,	-0.05,
	 0.056350832689629,	-0.056350832689629,	-0.0063508326896291,
	-0.44364916731037,	 0.44364916731037,	-0.39364916731037,
	 0.44364916731037,	 0.056350832689629,	-0.05,

	-0.056350832689629,	-0.25,				-0.028175416344815,
	 0.056350832689629,	-0.25,				-0.028175416344815,
	-0.44364916731037,	 0.25,				-0.221824583655185,
	 0.44364916731037,	 0.25,				-0.221824583655185,

	-0.056350832689629,	-0.056350832689629,	-0.0063508326896291,
	 0.056350832689629,	-0.44364916731037,	-0.05,
	-0.44364916731037,	 0.056350832689629,	-0.05,
	 0.44364916731037,	 0.44364916731037,	-0.39364916731037
};
