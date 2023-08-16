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
 *  This file contains the headers for converting matrix from CSRC to CSR format
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin                         
 *  Novosibirsk State Technical University,                 
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia       
 *  vdv_wk@mail.ru                                          
 *  Version 2.0 January 16, 2023                            
*/

#pragma once
class FormatConverter
{
public:
	void FromRSFToCSR_Real_1_Sym(int nb, int *ig, int *sz_iptr, int *sz_jptr);
	void From2x2ToCSR_Complex_1_Sym(int nb, int *ig, int *idi, int *ijg,int *sz_iptr, int *sz_jptr);
	void From2x2ToCSR_Complex_1_NS(int nb, int *ig, int *idi, int *ijg,int *sz_iptr, int *sz_jptr);
	void From2x2ToCSRComplex_2_Sym(int nb, int *ig, int *jg, int *idi, int *ijg,double *di_block, double *ggl_block,
		MKL_INT64 *iptr, MKL_INT64 *jptr, double *aelem);
	void FromRSFToCSR_Real_2_Sym(int nb, int *ig, int *jg, double *di, double *gg,MKL_INT64 *iptr, MKL_INT64 *jptr, double *aelem);
	void From2x2ToCSRComplex_2_NS(int nb, int *ig, int *jg, int *idi, int *ijg,double *di_block, double *ggl_block, double *ggu_block,
		MKL_INT64 *iptr, MKL_INT64 *jptr, double *aelem);
	FormatConverter();
	~FormatConverter();
};
