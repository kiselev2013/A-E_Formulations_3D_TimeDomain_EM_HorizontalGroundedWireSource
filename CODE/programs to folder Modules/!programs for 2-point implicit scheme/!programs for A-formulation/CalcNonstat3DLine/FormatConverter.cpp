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
 *  This file contains the code for converting matrix from CSRC to CSR format
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin 
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 January 16, 2023
*/

#include "stdafx.h"
#include "FormatConverter.h"
FormatConverter::FormatConverter()
{
}
FormatConverter::~FormatConverter()
{
}
void FormatConverter::From2x2ToCSR_Complex_1_Sym(int nb, int *ig, int *idi, int *ijg,int *sz_iptr, int *sz_jptr)
{
	*sz_iptr = nb+1;
	*sz_jptr = ig[nb] + nb;
}
void FormatConverter::FromRSFToCSR_Real_1_Sym(int nb, int *ig, int *sz_iptr, int *sz_jptr)
{
	*sz_iptr = nb+1;
	*sz_jptr = ig[nb] + nb;
}
void FormatConverter::From2x2ToCSR_Complex_1_NS(int nb, int *ig, int *idi, int *ijg,int *sz_iptr, int *sz_jptr)
{
	*sz_iptr = nb+1;
	*sz_jptr = ig[nb]*2 + nb;
}
void FormatConverter::FromRSFToCSR_Real_2_Sym(int nb, int *ig, int *jg, double *di, double *gg,MKL_INT64 *iptr, MKL_INT64 *jptr, double *aelem)
{
	int sz;
	int i, j, k, m;
	vector<MKL_INT64> col; 

	col.resize(nb, 0);

	for (i=0; i<nb; i++)
	{
		col[i] += 1; // 

		for (j=ig[i]; j<=ig[i+1]-1; j++)
		{
			k = jg[j];
			col[k]++;
		}
	}

	iptr[0] = 0;
	for (i=0; i<nb; i++)
		iptr[i+1] = iptr[i] + col[i];

	for (i=0; i<ig[nb] + nb; i++)
		aelem[i] = 0;

	for (i=0; i<nb; i++)
		col[i] = iptr[i]; //     

	for (i=0; i<nb; i++)
	{
		jptr[col[i]] = i;
		aelem[col[i]] = di[i];
		col[i]++;
	}

	for (i=0; i<nb; i++)
	{
		for (j=ig[i]; j<=ig[i+1]-1; j++)
		{
			k = jg[j];
			jptr[col[k]] = i;

			aelem[col[k]] = gg[j];

			col[k]++;
		}
	}	
}
void FormatConverter::From2x2ToCSRComplex_2_Sym(int nb, int *ig, int *jg, int *idi, int *ijg,double *di_block, double *ggl_block,
											MKL_INT64 *iptr, MKL_INT64 *jptr, double *aelem)
{
	int sz;
	int i, j, k, m;
	vector<MKL_INT64> col; 

	col.resize(nb, 0);

	for (i=0; i<nb; i++)
	{
		col[i] += 1; // 

		for (j=ig[i]; j<=ig[i+1]-1; j++)
		{
			k = jg[j];
			col[k]++;
		}
	}

	iptr[0] = 0;
	for (i=0; i<nb; i++)
		iptr[i+1] = iptr[i] + col[i];

	for (i=0; i<(ig[nb] + nb)*2; i++)
		aelem[i] = 0;

	for (i=0; i<nb; i++)
		col[i] = iptr[i]; //     

	for (i=0; i<nb; i++)
	{
		jptr[col[i]] = i;

		sz = idi[i+1] - idi[i];

		for (k=0; k<sz; k++)
			aelem[col[i]*2 + k] = di_block[idi[i] + k];

		col[i]++;
	}

	for (i=0; i<nb; i++)
	{
		for (j=ig[i]; j<=ig[i+1]-1; j++)
		{
			k = jg[j];
			jptr[col[k]] = i;
			sz = ijg[j+1] - ijg[j];

			for (m=0; m<sz; m++)
				aelem[col[k]*2 + m] = ggl_block[ijg[j] + m];

			col[k]++;
		}
	}	
}
void FormatConverter::From2x2ToCSRComplex_2_NS(int nb, int *ig, int *jg, int *idi, int *ijg,double *di_block, double *ggl_block, double *ggu_block,
											   MKL_INT64 *iptr, MKL_INT64 *jptr, double *aelem)
{
	int sz;
	int i, j, k, m;
	vector<int> col; 

	col.resize(nb, 0);

	for (i=0; i<nb; i++)
	{
		col[i] += ig[i+1] - ig[i] + 1; //   + 

		for (j=ig[i]; j<=ig[i+1]-1; j++)
		{
			k = jg[j];
			col[k]++;
		}
	}

	iptr[0] = 0;
	for (i=0; i<nb; i++)
		iptr[i+1] = iptr[i] + col[i];

	for (i=0; i<(ig[nb]*2 + nb)*2; i++)
		aelem[i] = 0;

	for (i=0; i<nb; i++)
		col[i] = iptr[i]; //     

	for (i=0; i<nb; i++)
	{
		for (j=ig[i]; j<=ig[i+1]-1; j++)
		{
			jptr[col[i]] = jg[j];

			sz = ijg[j+1] - ijg[j];

			for (k=0; k<sz; k++)
				aelem[col[i]*2 + k] = ggl_block[ijg[j] + k];

			col[i]++;
		}

		jptr[col[i]] = i;

		sz = idi[i+1] - idi[i];

		for (k=0; k<sz; k++)
			aelem[col[i]*2 + k] = di_block[idi[i] + k];

		col[i]++;
	}

	for (i=0; i<nb; i++)
	{
		for (j=ig[i]; j<=ig[i+1]-1; j++)
		{
			k = jg[j];
			jptr[col[k]] = i;
			sz = ijg[j+1] - ijg[j];

			for (m=0; m<sz; m++)
				aelem[col[k]*2 + m] = ggu_block[ijg[j] + m];

			col[k]++;
		}
	}	
}
