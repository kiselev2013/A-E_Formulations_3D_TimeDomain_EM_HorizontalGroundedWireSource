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
 *  This file contains the code for building portrait of SLAE matrix for 3D VFEM
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin                         
 *  Novosibirsk State Technical University,                 
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia       
 *  vdv_wk@mail.ru                                          
 *  Version 2.0 January 16, 2023                            
*/

#include "stdafx.h" 
#include "T_Portrait.h"
#include "Portret.h"
extern ofstream logfile;
extern void Memory_allocation_error(const char *var, const char *func);
T_Portrait::T_Portrait(int *ed, int n, int n_c, int n_elem)
{
	this->n = n;
	this->n_c = n_c;
	this->n_elem = n_elem;
	this->ed = (int(*)[25])ed;

	ig = NULL;
	jg = NULL;
	idi = NULL;
	ijg = NULL;
}
T_Portrait::~T_Portrait()
{
	if(ig) {delete [] ig; ig=NULL;}
	if(jg) {delete [] jg; jg=NULL;}
	if(idi) { delete [] idi; idi=NULL; }
	if(ijg) { delete [] ijg; ijg=NULL; }
}
void T_Portrait::Gen_Portrait()
{
	Portret P;
	int i, j, k;
	int tmp1, tmp2, e;
	int loc_size;
	_list *l=NULL, *s=NULL;
	std::vector<int> local, loc;

	logfile << "T_Portrait... ";

	s = new _list[n_c]; //    
	if(s==NULL)
		Memory_allocation_error("s", "T_Portrait::Gen_Portrait");

	for(i=0; i<n_c; i++)
		s[i].next = NULL;

	for(i=0; i<this->n_elem; i++) //  
	{
		for(j=0; j<12; j++) //   -   local   
		{                   //  <n_c,     -   -
			e = ed[i][j];
			local.push_back(e);
		}// j

		std::sort(local.begin(),local.end());
		std::unique_copy(local.begin(), local.end(), back_inserter(loc));

		loc_size = (int)loc.size();

		for(j=0; j<loc_size; j++) //    
		for(k=0; k<loc_size; k++)
		{
			tmp1 = loc[k];
			tmp2 = loc[j];
			if(tmp1 < tmp2) //   
			{
				P.Add_In_Ordered_List(&s[tmp2], tmp1);
			}
		}

		tmp1 = (int)local.size();// 
		for(j=0; j<tmp1; j++) local.pop_back();
		for(j=0; j<loc_size; j++) loc.pop_back();
	}// i

	this->size_jg = 0;
	for(i=0; i<this->n_c; i++)
		this->size_jg += P.Size_Of_List(s[i].next);

	ig = new int[this->n_c + 1];
	if(ig == 0)
		Memory_allocation_error("ig", "T_Portrait::Gen_Portrait");

	jg = new int[this->size_jg];
	if(jg == 0)
		Memory_allocation_error("jg", "T_Portrait::Gen_Portrait");

	i = 0;
	ig[0] = 0;
	for(k=0;k<this->n_c;k++)
	{
		l = s[k].next;
		while(l != NULL)
		{
			jg[i] = l->number;
			i++;
			l = l->next;
		}
		P.Clear_List(s[k].next);
		ig[k+1] = i;
	}

	if(s) {delete [] s; s=NULL;}

	logfile << "done.\n";
}
void T_Portrait::Gen_idi_ijg(int *nvkat, int (*nver)[14])
{
	int i, j, k, m, it, jt, ii, jj, i_mu, j_nu;
	int i2, j2;
	int t, t2;
	int a[12][12]; //       

	if((idi = new int[n_c+1])==0) Memory_allocation_error("idi", "T_Portrait::Gen_idi_ijg");
	if((ijg = new int[size_jg+1])==0) Memory_allocation_error("ijg", "T_Portrait::Gen_idi_ijg");

	for(i=0; i<n_c+1; i++)
		idi[i] = 0;

	for(i=0; i<size_jg+1; i++)
		ijg[i] = 0;


	for (i=0; i<n_elem; i++)
	{

		if(nvkat[i]==0)//   sigma=0 =>  -  =0
		{
			for(i2=0; i2<12; i2++)
				for(j2=0; j2<12; j2++)
					a[i2][j2] = 1;
		}
		else
		{
			for(i2=0; i2<12; i2++)
				for(j2=0; j2<12; j2++)
					a[i2][j2] = FILTER_MASS_MATRIX_VEC[i2][j2];
		}

		for(j=0; j<12; j++)
		{
			ii = ed[i][j];

			if(ii < n_c)
				Set_type_of_block(idi, ii, a[j][j]);

			for(k=0; k<12; k++)
			{
				jj = ed[i][k];
				if(jj < ii) 
				{
					for(m=ig[ii]; m<=ig[ii+1]-1; m++)
						if(jg[m]==jj)
						{
							Set_type_of_block(ijg, m, a[j][k]);
							break;
						}
				}
			}// k
		}// j
	}// i


	t = idi[0];
	idi[0] = 0;

	for (i=0; i<n_c; i++)
	{
		t2 = idi[i] + t;
		t = idi[i+1];
		idi[i+1] = t2;
	}

	t = ijg[0];
	ijg[0] = 0;

	for (i=0; i<size_jg; i++)
	{
		t2 = ijg[i] + t;
		t = ijg[i+1];
		ijg[i+1] = t2;
	}

}
void T_Portrait::Set_type_of_block(int *target_block, int adr, int type)
{
	if(target_block[adr]!=2)
		target_block[adr] = type;
}
