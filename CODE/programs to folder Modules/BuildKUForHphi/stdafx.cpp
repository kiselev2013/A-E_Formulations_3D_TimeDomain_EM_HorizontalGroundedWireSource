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
 *  This file contains code for utility functions
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 January 16, 2023
*/

#include "stdafx.h"

int NUMBEROFABDIPOLES;

int sig(const double &v)
{
	return ((v>0)? 1 : (v<0)? -1 : 0);
}

int FindIntervalInDoubleVec(vector<double> &vec,double elem)
{
	int i,j,k;
	i=0;
	k=(int)vec.size()-1;
	if(k==-1)return -1;
	if(elem<vec[0]-1e-6 || elem>vec[k]+1e-6)return -1;
	j=(i+k)/2;
	do{		
		if(elem>vec[j])i=j;
		else k=j;
		j=(i+k)/2;
	}while(j!=i);
	return i;
}

int FindIntervalInDoubleMas(double *vec,int size,double elem)
{
	int i,j,k;
	i=0;
	k=size-1;
	if(k==-1)return -1;
	if(elem<vec[0]-1e-6 || elem>vec[k]+1e-6)return -1;
	j=(i+k)/2;
	do{		
		if(elem>vec[j])i=j;
		else k=j;
		j=(i+k)/2;
	}while(j!=i);
	return i;
}

int FindIntervalInDoubleMas(double *vec, int beg, int end, double elem)
{
	int i, j, k;
	i = beg;
	k = end;
	if (k == -1)return -1;
	if (elem<vec[0] - 1e-6 || elem>vec[k] + 1e-6)return -1;
	j = (i + k) / 2;
	do{
		if (elem>vec[j])i = j;
		else k = j;
		j = (i + k) / 2;
	} while (j != i);
	return i;
}
