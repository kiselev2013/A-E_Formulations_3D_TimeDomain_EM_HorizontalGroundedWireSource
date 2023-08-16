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
 *  This file contains the code of function for generating a portrait of SLAE matrix for 3D scalar tasks
 *  
 *  
 *  Written by Prof. Yuri G. Soloveichik, Prof. Marina G. Persova, Ph.D. Denis V. Vagin 
 *  Novosibirsk State Technical University,                                             
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                   
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                    
 *  Version 2.0 January 16, 2023
*/

#include "stdafx.h"

using namespace std;

#include "Portret.h"

extern void CloseProgramm(int code);

void Portret::Add_In_Ordered_List(_list *s, int x) //  -   
{
	_list *head, *cur, *prev, *list_new;

	head = s->next;
	cur = head;
	prev = NULL;

	while(cur != NULL)
	{
   		if(x < cur->number)
      		break;
		prev = cur;
		cur = cur->next;
	}

	if(prev == NULL) //   
	{
	  	list_new = new _list;
		if(list_new == 0)
		{
			char var[] = {"list_new"};
			char func[] = {"Portret::Add_In_Ordered_List"};
			printf("Memory allocation error for variable \"%s\" in function \"%s\"\n",var,func);
			CloseProgramm(1);
		}

		list_new->number = x;

		list_new->next = head;
		s->next = list_new;
	}
	else if(prev->number!=x)//   
	{
	  	list_new = new _list;
		if(list_new == 0)
		{
			char var[] = {"list_new"};
			char func[] = {"Portret::Add_In_Ordered_List"};
			printf("Memory allocation error for variable \"%s\" in function \"%s\"\n",var,func);
			CloseProgramm(1);
		}

		list_new->number = x;

   		list_new->next = prev->next;
		prev->next = list_new;
	}
}
Portret::Portret(int (*nver)[14], int n_of_elements, int n_of_nodes, int n_c)
{
	this->nver = nver;
	this->n_of_elements = n_of_elements;
	this->is_mem_ig_jg_allocated = false;
	this->n = n_of_nodes;
	this->n_c = n_c;
}
Portret::Portret()
{
	this->is_mem_ig_jg_allocated = false;
}
void Portret::Clear_List(_list *s)
{
	_list *cur_pos, *next_pos;

	cur_pos = s;

	while(cur_pos != NULL)
	{
		next_pos = cur_pos->next;
		delete cur_pos;
		cur_pos = next_pos;
	}
}
Portret::~Portret()
{
	if(this->is_mem_ig_jg_allocated==true)
	{
		delete [] this->ig;
		delete [] this->jg;
	}
}
void Portret::Gen_T_Portrait2()
{
	Portret P;
	int i, j, k;
	int tmp1, tmp2, e;
	_list *l, *s;
	std::vector<int> local, loc;

	if ((s = new _list[n_c]) == 0) cout<<"s"<<"Gen_T_Portrait2"<<endl;

	for(i=0; i<n_c; i++)
		s[i].next = NULL;

	for(i=0; i<n_of_elements; i++) //  
	{
		for(j=0; j<8; j++) //   -   local   
		{                   //  <n_c,     -   -
			e = nver[i][j];
			local.push_back(e);
		}// j

		std::sort(local.begin(),local.end());
		std::unique_copy(local.begin(), local.end(), back_inserter(loc));

		for(j=0; j<(int)loc.size(); j++) //    
		{
			for(k=0; k<(int)loc.size(); k++)
			{
				tmp1 = loc[k];
				tmp2 = loc[j];

				if(tmp1 < tmp2) //   
					P.Add_In_Ordered_List(&s[tmp2], tmp1);					
			}
		}

		loc.clear();
		local.clear();
	}// i

	size_jg = 0;
	for(i=0; i<n_c; i++)
		size_jg += P.Size_Of_List(s[i].next);

	if ((ig = new int[n_c + 1]) == 0) cout<<"ig"<<"Gen_T_Portrait2"<<endl;
	if ((jg = new int[size_jg]) == 0) cout<<"jg"<<"Gen_T_Portrait2"<<endl;

	is_mem_ig_jg_allocated = true;

	i = 0;
	ig[0] = 0;
	for(k=0; k<n_c; k++)
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

	delete [] s;
}
int Portret::Size_Of_List(_list *s)
{
	int i;
	_list *cur_pos;

	i = 0;
	cur_pos = s;

	while(cur_pos != NULL)
	{
		i++;
		cur_pos = cur_pos->next;
	}

	return i;
}
