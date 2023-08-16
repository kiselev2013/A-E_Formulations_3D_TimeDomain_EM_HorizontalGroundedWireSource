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
 *  This file contains the code of functions for calculating the values of basic functions and their derivatives on an element
 *  
 *  
 *  Written by Prof. Yuri G. Soloveichik, Prof. Marina G. Persova                               
 *  Novosibirsk State Technical University,                                                                               
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                                                     
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)       
 *  Version 2.0 January 16, 2023                                                                                          
*/

#include "stdafx.h"
#include "T_Brick.h"
#include "gauss3.h"

extern void Mult_Plot(double *a, double *x, double *y, int n);
extern double Scal(const int &n, const double *v1, const double *v2);

double Scal(double *a, double *b, long n)
{
	long i;
	double sum = 0.0;

	for(i=0; i<n; i++)
		sum += a[i]*b[i];

	return sum;
}

T_Brick::T_Brick(double *x_coords, double *y_coords, double *z_coords)
{
	for (int i=0; i<8; i++)
	{
		this->x[i] = x_coords[i];
		this->y[i] = y_coords[i];
		this->z[i] = z_coords[i];
	}

	this->xk  = this->x[0];
	this->xk1 = this->x[7];
	this->yk  = this->y[0];
	this->yk1 = this->y[7];
	this->zk  = this->z[0];
	this->zk1 = this->z[7];

	this->hx = this->xk1 - this->xk;
	this->hy = this->yk1 - this->yk; 
	this->hz = this->zk1 - this->zk;
}
T_Brick::T_Brick(int num, int (*nver)[14], double (*xyz)[3])
{
	this->num = num;
	this->nver = nver;
	this->xyz = xyz;

	for (int i=0; i<8; i++)
	{
		this->x[i] = xyz[nver[num][i]][0];
		this->y[i] = xyz[nver[num][i]][1];
		this->z[i] = xyz[nver[num][i]][2];
	}

	this->xk  = this->x[0];
	this->xk1 = this->x[7];
	this->yk  = this->y[0];
	this->yk1 = this->y[7];
	this->zk  = this->z[0];
	this->zk1 = this->z[7];

	this->hx = this->xk1 - this->xk;
	this->hy = this->yk1 - this->yk; 
	this->hz = this->zk1 - this->zk;
}
T_Brick::~T_Brick()
{
}
void T_Brick::Calc_J(int n_of_point)
{
	Calc_J_in_parallelepiped();		
}
void T_Brick::Calc_J(double xi, double eta, double zeta)
{
	Calc_J_in_parallelepiped();
}

void T_Brick::Calc_V_Node_2(double *q,double xi, double eta, double zeta)
{
	int i;
	double dphi_x, dphi_y, dphi_z;

	for(i=0; i<3; i++)V[i] = 0.0;

	for(i=0; i<8; i++)
	{
		dphi_x = dPhi_node_2(i, 0, xi, eta, zeta);
		dphi_y = dPhi_node_2(i, 1, xi, eta, zeta);
		dphi_z = dPhi_node_2(i, 2, xi, eta, zeta);
		
		V[0] += q[i]*dphi_x; //  d_xi[i];
		V[1] += q[i]*dphi_y; //  d_eta[i];
		V[2] += q[i]*dphi_z; //  d_zeta[i];
	}
}

void T_Brick::Calc_J_Node_2(double xi, double eta, double zeta)
{
		J[0][0] = hx;
		J[1][1] = hy;
		J[2][2] = hz;
		J[0][1] = J[0][2] = J[1][0] = J[1][2] = J[2][0] = J[2][1] = 0.0;
		det_J = hx*hy*hz;
		det_J_abs = fabs(det_J);
		J_1_T[0][0] = J_1[0][0] = 1.0/hx;
		J_1_T[1][1] = J_1[1][1] = 1.0/hy;
		J_1_T[2][2] = J_1[2][2] = 1.0/hz;
		J_1_T[0][1] = J_1_T[0][2] = J_1_T[1][0] = J_1_T[1][2] = J_1_T[2][0] = J_1_T[2][1] =
		J_1[0][1] = J_1[0][2] = J_1[1][0] = J_1[1][2] = J_1[2][0] = J_1[2][1] =	0.0;
	}

void T_Brick::Calc_J_on_face(int n_of_point)
{
	Calc_J_in_parallelepiped();
}
double T_Brick::l0(double x)
{
	return (1.0 - x)*0.5;
}
double T_Brick::l1(double x)
{
	return (x + 1.0)*0.5;
}
double T_Brick::Phi_node(int i, double x, double y, double z)
{
	switch(i)
	{
	case 0:
		return l0(x)*l0(y)*l0(z);
	case 1:
		return l1(x)*l0(y)*l0(z);
	case 2:
		return l0(x)*l1(y)*l0(z);
	case 3:
		return l1(x)*l1(y)*l0(z);
	case 4:
		return l0(x)*l0(y)*l1(z);
	case 5:
		return l1(x)*l0(y)*l1(z);
	case 6:
		return l0(x)*l1(y)*l1(z);
	default:
		return l1(x)*l1(y)*l1(z);
	}
}

double T_Brick::dPhi_node_2(int i, int j, double xi, double eta, double zeta)
{
	switch(j)
	{
	case 0:
		switch(i)
		{
		case 0: return -(1.0-eta)*(1.0-zeta);
		case 1: return  (1.0-eta)*(1.0-zeta);
		case 2: return -(eta)*(1.0-zeta);
		case 3: return  (eta)*(1.0-zeta);
		case 4: return -(1.0-eta)*(zeta);
		case 5: return  (1.0-eta)*(zeta);
		case 6: return -(eta)*(zeta);
		default: return (eta)*(zeta);
		}
	case 1:
		switch(i)
		{
		case 0: return -(1.0-xi)*(1.0-zeta);
		case 1: return -(xi)*(1.0-zeta);
		case 2: return (1.0-xi)*(1.0-zeta);
		case 3: return (xi)*(1.0-zeta);
		case 4: return -(1.0-xi)*(zeta);
		case 5: return -(xi)*(zeta);
		case 6: return  (1.0-xi)*(zeta);
		default: return (xi)*(zeta);
		}
	default:
		switch(i)
		{
		case 0: return -(1.0-xi)*(1.0-eta);
		case 1: return -(xi)*(1.0-eta);
		case 2: return -(1.0-xi)*(eta);
		case 3: return -(xi)*(eta);
		case 4: return  (1.0-xi)*(1.0-eta);
		case 5: return  (xi)*(1.0-eta);
		case 6: return  (1.0-xi)*(eta);
		default: return (xi)*(eta);
		}
	}
}
double T_Brick::dPhi_node(int i, int j, double xi, double eta, double zeta)
{
	switch(j)
	{
	case 0:
		switch(i)
		{
		case 0: return -(-1.0+eta)*(-1.0+zeta)/8.0;
		case 1: return (-1.0+eta)*(-1.0+zeta)/8.0;
		case 2: return (eta+1.0)*(-1.0+zeta)/8.0;
		case 3: return -(eta+1.0)*(-1.0+zeta)/8.0;
		case 4: return (-1.0+eta)*(zeta+1.0)/8.0;
		case 5: return -(-1.0+eta)*(zeta+1.0)/8.0;
		case 6: return -(eta+1.0)*(zeta+1.0)/8.0;
		default: return (eta+1.0)*(zeta+1.0)/8.0;
		}
	case 1:
		switch(i)
		{
		case 0: return -(-1.0+xi)*(-1.0+zeta)/8.0;      
		case 1: return (xi+1.0)*(-1.0+zeta)/8.0;
		case 2: return (-1.0+xi)*(-1.0+zeta)/8.0;
		case 3: return -(xi+1.0)*(-1.0+zeta)/8.0;
		case 4: return (-1.0+xi)*(zeta+1.0)/8.0;
		case 5: return -(xi+1.0)*(zeta+1.0)/8.0;
		case 6: return -(-1.0+xi)*(zeta+1.0)/8.0;
		default: return (xi+1.0)*(zeta+1.0)/8.0;
		}
	default:
		switch(i)
		{
		case 0: return -(-1.0+xi)*(-1.0+eta)/8.0;      
		case 1: return (xi+1.0)*(-1.0+eta)/8.0;
		case 2: return (-1.0+xi)*(eta+1.0)/8.0;
		case 3: return -(xi+1.0)*(eta+1.0)/8.0;
		case 4: return (-1.0+xi)*(-1.0+eta)/8.0;
		case 5: return -(xi+1.0)*(-1.0+eta)/8.0;
		case 6: return -(-1.0+xi)*(eta+1.0)/8.0;
		default: return (xi+1.0)*(eta+1.0)/8.0;
		}
	}
}
void T_Brick::Mapping(double *in, double *out)
{
	int i;
	double temp;

	out[0] = out[1] = out[2] = 0.0;

	for(i=0; i<8; i++)
	{
		temp = Phi_node(i, in[0], in[1], in[2]);
		out[0] += x[i]*temp;
		out[1] += y[i]*temp;
		out[2] += z[i]*temp;
	}
}
void T_Brick::Calc_J_in_parallelepiped()
{

	J[0][0] = hx*0.5;
	J[1][1] = hy*0.5;
	J[2][2] = hz*0.5;

	J[0][1] = J[0][2] = J[1][0] = J[1][2] = J[2][0] = J[2][1] = 0.0;

	det_J = hx*hy*hz/8.0;

	det_J_abs = fabs(det_J);

	J_1_T[0][0] = J_1[0][0] = 2.0/hx;
	J_1_T[1][1] = J_1[1][1] = 2.0/hy;
	J_1_T[2][2] = J_1[2][2] = 2.0/hz;
	J_1_T[0][1] = J_1_T[0][2] = J_1_T[1][0] = J_1_T[1][2] = J_1_T[2][0] = J_1_T[2][1] =
		J_1[0][1] = J_1[0][2] = J_1[1][0] = J_1[1][2] = J_1[2][0] = J_1[2][1] =	0.0;
}
void T_Brick::Transformation_of_variables(const double *in, double *out)
{
	out[0] = Xi(in[0]);
	out[1] = Eta(in[1]);
	out[2] = Zeta(in[2]);
}
void T_Brick::Transformation_of_variables(double *x, double *y, double *z)
{
	x[0] = Xi(x[0]); 
	y[0] = Eta(y[0]); 
	z[0] = Zeta(z[0]); 
}
double T_Brick::Xi(double x)
{
	return (2.0*x - xk - xk1)/hx;
}
double T_Brick::Eta(double y)
{
	return (2.0*y - yk - yk1)/hy;
}
double T_Brick::Zeta(double z)
{
	return (2.0*z - zk - zk1)/hz;
}
double T_Brick::GetValueInHexCenter(double *q)
{
	return 0.125*(q[0]+q[1]+q[2]+q[3]+q[4]+q[5]+q[6]+q[7]);
}
void T_Brick::GetGradInHexCenter(double *q, double *out, int cj)
{
	int i,j;

	if(cj)Calc_J_Node_2(0.5,0.5,0.5);
	Calc_V_Node_2(q,0.5,0.5,0.5);
	out[0]=out[1]=out[2]=0;
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			out[i]+=J_1_T[i][j]*V[j];
		}
	}
}
double T_Brick::ScalarFieldOnPar(double x, double y, double z, double *ves)
{
	int i;
	double res = 0.0;
	double coordGlobal[3], coordLocal[3];

	coordGlobal[0] = x;
	coordGlobal[1] = y;
	coordGlobal[2] = z;
	Transformation_of_variables(coordGlobal, coordLocal);

	for (i=0; i<8; i++)
		res += Phi_node(i, coordLocal[0], coordLocal[1], coordLocal[2])*ves[i];

	return res;
}
double T_Brick::DxOfScalarFieldOnPar(double x, double y, double z, double *ves)
{
	int i;
	double res = 0.0;
	double coordGlobal[3], coordLocal[3];

	coordGlobal[0] = x;
	coordGlobal[1] = y;
	coordGlobal[2] = z;
	Transformation_of_variables(coordGlobal, coordLocal);

	for (i=0; i<8; i++)
		res += dPhi_node(i, 0, coordLocal[0], coordLocal[1], coordLocal[2])*ves[i];
	res = res*2.0/hx;
	
	return res;
}
double T_Brick::DyOfScalarFieldOnPar(double x, double y, double z, double *ves)
{
	int i;
	double res = 0.0;
	double coordGlobal[3], coordLocal[3];

	coordGlobal[0] = x;
	coordGlobal[1] = y;
	coordGlobal[2] = z;
	Transformation_of_variables(coordGlobal, coordLocal);

	for (i=0; i<8; i++)
		res += dPhi_node(i, 1, coordLocal[0], coordLocal[1], coordLocal[2])*ves[i];
	res = res*2.0/hy;

	return res;
}
double T_Brick::DzOfScalarFieldOnPar(double x, double y, double z, double *ves)
{
	int i;
	double res = 0.0;
	double coordGlobal[3], coordLocal[3];

	coordGlobal[0] = x;
	coordGlobal[1] = y;
	coordGlobal[2] = z;
	Transformation_of_variables(coordGlobal, coordLocal);

	for (i=0; i<8; i++)
		res += dPhi_node(i, 2, coordLocal[0], coordLocal[1], coordLocal[2])*ves[i];
	res = res*2.0/hz;

	return res;
}
