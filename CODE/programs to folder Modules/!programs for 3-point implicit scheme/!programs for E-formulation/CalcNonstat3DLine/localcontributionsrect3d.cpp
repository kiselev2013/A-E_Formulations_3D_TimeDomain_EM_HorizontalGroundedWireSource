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
 *  This file contains the code for calculating 3D local matirxes
 *  
 *  
 *  Written by Prof. Yuri G. Soloveichik, Prof. Marina G. Persova, Ph.D. Denis V. Vagin  
 *  Novosibirsk State Technical University,                                              
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                    
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                     
 *  Version 2.0 January 16, 2023                                                                          
*/

#include "stdafx.h"
#include <stdlib.h>

#include "TaskSP.h"


inline int indX(const int& i)
{
	return div(i, 2).rem;
}

inline int indY(const int& i)
{
	return div(div(i, 2).quot, 2).rem;
}

inline int indZ(const int& i)
{
	return div(i, 4).quot;
}


double G1d[2][2]={{1, -1}, 
				  {-1, 1}};
double M1d[2][2]={{2, 1}, 
				  {1, 2}};


/*!     */
int GetStiffnessMatrixForRect3D(
	const double &x1, const double &y1, const double &z1, // coordinates of local node 0
	const double &x2, const double &y2, const double &z2, // coordinates of local node 7
	double m[8][8])
{
	double 
		hx=x2-x1, hx1=1/hx,
		hy=y2-y1, hy1=1/hy,
		hz=z2-z1, hz1=1/hz;
	int i, j, ix, iy, iz, jx, jy, jz;
	for (i=0; i<8; i++)
	{
		ix=indX(i);
		iy=indY(i);
		iz=indZ(i);
		for (j=0; j<8; j++)
		{
			jx=indX(j);
			jy=indY(j);
			jz=indZ(j);
			m[i][j]=(G1d[ix][jx]*M1d[iy][jy]*M1d[iz][jz]*hx1*hy*hz+
					 M1d[ix][jx]*G1d[iy][jy]*M1d[iz][jz]*hx*hy1*hz+
					 M1d[ix][jx]*M1d[iy][jy]*G1d[iz][jz]*hx*hy*hz1)/36;
		}
	}
	return RETCODE_OK;
}

/*!     */
int GetMassMatrixForRect3D(
	const double &x1, const double &y1, const double &z1, // coordinates of local node 0
	const double &x2, const double &y2, const double &z2, // coordinates of local node 7
	double m[8][8])
{
	double 
		hx=x2-x1,
		hy=y2-y1,
		hz=z2-z1;
	int i, j, ix, iy, iz, jx, jy, jz;
	for (i=0; i<8; i++)
	{
		ix=indX(i);
		iy=indY(i);
		iz=indZ(i);
		for (j=0; j<8; j++)
		{
			jx=indX(j);
			jy=indY(j);
			jz=indZ(j);
			m[i][j]=M1d[ix][jx]*M1d[iy][jy]*M1d[iz][jz]*hx*hy*hz/216;
		}
	}
	return RETCODE_OK;
}

/*!     (Telma) */
void StiffnessMatrix(double Matr[8][8], double hx, double hy, double hz)
{			
	Matr[0][0] = 1./9*(hy*hz/hx+hx*hz/hy+hx*hy/hz);
	Matr[0][1] = -hy*hz/(9*hx)+hx*hz/(18*hy)+hx*hy/(18*hz);
	Matr[0][2] = hy*hz/(18*hx)-hx*hz/(9*hy)+hx*hy/(18*hz);
	Matr[0][3] = -hy*hz/(18*hx)-hx*hz/(18*hy)+hx*hy/(36*hz);
	Matr[0][4] = hy*hz/(18*hx)+hx*hz/(18*hy)-hx*hy/(9*hz);
	Matr[0][5] = -hy*hz/(18*hx)+hx*hz/(36*hy)-hx*hy/(18*hz);
	Matr[0][6] = hy*hz/(36*hx)-hx*hz/(18*hy)-hx*hy/(18*hz);
	Matr[0][7] = -hy*hz/(36*hx)-hx*hz/(36*hy)-hx*hy/(36*hz);

	Matr[1][1] = 1./9*(hy*hz/hx+hx*hz/hy+hx*hy/hz);
	Matr[1][2] = -hy*hz/(18*hx)-hx*hz/(18*hy)+hx*hy/(36*hz);
	Matr[1][3] = hy*hz/(18*hx)-hx*hz/(9*hy)+hx*hy/(18*hz);
	Matr[1][4] = -hy*hz/(18*hx)+hx*hz/(36*hy)-hx*hy/(18*hz);
	Matr[1][5] = hy*hz/(18*hx)+hx*hz/(18*hy)-hx*hy/(9*hz);
	Matr[1][6] = -hy*hz/(36*hx)-hx*hz/(36*hy)-hx*hy/(36*hz);
	Matr[1][7] = hy*hz/(36*hx)-hx*hz/(18*hy)-hx*hy/(18*hz);

	Matr[2][2] = 1./9*(hy*hz/hx+hx*hz/hy+hx*hy/hz);
	Matr[2][3] = -hy*hz/(9*hx)+hx*hz/(18*hy)+hx*hy/(18*hz); 
	Matr[2][4] = hy*hz/(36*hx)-hx*hz/(18*hy)-hx*hy/(18*hz);
	Matr[2][5] = -hy*hz/(36*hx)-hx*hz/(36*hy)-hx*hy/(36*hz);
	Matr[2][6] = hy*hz/(18*hx)+hx*hz/(18*hy)-hx*hy/(9*hz);
	Matr[2][7] = -hy*hz/(18*hx)+hx*hz/(36*hy)-hx*hy/(18*hz);

	Matr[3][3] = 1./9*(hy*hz/hx+hx*hz/hy+hx*hy/hz);
	Matr[3][4] = -hy*hz/(36*hx)-hx*hz/(36*hy)-hx*hy/(36*hz);
	Matr[3][5] = hy*hz/(36*hx)-hx*hz/(18*hy)-hx*hy/(18*hz);
	Matr[3][6] = -hy*hz/(18*hx)+hx*hz/(36*hy)-hx*hy/(18*hz);
	Matr[3][7] = hy*hz/(18*hx)+hx*hz/(18*hy)-hx*hy/(9*hz);

	Matr[4][4] = 1./9*(hy*hz/hx+hx*hz/hy+hx*hy/hz);
	Matr[4][5] = -hy*hz/(9*hx)+hx*hz/(18*hy)+hx*hy/(18*hz);
	Matr[4][6] = hy*hz/(18*hx)-hx*hz/(9*hy)+hx*hy/(18*hz);
	Matr[4][7] = -hy*hz/(18*hx)-hx*hz/(18*hy)+hx*hy/(36*hz);

	Matr[5][5] = 1./9*(hy*hz/hx+hx*hz/hy+hx*hy/hz);
	Matr[5][6] = -hy*hz/(18*hx)-hx*hz/(18*hy)+hx*hy/(36*hz);
	Matr[5][7] = hy*hz/(18*hx)-hx*hz/(9*hy)+hx*hy/(18*hz);

	Matr[6][6] = 1./9*(hy*hz/hx+hx*hz/hy+hx*hy/hz);
	Matr[6][7] = -hy*hz/(9*hx)+hx*hz/(18*hy)+hx*hy/(18*hz);

	Matr[7][7] = 1./9*(hy*hz/hx+hx*hz/hy+hx*hy/hz);


	for (int i=0; i<8; i++)
	{

		for (int j=i+1; j<8; j++)
			Matr[j][i] = Matr[i][j];
	}
}

/*!     (Telma) */
void MassMatrix(double Matr[8][8], double hx, double hy, double hz)
{			

	double hxhyhz = hx*hy*hz;

	Matr[0][0] = hxhyhz*(1./27);
	Matr[0][1] = hxhyhz*(1./54);
	Matr[0][2] = hxhyhz*(1./54);
	Matr[0][3] = hxhyhz*(1./108);
	Matr[0][4] = hxhyhz*(1./54);
	Matr[0][5] = hxhyhz*(1./108);
	Matr[0][6] = hxhyhz*(1./108);
	Matr[0][7] = hxhyhz*(1./216);

	Matr[1][1] = hxhyhz*(1./27);
	Matr[1][2] = hxhyhz*(1./108);
	Matr[1][3] = hxhyhz*(1./54);
	Matr[1][4] = hxhyhz*(1./108);
	Matr[1][5] = hxhyhz*(1./54);
	Matr[1][6] = hxhyhz*(1./216);
	Matr[1][7] = hxhyhz*(1./108);

	Matr[2][2] = hxhyhz*(1./27);
	Matr[2][3] = hxhyhz*(1./54);
	Matr[2][4] = hxhyhz*(1./108);
	Matr[2][5] = hxhyhz*(1./216);
	Matr[2][6] = hxhyhz*(1./54);
	Matr[2][7] = hxhyhz*(1./108);

	Matr[3][3] = hxhyhz*(1./27);
	Matr[3][4] = hxhyhz*(1./216);
	Matr[3][5] = hxhyhz*(1./108);
	Matr[3][6] = hxhyhz*(1./108);
	Matr[3][7] = hxhyhz*(1./54);

	Matr[4][4] = hxhyhz*(1./27);
	Matr[4][5] = hxhyhz*(1./54);
	Matr[4][6] = hxhyhz*(1./54);
	Matr[4][7] = hxhyhz*(1./108);

	Matr[5][5] = hxhyhz*(1./27);
	Matr[5][6] = hxhyhz*(1./108);
	Matr[5][7] = hxhyhz*(1./54);

	Matr[6][6] = hxhyhz*(1./27);
	Matr[6][7] = hxhyhz*(1./54);

	Matr[7][7] = hxhyhz*(1./27);


	for (int i=0; i<8; i++)
	{
		for (int j=i+1; j<8; j++)
			Matr[j][i] = Matr[i][j];

	}
}

/*!    grad(Phi) * (Phi) */
void GradMatrix(double Dx[8][8], double Dy[8][8], double Dz[8][8], double hx, double hy, double hz)
{			
	double h1,h2,h3; 
	h1=hy*hz/72; h2=hx*hz/72; h3=hy*hx/72;
	Dx[0][0]=-4*h1; Dx[0][1]=-4*h1; Dx[0][2]=-2*h1;Dx[0][3]=-2*h1; 
	Dx[0][4]=-2*h1; Dx[0][5]=-2*h1; Dx[0][6]=-h1; Dx[0][7]=-h1; 

	int i;
	for (i=0; i<8; i++)
		Dx[1][i] = - Dx[0][i];


	Dx[2][0]=-2*h1; Dx[2][1]=-2*h1; Dx[2][2]=-4*h1;Dx[2][3]=-4*h1; 
	Dx[2][4]=-h1; Dx[2][5]=-h1; Dx[2][6]=-2*h1; Dx[2][7]=-2*h1; 

	for ( i=0; i<8; i++)
		Dx[3][i] = - Dx[2][i];


	Dx[4][0]=-2*h1; Dx[4][1]=-2*h1; Dx[4][2]=-h1;Dx[4][3]=-h1; 
	Dx[4][4]=-4*h1; Dx[4][5]=-4*h1; Dx[4][6]=-2*h1; Dx[4][7]=-2*h1; 
	for ( i=0; i<8; i++)
		Dx[5][i] = - Dx[4][i];

	Dx[6][0]=-h1; Dx[6][1]=-h1; Dx[6][2]=-2*h1; Dx[6][3]=-2*h1; 
	Dx[6][4]=-2*h1; Dx[6][5]=-2*h1; Dx[6][6]=-4*h1; Dx[6][7]=-4*h1; 
	for ( i=0; i<8; i++)
		Dx[7][i] = - Dx[6][i];


	Dy[0][0]=-4*h2; Dy[0][1]=-2*h2; Dy[0][2]=-4*h2;Dy[0][3]=-2*h2; 
	Dy[0][4]=-2*h2; Dy[0][5]=-h2; Dy[0][6]=-2*h2; Dy[0][7]=-h2; 

	Dy[1][0]=-2*h2; Dy[1][1]=-4*h2; Dy[1][2]=-2*h2;Dy[1][3]=-4*h2; 
	Dy[1][4]=-h2; Dy[1][5]=-2*h2; Dy[1][6]=-h2; Dy[1][7]=-2*h2; 

	for ( i=0; i<8; i++)
	{  Dy[2][i] = - Dy[0][i];
	Dy[3][i] = - Dy[1][i];
	}


	Dy[4][0]=-2*h2; Dy[4][1]=-h2; Dy[4][2]=-2*h2;Dy[4][3]=-h2; 
	Dy[4][4]=-4*h2; Dy[4][5]=-2*h2; Dy[4][6]=-4*h2; Dy[4][7]=-2*h2; 

	Dy[5][0]=-h2; Dy[5][1]=-2*h2; Dy[5][2]=-h2;Dy[5][3]=-2*h2; 
	Dy[5][4]=-2*h2; Dy[5][5]=-4*h2; Dy[5][6]=-2*h2; Dy[5][7]=-4*h2; 

	for ( i=0; i<8; i++)
	{		Dy[6][i] = - Dy[4][i];
	Dy[7][i] = - Dy[5][i];
	}

	Dz[0][0]=-4*h3; Dz[0][1]=-2*h3; Dz[0][2]=-2*h3; Dz[0][3]=-h3; 
	Dz[0][4]=-4*h3; Dz[0][5]=-2*h3; Dz[0][6]=-2*h3; Dz[0][7]=-h3; 

	Dz[1][0]=-2*h3; Dz[1][1]=-4*h3; Dz[1][2]=-h3; Dz[1][3]=-2*h3; 
	Dz[1][4]=-2*h3; Dz[1][5]=-4*h3; Dz[1][6]=-h3; Dz[1][7]=-2*h3; 

	Dz[2][0]=-2*h3; Dz[2][1]=-h3; Dz[2][2]=-4*h3; Dz[2][3]=-2*h3; 
	Dz[2][4]=-2*h3; Dz[2][5]=-h3; Dz[2][6]=-4*h3; Dz[2][7]=-2*h3; 

	Dz[3][0]=-h3; Dz[3][1]=-2*h3; Dz[3][2]=-2*h3; Dz[3][3]=-4*h3; 
	Dz[3][4]=-h3; Dz[3][5]=-2*h3; Dz[3][6]=-2*h3; Dz[3][7]=-4*h3; 

	for ( i=0; i<8; i++)
	{		Dz[4][i] = - Dz[0][i];
	Dz[5][i] = - Dz[1][i];
	Dz[6][i] = - Dz[2][i];
	Dz[7][i] = - Dz[3][i];
	}

}
