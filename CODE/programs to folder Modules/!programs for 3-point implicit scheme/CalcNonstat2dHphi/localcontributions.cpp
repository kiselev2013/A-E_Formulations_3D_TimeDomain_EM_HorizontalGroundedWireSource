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
 *  This file contains the code for calculating 2D local matirxes
 *  
 *  
 *  Written by Prof. Yuri G. Soloveichik, Prof. Marina G. Persova, Ph.D. Denis V. Vagin  
 *  Novosibirsk State Technical University,                                              
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                    
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                     
*/

#include "stdafx.h"
#include <math.h>
#include "TaskSP.h"


/*!     rz */
int GetStiffnessMatrix(const double &r0, const double &z0, // coordinates of local node 2
					   const double &r3, const double &z3, // coordinates of local node 1
					   double m[4][4])
{

	double h1=r3-r0, h2=z0-z3, rk=r0;
	double r11, r12, r22, z11, z12, z22;
	double rp11, rp12, rp22;

	r11=(h1/3)*(rk+h1/4); r12=(h1/6)*(rk+h1/2); r22=(h1/3)*(rk+3*h1/4);
	z11=h2/3; z12=h2/6; z22=h2/3;

	rp11=(1+rk/h1)*(1+rk/h1)*log(1+h1/rk)-rk/h1-1.5; 
	rp12=-(1+rk/h1)*(rk/h1)*log(1+h1/rk)+rk/h1+0.5; 
	rp22=0.5-rk/h1+(rk/h1)*(rk/h1)*log(1+h1/rk);



	m[2][2]=		(1./h1)*(rk+h1/2)*z11+(1./h2)*r11	+	rp11*z11;
	m[2][3]=m[3][2]=-(1./h1)*(rk+h1/2)*z11+(1./h2)*r12	+	rp12*z11;
	m[2][0]=m[0][2]=(1./h1)*(rk+h1/2)*z12-(1./h2)*r11	+	rp11*z12;
	m[2][1]=m[1][2]=-(1./h1)*(rk+h1/2)*z12-(1./h2)*r12	+	rp12*z12;
	
	m[3][3]=		(1./h1)*(rk+h1/2)*z11+(1./h2)*r22	+	rp22*z11;
	m[3][0]=m[0][3]=-(1./h1)*(rk+h1/2)*z12-(1./h2)*r12	+	rp12*z12;
	m[3][1]=m[1][3]=(1./h1)*(rk+h1/2)*z12-(1./h2)*r22	+	rp22*z12;
	
	m[0][0]=		(1./h1)*(rk+h1/2)*z22+(1./h2)*r11	+	rp11*z22;
	m[0][1]=m[1][0]=-(1./h1)*(rk+h1/2)*z22+(1./h2)*r12	+	rp12*z22;
	
	m[1][1]=		(1./h1)*(rk+h1/2)*z22+(1./h2)*r22	+	rp22*z22;


	return RETCODE_OK;
}

/*!     rz */
int GetMassMatrix(const double &r0, const double &z0, // coordinates of local node 2
				   const double &r3, const double &z3, // coordinates of local node 1
				   double m[4][4])
{

	double h1=r3-r0, h2=z0-z3, rk=r0;
	double r11, r12, r22, z11, z12, z22;
	
	r11=(h1/3)*(rk+h1/4); r12=(h1/6)*(rk+h1/2); r22=(h1/3)*(rk+3*h1/4);
	z11=h2/3; z12=h2/6; z22=h2/3;

	


	m[2][2]=		r11*z11;
	m[2][3]=m[3][2]=r12*z11;
	m[2][0]=m[0][2]=r11*z12;
	m[2][1]=m[1][2]=r12*z12;

	m[3][3]=		r22*z11;
	m[3][0]=m[0][3]=r12*z12;
	m[3][1]=m[1][3]=r22*z12;

	m[0][0]=		r11*z22;
	m[0][1]=m[1][0]=r12*z22;

	m[1][1]=		r22*z22;

	return RETCODE_OK; 
}

/*!     rz  1/r^2 */
int GetStiffnessMatrix0(const double &r0, const double &z0, // coordinates of local node 2
					   const double &r3, const double &z3, // coordinates of local node 1
					   double m[4][4])
{

	double h1=r3-r0, h2=z0-z3, rk=r0;
	double r11, r12, r22, z11, z12, z22;
	
	r11=(h1/3)*(rk+h1/4); r12=(h1/6)*(rk+h1/2); r22=(h1/3)*(rk+3*h1/4);
	z11=h2/3; z12=h2/6; z22=h2/3;




	m[2][2]=		(1./h1)*(rk+h1/2)*z11+(1./h2)*r11;
	m[2][3]=m[3][2]=-(1./h1)*(rk+h1/2)*z11+(1./h2)*r12;
	m[2][0]=m[0][2]=(1./h1)*(rk+h1/2)*z12-(1./h2)*r11;
	m[2][1]=m[1][2]=-(1./h1)*(rk+h1/2)*z12-(1./h2)*r12;

	m[3][3]=		(1./h1)*(rk+h1/2)*z11+(1./h2)*r22;
	m[3][0]=m[0][3]=-(1./h1)*(rk+h1/2)*z12-(1./h2)*r12;
	m[3][1]=m[1][3]=(1./h1)*(rk+h1/2)*z12-(1./h2)*r22;

	m[0][0]=		(1./h1)*(rk+h1/2)*z22+(1./h2)*r11;
	m[0][1]=m[1][0]=-(1./h1)*(rk+h1/2)*z22+(1./h2)*r12;

	m[1][1]=		(1./h1)*(rk+h1/2)*z22+(1./h2)*r22;


	return RETCODE_OK;
}

/*!    d(Phi)/dr * Phi */
int GetMatrixDr(const double &r0, const double &z0, // coordinates of local node 2
			   const double &r3, const double &z3, // coordinates of local node 1
			   double m[4][4])
{
	double hrk=r3-r0, hzk=z0-z3, rk=r0;
	double i1=(3*rk+hrk)*hzk/18., i2=(3*rk+2*hrk)*hzk/18.;
	
	m[0][0]=    -i1;
	m[0][1]=     i1;
	m[0][2]=-0.5*i1;
	m[0][3]= 0.5*i1;

	m[1][0]=    -i2;
	m[1][1]=     i2;
	m[1][2]=-0.5*i2;
	m[1][3]= 0.5*i2;

	m[2][0]=-0.5*i1;
	m[2][1]= 0.5*i1;
	m[2][2]=    -i1;
	m[2][3]=     i1;

	m[3][0]=-0.5*i2;
	m[3][1]= 0.5*i2;
	m[3][2]=    -i2;
	m[3][3]=     i2;

	return RETCODE_OK; 
}

/*!    d(Phi)/dz * Phi */
int GetMatrixDz(const double &r0, const double &z0, // coordinates of local node 2
			   const double &r3, const double &z3, // coordinates of local node 1
			   double m[4][4])
{
	double hrk=r3-r0, hzk=z0-z3, rk=r0;
	double i1=hrk*(4*rk+hrk)/24., i2=hrk*(2*rk+hrk)/24., i3=hrk*(4*rk+3*hrk)/24.;
	
	m[0][0]=m[2][0]=-i1;
	m[0][1]=m[2][1]=-i2;
	m[0][2]=m[2][2]= i1;
	m[0][3]=m[2][3]= i2;

	m[1][0]=m[3][0]=-i2;
	m[1][1]=m[3][1]=-i3;
	m[1][2]=m[3][2]= i2;
	m[1][3]=m[3][3]= i3;

	return RETCODE_OK; 
}

/*!   (dPi/dz)*(dPj/dz)   rz */
int GetStiffnessMatrixDz(const double &r0, const double &z0, // coordinates of local node 2
					   const double &r3, const double &z3, // coordinates of local node 1
					   double m[4][4])
{

	double h1=r3-r0, h2=z0-z3, rk=r0;
	double r11, r12, r22, z11, z12, z22;
	double rp11, rp12, rp22;

	r11=(h1/3)*(rk+h1/4); r12=(h1/6)*(rk+h1/2); r22=(h1/3)*(rk+3*h1/4);
	z11=h2/3; z12=h2/6; z22=h2/3;






	m[2][2]=(1./h2)*r11;
	m[2][3]=m[3][2]=(1./h2)*r12;
	m[2][0]=m[0][2]=-(1./h2)*r11;
	m[2][1]=m[1][2]=-(1./h2)*r12;
	
	m[3][3]=(1./h2)*r22;
	m[3][0]=m[0][3]=-(1./h2)*r12;
	m[3][1]=m[1][3]=-(1./h2)*r22;
	
	m[0][0]=(1./h2)*r11;
	m[0][1]=m[1][0]=(1./h2)*r12;
	
	m[1][1]=(1./h2)*r22;

	return RETCODE_OK;
}

/*!  (dPi/dr)*(dPj/dr) + ...   rz */
int GetStiffnessMatrixDr(const double &r0, const double &z0, // coordinates of local node 2
					   const double &r3, const double &z3, // coordinates of local node 1
					   double m[4][4])
{

	double h1=r3-r0, h2=z0-z3, rk=r0;
	double r11, r12, r22, z11, z12, z22;
	double rp11, rp12, rp22;

	r11=(h1/3)*(rk+h1/4); r12=(h1/6)*(rk+h1/2); r22=(h1/3)*(rk+3*h1/4);
	z11=h2/3; z12=h2/6; z22=h2/3;

	rp11=(1+rk/h1)*(1+rk/h1)*log(1+h1/rk)-rk/h1-1.5; 
	rp12=-(1+rk/h1)*(rk/h1)*log(1+h1/rk)+rk/h1+0.5; 
	rp22=0.5-rk/h1+(rk/h1)*(rk/h1)*log(1+h1/rk);



	m[2][2]=		(1./h1)*(rk+h1/2)*z11	+	rp11*z11;
	m[2][3]=m[3][2]=-(1./h1)*(rk+h1/2)*z11	+	rp12*z11;
	m[2][0]=m[0][2]=(1./h1)*(rk+h1/2)*z12	+	rp11*z12;
	m[2][1]=m[1][2]=-(1./h1)*(rk+h1/2)*z12	+	rp12*z12;
	
	m[3][3]=		(1./h1)*(rk+h1/2)*z11	+	rp22*z11;
	m[3][0]=m[0][3]=-(1./h1)*(rk+h1/2)*z12	+	rp12*z12;
	m[3][1]=m[1][3]=(1./h1)*(rk+h1/2)*z12	+	rp22*z12;
	
	m[0][0]=		(1./h1)*(rk+h1/2)*z22	+	rp11*z22;
	m[0][1]=m[1][0]=-(1./h1)*(rk+h1/2)*z22	+	rp12*z22;
	
	m[1][1]=		(1./h1)*(rk+h1/2)*z22	+	rp22*z22;

	return RETCODE_OK;
}

