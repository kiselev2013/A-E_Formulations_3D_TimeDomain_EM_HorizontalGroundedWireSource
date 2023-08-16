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
 *  This file contains the code for work with 3D points
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin                                
 *  Novosibirsk State Technical University,                        
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia              
 *  Corresponding author: vdv_wk@mail.ru                           
 *  Version 2.0 January 16, 2023                                   
*/


#include "stdafx.h"
#include "PointXYZ.h"

PointXYZ::PointXYZ() 
{ 
	x=y=z=0; 
}
PointXYZ::PointXYZ(const double& _x, const double& _y, const double& _z) 
{
	x=_x;
	y=_y;
	z=_z; 
}
const PointXYZ& PointXYZ::operator+=(const PointXYZ& p)
{
	x+=p.x;
	y+=p.y;
	z+=p.z;
	return *this;
}
const PointXYZ& PointXYZ::operator-=(const PointXYZ& p)
{
	x-=p.x;
	y-=p.y;
	z-=p.z;
	return *this;
}


PointXYZ operator*(const PointXYZ& p, const double& a)
{
	return PointXYZ(p.x*a, p.y*a, p.z*a);
}

PointXYZ operator+(const PointXYZ& p1, const PointXYZ& p2)
{
	return PointXYZ(p1.x+p2.x, p1.y+p2.y, p1.z+p2.z);
}

ifstream& operator>>(ifstream& inf, PointXYZ& p)
{
	inf >> p.x >> p.y >> p.z;
	return inf;
}

ifstream& operator>(ifstream& inf, PointXYZ& p)
{
	inf.read((char *)&(p.x),sizeof(double));
	inf.read((char *)&(p.y),sizeof(double));
	inf.read((char *)&(p.z),sizeof(double));
	return inf;
}

ofstream& operator<<(ofstream& outf, const PointXYZ& p)
{
	outf << p.x << " " << p.y << " " << p.z;
	return outf;
}

ofstream& operator<(ofstream& outf, const PointXYZ& p)
{
	outf.write((char *)&(p.x),sizeof(double));
	outf.write((char *)&(p.y),sizeof(double));
	outf.write((char *)&(p.z),sizeof(double));
	return outf;
}
