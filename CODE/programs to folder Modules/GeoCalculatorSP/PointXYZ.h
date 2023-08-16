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
 *  This file contains the headers for work with 3D points
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin                                
 *  Novosibirsk State Technical University,                        
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia              
 *  Corresponding author: vdv_wk@mail.ru                           
 *  Version 2.0 January 16, 2023                                   
*/

#pragma once

// The class contains the node data of a 3D grid
struct PointXYZ
{
	double x, y, z; //!<  
	PointXYZ();
	PointXYZ(const double& _x, const double& _y, const double& _z);
	const PointXYZ& operator+=(const PointXYZ& p);
	const PointXYZ& operator-=(const PointXYZ& p);
};

PointXYZ operator*(const PointXYZ& p, const double& a);
PointXYZ operator+(const PointXYZ& p1, const PointXYZ& p2);
ifstream& operator>>(ifstream& inf, PointXYZ& p);
ifstream& operator>(ifstream& inf, PointXYZ& p);
ofstream& operator<<(ofstream& outf, const PointXYZ& p);
ofstream& operator<(ofstream& outf, const PointXYZ& p);
