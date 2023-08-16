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
 *  This file contains the headers of auxiliary utilities for outputting with smoothing
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin                         
 *  Novosibirsk State Technical University,                 
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia       
 *  vdv_wk@mail.ru                                          
 *  Version 2.0 January 16, 2023                            
*/

#pragma once
#include "T_Mapping.h"


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>


using namespace std;

#define OUTPUT
class PointRes2
{
public:
	int num; //  
	double point[3]; //  

	PointRes2()	{}
	~PointRes2() {}
};
struct PointRes_less_x : public binary_function<PointRes2, PointRes2, bool>
{
	bool operator() (const PointRes2& a, const PointRes2& b) const {return a.point[0] < b.point[0];}
};
struct PointRes_less_y : public binary_function<PointRes2, PointRes2, bool>
{
	bool operator() (const PointRes2& a, const PointRes2& b) const {return a.point[1] < b.point[1];}
};
struct PointRes_less_z : public binary_function<PointRes2, PointRes2, bool>
{
	bool operator() (const PointRes2& a, const PointRes2& b) const {return a.point[2] < b.point[2];}
};
class Long_double_eq : public unary_function<long_double, bool>
{
	long_double a;
public:
	explicit Long_double_eq(const long_double& aa) : a(aa) {}
	bool operator() (const long_double& b) const {return fabs(b.d - a.d) < 1e-6;}
};
class Long_double_num_eq : public unary_function<long_double, bool>
{
	long_double a;
public:
	explicit Long_double_num_eq(const long_double& aa) : a(aa) {}
	bool operator() (const long_double& b) const {return a.i == b.i;}
};
struct Long_double_less : public binary_function<long_double, long_double, bool>
{
	bool operator() (const long_double& a, const long_double& b) const {return a.d < b.d;}
};
struct Long_double_num_less : public binary_function<long_double, long_double, bool>
{
	bool operator() (const long_double& a, const long_double& b) const {return a.i < b.i;}
};

class Output3dArbitrary
{
protected:
	int withSpline3d;
	int withSpline2d;
	int zeroPlane; //    z=0       -
	int kpar;
	int kuzlov;
	int n_pointres;
	double (*pointres)[3];
	double (*xyz)[3];

	int (*nver)[14];
	int (*nvtr)[8];

	vector<int> elemForPoint; //     -,   
	vector< vector<int> > PointresForElem; //   -   , -   
	vector<long_double> PointresXsorted;
	vector<long_double> PointresYsorted;
	vector<long_double> PointresZsorted;
	
public:
	Output3dArbitrary(int withSpline3d, int withSpline2d, int zeroPlane, int kuzlov, int kpar, int n_pointres,
		double (*pointres)[3], double (*xyz)[3]);
	~Output3dArbitrary();

	void PointresSort();
	int FindElemForReceivers();

	int GetGlobalVertex(int nElem, int nLocVertex);
};

class OutputNode3d : public Output3dArbitrary
{
public:
	OutputNode3d(int withSpline3d, int withSpline2d, int zeroPlane,
		int kuzlov, int kpar, int n_pointres,
		double (*pointres)[3], double (*xyz)[3], int (*nvtr)[8]);

	~OutputNode3d();

	int OutputFieldAtReceiversHarm(double *v3, vector<vector<double>> &res, double w, int is_b);
	int OutputFieldAtReceivers(double *v3, vector<double> &result, int derive);
};

