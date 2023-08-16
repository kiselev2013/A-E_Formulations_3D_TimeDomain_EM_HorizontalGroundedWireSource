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
 *  This file contains the headers for working with a 2D mesh and solution
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 January 16, 2023
*/

#pragma once

struct PointRZ
{
	double r, z;
};

struct Rect
{
	int nodes[5];
	int rtype;
	int mtr;
};

bool inTriangle(PointRZ *trg,const PointRZ &ptmp);
bool CalcBaseFunctions(PointRZ *trg,double a[3][3]);
double scal(const PointRZ &p1,const PointRZ &p2);
double dist(const PointRZ &p1,const PointRZ &p2);

struct task2d
{
	char path[256];

	int kpnt, krect, nc;
	PointRZ* pnt;
	Rect* rect;

	double *v2s,*v2c;

	int nmat3d;
	int *mtr3d2d;

	int nmat2d;
	double *sigma;

	int nreg, qr, qz;
	int *reg;
	double *rm, *zm;

	int *indr, *indz;

	task2d(char *_path);
	~task2d();

	int Read();
	int ReadRZind();
	int ReadSolution(int npls);

	double rmin,rmax,zmin,zmax;

	void GenKU(int nGels,int nVels,vector<double> &zgel,vector<double> &zvela,vector<double> &zvelb);
};
