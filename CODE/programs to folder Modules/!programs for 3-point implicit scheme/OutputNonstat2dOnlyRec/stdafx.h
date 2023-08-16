#pragma once

#include <io.h>
#include <fcntl.h>
#include <share.h>
#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <direct.h>
#include <math.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include <omp.h>

#define RETCODE_OK				0x0000
#define RETCODE_NOMEM			0x0001
#define RETCODE_NOFILE			0x0002
#define RETCODE_OUTOFRANGE		0x0004
#define RETCODE_SQFROMNEG		0x0008
#define RETCODE_DEVBYZERO		0x0010
#define RETCODE_NOTINIT			0x0020
#define RETCODE_BADFILE			0x0040
#define RETCODE_ERROR			0x0080
#define RETCODE_NOANOMALOBJECTS	0x0100

#define PI 3.1415926535897932
#define MU0 1.25663706143592e-6

const double d_eps=1e-9;
const double d_eps_1=1e-14;
const double d_eps_2=1e-6;
const double d_eps_3=1e-3;
const double line_eps=1e-3;

extern int NUMBEROFABDIPOLES;

using namespace std;

__forceinline ostream& operator < (ostream& file,const double& data)
{
	file.write((char*)&data,sizeof(data));
	return file;
}
__forceinline ostream& operator < (ostream& file,const int& data)
{
	file.write((char*)&data,sizeof(data));
	return file;
}

__forceinline istream& operator > (istream& file,double&  data)
{
	file.read((char*)&data,sizeof(data));
	return file;
}

__forceinline istream& operator > (istream& file,int&  data)
{
	file.read((char*)&data,sizeof(data));
	return file;
}

struct PointXYZ
{
	double x,y,z;
	PointXYZ(){x=y=z=0.0;}
	void AddPoint(PointXYZ &p,double cff)
	{
		x+=cff*p.x;
		y+=cff*p.y;
		z+=cff*p.z;
	}
	void Clear(){x=y=z=0.0;}
};

struct PointRZ
{
	double r,z;
};

struct Rect
{
	int nodes[5];	//!<  
	int rtype;		//!<  
	int mtr;		//!<  
};

int sig(const double &v);
int FindInDoubleVec(vector<double> &vec,double elem);
int FindIntervalInDoubleMas(double *vec,int size,double elem);
int FindValueInDoubleMas(double *vec,int size,double elem);
int FindIntervalInDoubleMas(double *vec, int beg, int end, double elem);
