#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include<iostream>
#include<fstream>
#include<iomanip>

#include<vector>

/*!     FEM */
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

using namespace std;

__forceinline ostream& operator < (ostream& file,const double& data)
{
	file.write((char*)&data,sizeof(data));
	return file;
}

__forceinline istream& operator > (istream& file,double&  data)
{
	file.read((char*)&data,sizeof(data));
	return file;
}

__forceinline ostream& operator < (ostream& file,const float& data)
{
	file.write((char*)&data,sizeof(data));
	return file;
}

__forceinline istream& operator > (istream& file,float&  data)
{
	file.read((char*)&data,sizeof(data));
	return file;
}

__forceinline ostream& operator < (ostream& file,const int& data)
{
	file.write((char*)&data,sizeof(data));
	return file;
}

__forceinline istream& operator > (istream& file,int&  data)
{
	file.read((char*)&data,sizeof(data));
	return file;
}

int sig(const double &v);

int FindInDoubleVec(vector<double> &vec,double elem);
int FindIntervalInDoubleMas(double *vec,int size,double elem);
int FindIntervalInDoubleMas(double *vec,int beg, int end, double elem);

#define PI 3.1415926535897932

struct PointXYZ
{
	double x,y,z;
};

extern int NUMBEROFABDIPOLES;

const int size_i=sizeof(int);
