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


const double PI=3.1415926535897932;

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
