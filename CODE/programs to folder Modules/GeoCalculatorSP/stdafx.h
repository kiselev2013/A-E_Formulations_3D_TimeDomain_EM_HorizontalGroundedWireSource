
#pragma once

#include <omp.h>

#include <vector>
#include <queue>
#include <list>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <stack>
#include <algorithm>
#include <cmath>
#include <memory>
#include <sstream>
#include <iomanip>
#include <sstream>
#include <numeric>
#include <stdio.h>
#include <conio.h>
#include <math.h>
#include <direct.h>
#include <stdlib.h>
#include <string.h>
#include <windows.h>
#include <time.h>

using namespace std;

#include <sys/stat.h>

const int size_i=sizeof(int);
const int size_d=sizeof(double);

extern stringstream stringbuffer;

#define EpsRect3D 1e-6
/*!    */
#define MAXNUMBEROFTIMES 10000
/*!    */
#define MAXNUMBEROFMATERIALS 10000

#define _PI_ 3.14159265358979323846
#define MU0 4e-7*_PI_

#define strcpy__m(str_d, str_s) str_d=str_s
#define strcat__m(str_d, str_s) (str_d+str_s).c_str()

#define CHECKSTRING0(s, n) if (int(sizeof(s)/sizeof(s[0])-1-n)<0) CloseProgramm(1);
#define CHECKSTRING0WITHEXCEPTION(s, n) if (int(sizeof(s)/sizeof(s[0])-1-n)<0) throw logic_error("too int string");

#define DIAGNULL 1e-60
#define DENOMNULL 1e-60
