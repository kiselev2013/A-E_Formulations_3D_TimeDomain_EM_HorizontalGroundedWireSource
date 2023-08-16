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
#include <vector>
#include "ListV.h"
#include "mesh3d.h"

#define MAXNUMBEROFTIMES 10000

#define MAXNUMBEROFMATERIALS 10000

// The class contains the node data of a 2D grid
struct PointRZ
{
	double r, z; //!<  
	PointRZ() 
	{
	}
	PointRZ(const double& _r, const double& _z) 
	{ 
		r=_r; 
		z=_z; 
	}
	PointRZ operator+(const PointRZ& prz) const
	{
		return PointRZ(r+prz.r, z+prz.z);
	}
};

ifstream& operator>(ifstream& inf, PointRZ& p);

#ifdef USEREGINFO

/*!         .   */
struct ResCalcElInfo
{
	bool isResCalc;
	ListOfValues<int> ResCalcPoints;
	ResCalcElInfo()	{ isResCalc=false; }

};

/*!    1d  */
int GetRegElem1DNum(const int& _qc0, const int& _qc, const double* _cm, const double& c, const double& eps);

/*!    2d  */
int GetRegElem2DNum(const int& _qr, const double* _rm,
					const int& _qz, const double* _zm,
					const PointRZ& p, const double& eps);

#endif

// The class contains the data of a 2D grid element
struct Rect
{
	int nodes[5];	//!<  
	int rtype;		//!<  
	int mtr;		//!<  

	/*!     */
	bool in(const PointRZ &p, const PointRZ *pntArray) const
	{
		const PointRZ& p0=pntArray[nodes[0]-1];
		const PointRZ& p3=pntArray[nodes[3]-1];
		return (p0.r<=p.r&&p.r<=p3.r&&
				p0.z<=p.z&&p.z<=p3.z);
	}
	/*!    , +-eps */
	bool in(const PointRZ &p, const PointRZ *pntArray, const double& eps) const
	{
		const PointRZ& p0=pntArray[nodes[0]-1];
		const PointRZ& p3=pntArray[nodes[3]-1];
		return ( (p0.r-eps<=p.r)&&(p.r<=p3.r+eps) &&
				 (p0.z-eps<=p.z)&&(p.z<=p3.z+eps) );
	}
	/*!        */
	void calcH(const PointRZ &p, const PointRZ *pntArray, double& ra, double& rb, double& za, double& zb) const
	{
		const PointRZ& p0=pntArray[nodes[0]-1];
		const PointRZ& p3=pntArray[nodes[3]-1];
		ra=p.r-p0.r,
		rb=p3.r-p.r,
		za=p.z-p0.z,
		zb=p3.z-p.z;
	}
};

ifstream& operator>(ifstream& inf, Rect& r);

/*!    */
enum MeshOptions
{
	moForSpider			=0x0001,	// Afi
	moForDipole			=0x0002,	// Hfi
	moUnloadSolutions	=0x0004,
	moReadWkInfo		=0x0008,
	moDontReadResPoints =0x0010,
	moForLine3d			=0x0020,
	moUsePntMtrInfo		=0x0040,
	moForPolygon		=0x0080,
	moCalcConfResB		=0x0100,
	moForDisk			=0x0200,
	moForVEL			=0x0400,
	moForMTZ			=0x0800,
	moForHarmLoop		=0x1000,
	moForCED			=0x2000,
	moForPED			=0x4000
};

/*!     */
struct PointForResultCalculationInMesh3D
{
	PointXYZ p;
	PointXYZ* res;
	double xa, xb, ya, yb, za, zb;
	int rnum; // from 0 !
	bool isfound;
	PointForResultCalculationInMesh3D()
	{
		isfound=false;
		res=NULL;
	}
	~PointForResultCalculationInMesh3D()
	{
		if (res) delete [] res;
	}
};

/*!     */
struct PointForResultCalculation
{
	PointXYZ p;					//!<   xyz
	PointXYZ* res;				//!<   
	double r, ra, rb, za, zb;	//!< r-        
	double sinfi, cosfi;		//!< sin, cos
	int rnum;					//!<  ,  0
	bool isfound;				//!<    
	PointForResultCalculation()
	{
		isfound=false;
		res=NULL;
	}
	~PointForResultCalculation()
	{
		if (res) delete [] res;
	}
};

/*!     */
struct RectForResultCalculationInMesh3D
{
	int rnum;		// number of rect from 0 !
	double qx[8];	// weights of basic x funcs
	double qy[8];	// weights of basic y funcs
	double qz[8];	// weights of basic z funcs
	double qv[8];	// weights of basic v funcs
	double hx, hy, hz;	
};


/*!     */
struct RectForResultCalculation
{
	int rnum;		//!<  ,  0
	double q[4];	//!<   
	double q_dr[4];	//!<        r
	double q_dz[4];	//!<        z
	double hr, hz;	//!<  
};

/*!        */
struct ResCalcPointBounds
{
	int beginP, endP;
};

// The class contains a grid for 2D scalar problems
struct MeshForVP
{
	const char* taskpath;	//!<     

	int meshoptions;		//!<  	

	PointXYZ A, B;			//!<   
	vector<PointXYZ> vB;	//!<     
	int NvB;				//!<     

	double cedr;

	int kpnt, krect, kt1, nc;	//!< - , ,   1,  
	PointRZ* pnt;				//!<  
	Rect* rect;					//!<  
	int* l1;					//!<    1

	int c1, c2, c1b, c2b;		//!<     2

	double Rbak;				//!< 		

	int ntimes;						//!< - 
	double times[MAXNUMBEROFTIMES];	//!< 

	double sigma[MAXNUMBEROFMATERIALS];	//!<  
	double currc;	//!<  

	int qMN;		//!< -  

	int ResCalcPointNum;							//!< -    
	PointForResultCalculation* ResCalcPointArray;	//!<     
	int ResCalcRectNum;								//!< -    
	RectForResultCalculation* ResCalcRectArray;		//!<     
#ifndef USEREGINFO
	PointRZ gResCalcPoint1, gResCalcPoint2;			//!<   ,    
#endif
	ResCalcPointBounds boundsVA, boundsVB;						//!<     
	ResCalcPointBounds boundsVforMNfromA, boundsVforMNfromB;	//!<     
	ResCalcPointBounds boundsVforEfromA, boundsVforEfromB;		//!<     

	double *qv;		//!< 
	double *q_all;

	int nreg, qr, qz;
	int *reg;
	double *rm, *zm;

	MeshForVP(int mOptions);
	
	~MeshForVP();

	int Read(const char* path, bool read_times=false);
	int ReadRecievers(const char* path, int pntVsize, const PointXYZ* pntV);

	/*!      */
	int GetWeights(const double* q, bool storeq=false);

	/*!     */
	int ReadMtrCatalog(const char* fname, double* cat);

	/*!   xyzMN */
	int ReadResMNCalcFile(const char* fname, ListOfValues<PointForResultCalculation>& _listP, ResCalcPointBounds& _boundsP,
#ifndef USEREGINFO
		PointXYZ& gp1, PointXYZ& gp2, // gabarit points
		const bool& first, 
#endif
		const PointXYZ& pEL);

	/*!   xyzVector* */
	int ReadResCalcFile(const char* fname, ListOfValues<PointForResultCalculation>& _listP, ResCalcPointBounds& _boundsP,
#ifndef USEREGINFO
		PointXYZ& gp1, PointXYZ& gp2, // gabarit points
		const bool& first, 
#endif
		const PointXYZ& pEL);

	/*!    */
	int ReadTimesFile(const char* fname);
};
