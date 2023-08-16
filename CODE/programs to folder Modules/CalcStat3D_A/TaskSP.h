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
 *  This file contains the headers for calculation of secondary part of field A for stationary task
 *  
 *  
 *  Written by Prof. Yuri G. Soloveichik, Prof. Marina G. Persova, Ph.D. Denis V. Vagin 
 *  Novosibirsk State Technical University,                                             
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                   
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                    
 *  Version 2                                                                           
*/

#pragma once
#include "PointVector.h"
#include "ListV.h"
#include "ArrayOf.h"
#include "gauss3.h"

extern bool CheckStop(void);

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

#define UngroundedLineAsPoints

#define CHECKSTRING0WITHEXCEPTION(s, n) if (int(sizeof(s)/sizeof(s[0])-1-n)<0) throw logic_error("too int string");

#ifndef UngroundedLineAsPoints 
#define NUMBER_OF_DERIVATIVES_FOR_UNGROUNDED_LINE 4
#else
#define NUMBER_OF_DERIVATIVES_FOR_UNGROUNDED_LINE 8
#endif

#define EpsRect3D 1e-6
/*!    */
#define MAXNUMBEROFTIMES 10000
/*!    */
#define MAXNUMBEROFMATERIALS 10000

#define _PI_ 3.14159265358979323846
#define MU0 4e-7*_PI_

#define USEREGINFO

#define strcpy__m(str_d, str_s) str_d=str_s
#define strcat__m(str_d, str_s) (str_d+str_s).c_str()

#define CHECKSTRING0(s, n) if (int(sizeof(s)/sizeof(s[0])-1-n)<0) return RETCODE_OUTOFRANGE;

extern ofstream logfile;

void mult_matrix(double a[][3],double b[][3],double c[][3]);

double GetDeterminant33(double m[][3]);
void TransposeMatrix33(double a[][3]);
void InverseMatrix33(double m[][3],double m_1[][3]);

void loc_correct(double &ksi,double &eta,double &dzeta);

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

// The class contains the node data of a 3D grid
struct PointXYZ
{
	double x, y, z; //!<  
	PointXYZ() 
	{ 
		x=y=z=0; 
	}
	PointXYZ(const double& _x, const double& _y, const double& _z) 
	{
		x=_x;
		y=_y;
		z=_z; 
	}
	const PointXYZ& operator+=(const PointXYZ& p)
	{
		x+=p.x;
		y+=p.y;
		z+=p.z;
		return *this;
	}
	const PointXYZ& operator-=(const PointXYZ& p)
	{
		x-=p.x;
		y-=p.y;
		z-=p.z;
		return *this;
	}
};

PointXYZ operator*(const PointXYZ& p, const double& a);
PointXYZ operator+(const PointXYZ& p1, const PointXYZ& p2);
ifstream& operator>>(ifstream& inf, PointXYZ& p);
ifstream& operator>(ifstream& inf, PointXYZ& p);
ofstream& operator<<(ofstream& outf, const PointXYZ& p);
ofstream& operator<(ofstream& outf, const PointXYZ& p);

// The class contains the data of a 3D grid element
struct Rect3D
{
	int nodes[13];
	int rtype;
	int mtr;
	bool in(const PointXYZ &p, const PointXYZ *pntArray) const
	{
		const PointXYZ& p0=pntArray[nodes[0]-1];
		const PointXYZ& p7=pntArray[nodes[7]-1];
		return ((p0.x-EpsRect3D)<=p.x&&p.x<=(p7.x+EpsRect3D)&&
				(p0.y-EpsRect3D)<=p.y&&p.y<=(p7.y+EpsRect3D)&&
				(p0.z-EpsRect3D)<=p.z&&p.z<=(p7.z+EpsRect3D));
	}
	void calcH(const PointXYZ &p, const PointXYZ *pntArray, 
				double& xa, double& xb, double& ya, double& yb, double& za, double& zb) const
	{
		const PointXYZ& p0=pntArray[nodes[0]-1];
		const PointXYZ& p7=pntArray[nodes[7]-1];
		xa=p.x-p0.x,
		xb=p7.x-p.x,
		ya=p.y-p0.y,
		yb=p7.y-p.y,
		za=p.z-p0.z,
		zb=p7.z-p.z;
	}
	PointXYZ GetMassCenter(const PointXYZ *pntArray) const
	{
		PointXYZ mc(0, 0, 0);
		mc+=pntArray[nodes[0]-1];
		mc+=pntArray[nodes[7]-1];
		mc.x*=0.5; 
		mc.y*=0.5; 
		mc.z*=0.5;
		return mc;
	}
	double GetVolume(const PointXYZ *pntArray) const
	{
		const PointXYZ& p0=pntArray[nodes[0]-1];
		const PointXYZ& p7=pntArray[nodes[7]-1];
		return (p7.x-p0.x)*(p7.y-p0.y)*(p7.z-p0.z);
	}
};

ifstream& operator>(ifstream& inf, Rect3D& r);

#ifdef USEREGINFO

/*!         .   */
struct ResCalcElInfo
{
	bool isResCalc;
	ListOfValues<int> ResCalcPoints;
	ResCalcElInfo()	{ isResCalc=false; }

};

#endif

void GetNodeSolutionByLocalcCoords(double *lc,double *q,double &Val);
void GetNodeCoefficientsByLocalcCoords(double *lc,double *cff);
void GetGlobalCoordinates(PointXYZ *HexPnt,double *lc,double *gc);
void GetSolutionOnHex(PointXYZ &R,PointXYZ *CntHex,double *q,double &Val);
void FindLocalCoordinates(PointXYZ &R,PointXYZ *HexPnt,double *lc);
bool CheckPointInHex(PointXYZ &R,PointXYZ *HexPnt,double &mdist);

/*!    */
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

/*! ,   - */
struct DeltaFunction
{
	/*!     */
	struct DFValue
	{
		double t, v; //!<    
	};

	PointRZ dfPoint;	//!< 
	int dfNode;			//!<     ,  0
	DFValue* dfArray;	//!<  
	int dfArraySize;	//!< - 

	DeltaFunction()
	{
		dfArray=NULL;
	}
	~DeltaFunction()
	{
		if (dfArray) delete [] dfArray;
	}

	int Read(const char* path)
	{
		int i, sz;
		char buf[1024];
		string fname=path;
		ifstream inf((fname+"\\DeltaFunction").c_str());
		if (!inf) return RETCODE_NOFILE;
		inf>>dfPoint.r/*>>dfPoint.z*/;
		inf>>sz;
		dfArray=new DFValue[dfArraySize=sz];
		if (dfArray==NULL) return RETCODE_NOMEM;
		for (i=0; i<sz; i++)
		{
			DFValue& dfValuei=dfArray[i];
			if (!inf.good()) return RETCODE_BADFILE;
			inf>>dfValuei.t>>dfValuei.v;
		}
		inf.close();
		inf.clear();
		if (strcmp(fname.c_str(), ".")!=0)
			fname=".";
		inf.open((fname+"\\geoprep.dat").c_str());
		if (!inf) return RETCODE_NOFILE;
		inf.getline(buf, 1000, '\n');
		inf.getline(buf, 1000, '\n');
		inf.getline(buf, 1000, '\n');
		inf.getline(buf, 1000, '\n');
		inf.getline(buf, 1000, '\n');
		inf.getline(buf, 1000, '\n');
		inf.getline(buf, 1000, '\n');
		inf>>dfPoint.z;
		inf.close();
		inf.clear();
		return RETCODE_OK;
	}
};

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
	moForCED			=0x2000
};

/*!     */
struct PointForResultCalculationInMesh3D
{
	PointXYZ p;
	PointXYZ* res;
	double xa, xb, ya, yb, za, zb, ksi, eta, zet;
	int rnum; // from 0 !
	bool isfound,ishex;
	PointForResultCalculationInMesh3D()
	{
		ishex=false;
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

bool CompareFiles(const char* _f1, const char* _f2);

#define DIAGNULL 1e-60
#define DENOMNULL 1e-60

int MatrixOnVector(
	const double *x1,	// z=A(x1+x2)// IN 
	const double *x2,	// 			 //
	const double *adiag,// matrix	 //
	const double *altr,	//			 //
	const double *autr,	//			 //
	const int *iptr,	// portrait	 //
	const int *jptr,	//			 //
	const int &n,		// size		 //
	double *z);			//			 // OUT

int MatrixOnVector(
	const double *x,	// z=Ax		// IN 
	const double *adiag,// matrix	//
	const double *altr,	//			//
	const double *autr,	//			//
	const int *iptr,	// portrait	//
	const int *jptr,	//			//
	const int &n,		// size		//
	double *z);			//			// OUT	

void CopyV(const int &n, const double *from, double *to);
double Scal(const int &n, const double *v1, const double *v2);

// The class contains a finite element SLAE
struct SparseSLAE
{
	int n;
	int jsize;
	int* iptr;
	int* jptr;

	double *ldiag, *ltr, *utr;
	double *xk, *rk, *zk, *pk, *vk, *Urk, *LAUrk;
	double *z_omp;
	
	struct SparseMatrix
	{
		double* adiag;
		double* altr;
		double* autr;
		SparseMatrix()
		{
			adiag=NULL;
			altr=NULL;
			autr=NULL;
		}
		~SparseMatrix()
		{
			if (adiag) delete [] adiag;
			if (altr) delete [] altr;
			if (autr) delete [] autr;
		}
		int Init(const int& size, const int& trsize)
		{
			if ((adiag=new double[size])==NULL) return RETCODE_NOMEM;
			if ((altr=new double[trsize])==NULL) return RETCODE_NOMEM;
			if ((autr=new double[trsize])==NULL) return RETCODE_NOMEM;
			Clear(size, trsize);
			return RETCODE_OK;
		}
		void Clear(const int& size, const int& trsize)
		{
			int i;
			for (i=0; i<size; i++) 
				adiag[i]=0;
			for (i=0; i<trsize; i++) 
				altr[i]=autr[i]=0;
		}
		void AddToMatrix(const int& size, const int& trsize, const SparseMatrix& Matrix, const double& coef=1)
		{
			int i;
			for (i=0; i<size; i++) 
				adiag[i]+=coef*Matrix.adiag[i];
			for (i=0; i<trsize; i++)
			{
				altr[i]+=coef*Matrix.altr[i];
				autr[i]+=coef*Matrix.autr[i];
			}
		}
	};

	struct VectorForSparseSLAE
	{
		double *v;
		VectorForSparseSLAE()
		{
			v=NULL;
		}
		~VectorForSparseSLAE()
		{
			if (v) delete [] v;
		}
		int Init(const int& size)
		{
			if ((v=new double[size])==NULL) return RETCODE_NOMEM;
			Clear(size);
			return RETCODE_OK;
		}
		void Clear(const int& size)
		{
			for (int i=0; i<size; i++) 
				v[i]=0;
		}
		void AddToVector(const int& size, const VectorForSparseSLAE& Vec, const double& coef=1)
		{
			for (int i=0; i<size; i++) 
				v[i]+=coef*Vec.v[i];
		}
	};

	SparseMatrix A, B, C, C0, Dr, Dz;
	VectorForSparseSLAE F, G, tmpV; 
	VectorForSparseSLAE U, F3; // U(j), U(j-1), U(j-2)
	VectorForSparseSLAE *pU; // pointers to U(j), U(j-1), U(j-2)

	SparseSLAE()
	{
		n=jsize=0;
		iptr=NULL;
		jptr=NULL;
		pU=&U;

		ldiag=NULL; 
		ltr=NULL;
		utr=NULL;
		xk=NULL;
		rk=NULL;
		zk=NULL;
		pk=NULL;
		vk=NULL; 
		Urk=NULL; 
		LAUrk=NULL;
		z_omp=NULL;
	}

	~SparseSLAE()
	{
		if (iptr) delete [] iptr;
		if (jptr) delete [] jptr;
		
		if (ldiag) delete [] ldiag;
		if (ltr) delete [] ltr;
		if (utr) delete [] utr;
		if (xk) delete [] xk;
		if (rk) delete [] rk;
		if (zk) delete [] zk;
		if (pk) delete [] pk;
		if (vk) delete [] vk;
		if (Urk) delete [] Urk;
		if (LAUrk) delete [] LAUrk;
		if (z_omp) {delete [] z_omp; z_omp=NULL;}
	}

	int AddToMatrix(const int &i, const int &j, const double &val, SparseMatrix& Matrix)
	{
		if (i==j) 
		{ 
			Matrix.adiag[i]+=val; 
			return RETCODE_OK;
		}
		if (i>j) 
		{
			for(int k=iptr[i]; k<iptr[i+1]; k++) 
				if (jptr[k-1]==j+1) 
				{ 
					Matrix.altr[k-1]+=val; 
					return RETCODE_OK;
				}
		}
		else
		{
			for(int k=iptr[j]; k<iptr[j+1]; k++) 
				if (jptr[k-1]==i+1) 
				{ 
					Matrix.autr[k-1]+=val; 
					return RETCODE_OK;
				}
		}
		return RETCODE_OUTOFRANGE;
	}

	void SetMatrixDiagElement(const int &i, const double &val, SparseMatrix& Matrix)
	{
		Matrix.adiag[i]=val;
	}

	bool ZeroMatrixDiagElement(const int &i, SparseMatrix& Matrix)
	{
		return fabs(Matrix.adiag[i])<1e-6;
	}

	void AddToVector(const int &i, const double &val, VectorForSparseSLAE& Vec) 
	{ 
		Vec.v[i]+=val;
	}

	void SetVectorElement(const int &i, const double &val, VectorForSparseSLAE& Vec)
	{
		Vec.v[i]=val;
	}

	void CorrectDiagonal()
	{
		for (int i=0; i<n; i++)
			if (fabs(A.adiag[i])<DIAGNULL) { A.adiag[i]=1e30; F.v[i]=0; }
	}
	void CorrectOnlyDiagonal()
	{
		for (int i=0; i<n; i++)
			if (fabs(A.adiag[i])<DIAGNULL) { A.adiag[i]=1; }
	}

	int SimmetrifyA()
	{
		for (int i=0; i<jsize; i++)
		{
			double &altri=A.altr[i];
			double &autri=A.autr[i];
			if (fabs(altri)<1e-13)
				altri=0;
			if (fabs(autri)<1e-13)
				autri=0;
		}
		return RETCODE_OK;
	}

	void Nulify(const int& k)
	{
		int i;
		for (i=iptr[k-1]-1; i<iptr[k]-1; i++)
			A.altr[i]=0.;
		for (i=0; i<jsize; i++)
			if (jptr[i]==k)
				A.autr[i]=0.;
	}

};

struct LocalMatrix3D
{
	double m[8][8];
	LocalMatrix3D()
	{
		int i,j;
		for(i=0;i<8;i++)
		{
			for(j=0;j<8;j++)
			{
				m[i][j]=0.0;
			}
		}
	}
};

// The class contains a grid for 3D scalar problems
class Mesh3DForVP
{
public:

	const char* taskpath;	//!<     

	int meshoptions;	//!<  

	int kpnt, krect3d, kt1, nc;	//!< - , ,   1,  
	PointXYZ* pnt;				//!<  
	Rect3D* rect3d;				//!<  
	int* l1;					//!<    1

	int ntimes;						//!< - 
	double times[MAXNUMBEROFTIMES];	//!< 

	double sig[MAXNUMBEROFMATERIALS];	//!<   3d
	double sig0[MAXNUMBEROFMATERIALS];	//!<   2d

	int ResCalcPointNum;									//!< -    
	PointForResultCalculationInMesh3D* ResCalcPointArray;	//!<     
	int ResCalcRectNum;										//!< -    
	RectForResultCalculationInMesh3D* ResCalcRectArray;		//!<     
	PointXYZ gResCalcPoint1, gResCalcPoint2;				//!<   ,    
	
	int qMN;												//!< -  
	ResCalcPointBounds boundsV;								//!<     ,   MN
	ResCalcPointBounds boundsVforE;							//!<     
	ResCalcPointBounds boundsB;								//!<     
	ResCalcPointBounds boundsA;								//!<     

	double *q_all;		//!< 
	double *Ax;
	double *Ay;
	double *Az;

	double *qva;
	double *qvb;
	double *qv;
	double *qv_sum;

	double *q_res1;											//!<   
	double *q_res2;											//!<   
	
	vector<int> regular;									//!<      

	int n_mesh_regular_x;									
	int n_mesh_regular_y;
	int n_mesh_regular_z;
	double *mesh_regular_x;
	double *mesh_regular_y;
	double *mesh_regular_z;

	vector<double> status;

	int npls;

	Mesh3DForVP(int mOptions=0):
		taskpath("."),
		meshoptions(mOptions)
	{
		pnt=NULL;
		q_all=Ax=Ay=Az=NULL;
		qva=NULL; qvb=NULL; qv=NULL; qv_sum=NULL;
		rect3d=NULL;
		l1=NULL;
		ResCalcPointArray=NULL;
		ResCalcRectArray=NULL;
		ResCalcPointNum=0;
		ResCalcRectNum=0;

		mesh_regular_x=NULL;
		mesh_regular_y=NULL;
		mesh_regular_z=NULL;

		q_res1=NULL;
		q_res2=NULL;

		nc=0;
	}
	~Mesh3DForVP()
	{
		if (pnt) delete [] pnt;
		if (q_all) delete [] q_all;
		if (Ax) delete [] Ax;
		if (Ay) delete [] Ay;
		if (Az) delete [] Az;
		if (qva) delete [] qva;
		if (qvb) delete [] qvb;
		if (qv) delete [] qv;
		if (qv_sum) delete [] qv_sum;
		if (rect3d) delete [] rect3d;
		if (l1) delete [] l1;
		if (ResCalcPointArray) delete [] ResCalcPointArray;
		if (ResCalcRectArray) delete [] ResCalcRectArray;

		if (mesh_regular_x) delete [] mesh_regular_x;
		if (mesh_regular_y) delete [] mesh_regular_y;
		if (mesh_regular_z) delete [] mesh_regular_z;
	}

	int GetNearestElement(double *gc,double *lc,int fuseair);

	int Read(const char* path, bool read_times=false)
	{
		int i, j, retc, qrptotal=0;
		string buf;
		ifstream inf, rzf, nvtrf, nvkatf, l1f, wkf;
		int fuseair;

		taskpath=path;

		strcpy__m(buf, path);
		inf.open(strcat__m(buf, "\\inftry.dat"));
		if (!inf) return RETCODE_NOFILE;
		inf.ignore(1000, '\n');
		inf.ignore(1000, '='); inf>>kpnt;
		inf.ignore(1000, '='); inf>>krect3d;
		inf.ignore(1000, '='); inf>>kt1;
		inf.close();
		inf.clear();

		strcpy__m(buf, path);
		rzf.open(strcat__m(buf, "\\xyz.dat"), ios::binary);
		if (!rzf) return RETCODE_NOFILE;
		pnt=new PointXYZ[kpnt];
		if (!pnt) return RETCODE_NOMEM;
		for (i=0; i<kpnt; i++)
			rzf > pnt[i];
		rzf.close();

		q_all=new double[kpnt];
		if (q_all==NULL) return RETCODE_NOMEM;
		Ax=new double[kpnt];
		if (Ax==NULL) return RETCODE_NOMEM;
		Ay=new double[kpnt];
		if (Ay==NULL) return RETCODE_NOMEM;
		Az=new double[kpnt];
		if (Az==NULL) return RETCODE_NOMEM;

		strcpy__m(buf, path);
		nvtrf.open(strcat__m(buf, "\\nver.dat"), ios::binary);
		if (!nvtrf) return RETCODE_NOFILE;
		rect3d=new Rect3D[krect3d];
		if (!rect3d) return RETCODE_NOMEM;
		for (i=0; i<krect3d; i++)
			nvtrf > rect3d[i];
		nvtrf.close();

		strcpy__m(buf, path);
		nvkatf.open(strcat__m(buf, "\\nvkat.dat"), ios::binary);
		if (!nvkatf) return RETCODE_NOFILE;
		for (i=0; i<krect3d; i++)
			nvkatf > rect3d[i].mtr;
		nvkatf.close();

		strcpy__m(buf, path);
		l1f.open(strcat__m(buf, "\\l13d.dat"), ios::binary);
		if (!l1f) return RETCODE_NOFILE;
		l1=new int[kt1];
		if (!l1) return RETCODE_NOMEM;
		for (i=0; i<kt1; i++)
			l1f > l1[i];
		l1f.close();

		nc=kpnt;

		if ((qva=new double[kpnt/*nc*/*npls])==NULL) return RETCODE_NOMEM;
		if ((qvb=new double[kpnt/*nc*/*npls])==NULL) return RETCODE_NOMEM;
		if ((qv=new double[kpnt/*nc*/])==NULL) return RETCODE_NOMEM;
		if ((qv_sum=new double[kpnt/*nc*/])==NULL) return RETCODE_NOMEM;


		strcpy__m(buf, path);
		if ((retc=ReadMtrCatalog(strcat__m(buf, "\\sig3d"), sig, sig0))!=0)
			return retc;


		qMN=0;

		ListOfValues<PointForResultCalculationInMesh3D> PointList;

		strcpy__m(buf, path);
		if ((retc=ReadResMNCalcFile(strcat__m(buf, "\\xyzMN"), 
								PointList, 
								boundsV, 
								gResCalcPoint1, 
								gResCalcPoint2, 
								true))!=0)
			return retc;

		strcpy__m(buf, path);
		if ((retc=ReadResCalcFile(strcat__m(buf, "\\xyzvectorE"), 
								PointList, 
								boundsVforE, 
								gResCalcPoint1, 
								gResCalcPoint2, 
								boundsV.beginP==boundsV.endP))!=0)
			return retc;

		strcpy__m(buf, path);
		if ((retc=ReadResCalcFile(strcat__m(buf, "\\xyzvectorB"), 
								PointList, 
								boundsB, 
								gResCalcPoint1, 
								gResCalcPoint2, 
								(boundsV.beginP==boundsV.endP)&&
								(boundsVforE.beginP==boundsVforE.endP)))!=0)
			return retc;

		if ((retc=ReadResCalcFile(strcat__m(buf, "\\xyzvectorA"), 
								PointList, 
								boundsA, 
								gResCalcPoint1, 
								gResCalcPoint2, 
								(boundsV.beginP==boundsV.endP)&&
								(boundsVforE.beginP==boundsVforE.endP)&&
								(boundsB.beginP==boundsB.endP)))!=0)
			return retc;


			gResCalcPoint1=PointXYZ(pnt[0].x-1e-3, pnt[0].y-1e-3, pnt[0].z-1e-3);
			gResCalcPoint2=PointXYZ(pnt[0].x+1e-3, pnt[0].y+1e-3, pnt[0].z+1e-3);
			for(i=1;i<kpnt;i++){
				if(pnt[i].x<gResCalcPoint1.x)gResCalcPoint1.x=pnt[i].x-1e-3;
				if(pnt[i].y<gResCalcPoint1.y)gResCalcPoint1.y=pnt[i].y-1e-3;
				if(pnt[i].z<gResCalcPoint1.z)gResCalcPoint1.z=pnt[i].z-1e-3;
				if(pnt[i].x>gResCalcPoint2.x)gResCalcPoint2.x=pnt[i].x+1e-3;
				if(pnt[i].y>gResCalcPoint2.y)gResCalcPoint2.y=pnt[i].y+1e-3;
				if(pnt[i].z>gResCalcPoint2.z)gResCalcPoint2.z=pnt[i].z+1e-3;
			}

		if (ResCalcPointNum=PointList.GetListLength())
		{
			ResCalcPointArray=new PointForResultCalculationInMesh3D[ResCalcPointNum];
			if (!ResCalcPointArray) return RETCODE_NOMEM;
			PointList.LoadList(ResCalcPointArray);
		}

		for (i=0; i<ResCalcPointNum; i++)
		{
			PointForResultCalculationInMesh3D& Pnti=ResCalcPointArray[i];
			Pnti.res=new PointXYZ[1]; //   !!!
			if (!Pnti.res) return RETCODE_NOMEM;
		}

		ListOfValues<RectForResultCalculationInMesh3D> RectList;

		for (j=0; j<ResCalcPointNum; j++)
		{
			PointForResultCalculationInMesh3D& pj=ResCalcPointArray[j];
			double lc[3],gc[3];

			gc[0]=pj.p.x;
			gc[1]=pj.p.y;
			gc[2]=pj.p.z;

			fuseair=((j>=boundsB.beginP && j<boundsB.endP) || (j>=boundsA.beginP && j<boundsA.endP));

			i=GetNearestElement(gc,lc,fuseair);
			if(i!=-1)
			{
				int t;
				PointXYZ P[8];
				const Rect3D &r=rect3d[i];
				const PointXYZ& p0=pnt[r.nodes[0]-1];
				const PointXYZ& p7=pnt[r.nodes[7]-1];
				loc_correct(lc[0],lc[1],lc[2]);
				pj.isfound=true;
				pj.ksi=lc[0];
				pj.eta=lc[1];
				pj.zet=lc[2];
				for(t=0;t<8;t++){P[t]=pnt[r.nodes[t]-1];}
				GetGlobalCoordinates(P,lc,gc);
				pj.p.x=gc[0];
				pj.p.y=gc[1];
				pj.p.z=gc[2];
				pj.ishex=true;
				pj.rnum=RectList.GetListLength();
				RectForResultCalculationInMesh3D RectForCalc;
				RectForCalc.rnum=i;		
				RectForCalc.hx=p7.x-p0.x;
				RectForCalc.hy=p7.y-p0.y;
				RectForCalc.hz=p7.z-p0.z;
				RectList.AddToList(RectForCalc);
				qrptotal+=1;
			}
			else
			{
				logfile<<j+1<<" "<<pj.p.x<<" "<<pj.p.y<<" "<<pj.p.z<<" "<<endl;
				cout<<j+1<<" "<<pj.p.x<<" "<<pj.p.y<<" "<<pj.p.z<<" "<<endl;
				return RETCODE_ERROR;
			}
		}

		if (ResCalcRectNum=RectList.GetListLength())
		{
			ResCalcRectArray=new RectForResultCalculationInMesh3D[ResCalcRectNum];
			if (ResCalcRectArray==NULL) return RETCODE_NOMEM;
			RectList.LoadList(ResCalcRectArray);
		}

		return RETCODE_OK;
	}

	int ReadOnlyPolarization(const char* path)
	{
		int retc, qrptotal=0;
		string buf;
		ifstream inf, rzf, nvtrf, nvkatf, l1f, wkf;

		strcpy__m(buf, path);
		if ((retc=ReadTimesFile(strcat__m(buf, "\\t")))!=0)
			return retc;

		qMN=0;

		ListOfValues<PointForResultCalculationInMesh3D> PointList;

		strcpy__m(buf, path);
		if ((retc=ReadResMNCalcFile(strcat__m(buf, "\\xyzMN"), 
								PointList, 
								boundsV, 
								gResCalcPoint1, 
								gResCalcPoint2, 
								true))!=0)
			return retc;


		return RETCODE_OK;
	}

	/*!      */
	int GetWeights(const double* q, const char* fname)
	{
		int i;
		for (i=0; i<nc; i++)
			q_all[i]=q[i];
		if (meshoptions&moUnloadSolutions)
		{
			char buf[256];
			ofstream outf;
			CHECKSTRING0(buf, strlen(taskpath)+32);
			sprintf_s(buf, "%s\\%s", taskpath, fname);
			outf.open(buf, ios::binary);
			for (i=0; i<kpnt; i++)
				outf < q_all[i];
			outf.close();
		}
		return RETCODE_OK;
	}

	/*!      */
	int GetWeights(const double* q, bool storeq=false)
	{
		int i;
		for (i=0; i<nc; i++)
		{
			q_all[i]=q[i];
			if (storeq)
				qv[i]=q_all[i];
		}

		return RETCODE_OK;
	}

	/*!     */
	int ReadMtrCatalog(const char* fname, double* cat, double* cat0, double vcoeff=1)
	{
		ifstream inf(fname);
		if (!inf) return RETCODE_NOFILE;
		int i;
		double v, v0;
		while (inf&&!inf.eof())
		{
			if (!inf.good()) return RETCODE_BADFILE;
			inf>>i>>v>>v0;
			inf.ignore(1000,'\n');
			if (!inf.eof()) { cat[i]=v*vcoeff; cat0[i]=v0*vcoeff; }
		}
		inf.close();
		return RETCODE_OK;
	}

	/*! c V0 */
	void CalcV0(int ipls)
	{
#ifndef FDLINE
		for (int i=0; i<nc; i++)
			qv_sum[i]=qva[i+kpnt*ipls]+qvb[i+kpnt*ipls]+qv[i];
#else
		for (int i=0; i<nc; i++)
			qv_sum[i]=qv[i];
#endif
	}

	/*!   xyzVector* */
	int ReadResCalcFile(const char* fname, ListOfValues<PointForResultCalculationInMesh3D>& _listP, ResCalcPointBounds& _boundsP,
						PointXYZ& gp1, PointXYZ& gp2, // gabarit points
						const bool& first=false)
	{
		int i, sz;
		ifstream xyzf(fname);
		_boundsP.beginP=_listP.GetListLength();
		if (xyzf)
		{
			xyzf>>sz;
			for (i=0; i<sz; i++)
			{
				PointForResultCalculationInMesh3D Pnti;
				if (!xyzf.good()) 
					return RETCODE_BADFILE;
				xyzf>>Pnti.p;
				if (first&&i==0)
				{
					gp1.x=gp2.x=Pnti.p.x;
					gp1.y=gp2.y=Pnti.p.y;
					gp1.z=gp2.z=Pnti.p.z;
				}
				else
				{
					gp1.x=min(gp1.x, Pnti.p.x);
					gp2.x=max(gp2.x, Pnti.p.x);
					gp1.y=min(gp1.y, Pnti.p.y);
					gp2.y=max(gp2.y, Pnti.p.y);
					gp1.z=min(gp1.z, Pnti.p.z);
					gp2.z=max(gp2.z, Pnti.p.z);
				}
				if (!_listP.AddToList(Pnti)) 
					return RETCODE_NOMEM;
			}
			xyzf.close();
		}
		_boundsP.endP=_listP.GetListLength();
		return RETCODE_OK;
	}

	/*!   xyzMN */
	int ReadResMNCalcFile(const char* fname, ListOfValues<PointForResultCalculationInMesh3D>& _listP, ResCalcPointBounds& _boundsP,
						PointXYZ& gp1, PointXYZ& gp2, // gabarit points
						const bool& first=false)
	{
		int i, sz;
		ifstream xyzf(fname);
		_boundsP.beginP=_listP.GetListLength();
		if (xyzf)
		{
			xyzf>>sz;
			qMN=sz;
			for (i=0; i<sz; i++)
			{
				PointForResultCalculationInMesh3D PntM, PntN, PntO;
				if (!xyzf.good()) 
					return RETCODE_BADFILE;
				xyzf>>PntM.p>>PntN.p;
				PntO.p=(PntM.p+PntN.p)*0.5;
				if (first&&i==0)
				{
					gp1.x=gp2.x=PntM.p.x;
					gp1.y=gp2.y=PntM.p.y;
					gp1.z=gp2.z=PntM.p.z;
					gp1.x=min(gp1.x, PntN.p.x);
					gp2.x=max(gp2.x, PntN.p.x);
					gp1.y=min(gp1.y, PntN.p.y);
					gp2.y=max(gp2.y, PntN.p.y);
					gp1.z=min(gp1.z, PntN.p.z);
					gp2.z=max(gp2.z, PntN.p.z);
				}
				else
				{
					gp1.x=min(gp1.x, PntM.p.x);
					gp2.x=max(gp2.x, PntM.p.x);
					gp1.y=min(gp1.y, PntM.p.y);
					gp2.y=max(gp2.y, PntM.p.y);
					gp1.z=min(gp1.z, PntM.p.z);
					gp2.z=max(gp2.z, PntM.p.z);
					gp1.x=min(gp1.x, PntN.p.x);
					gp2.x=max(gp2.x, PntN.p.x);
					gp1.y=min(gp1.y, PntN.p.y);
					gp2.y=max(gp2.y, PntN.p.y);
					gp1.z=min(gp1.z, PntN.p.z);
					gp2.z=max(gp2.z, PntN.p.z);
				}
				if (!_listP.AddToList(PntM)) 
					return RETCODE_NOMEM;
				if (!_listP.AddToList(PntN)) 
					return RETCODE_NOMEM;
				if (!_listP.AddToList(PntO)) 
					return RETCODE_NOMEM;
			}
			xyzf.close();
		}
		_boundsP.endP=_listP.GetListLength();
		return RETCODE_OK;
	}

	/*!    */
	int ReadTimesFile(const char* fname)
	{
		ifstream inf(fname);
		int i, sz;
		if (!inf) return RETCODE_NOFILE;
		inf >> sz; inf.ignore(1000, '\n');
		ntimes=sz;
		for (i=0; i<sz; i++)
		{
			if (!inf.good()) return RETCODE_BADFILE;
			inf >> times[i];
		}
		inf.close();
		return RETCODE_OK;
	}

	vector<LocalMatrix3D> vGx,vGy,vGz;
};
