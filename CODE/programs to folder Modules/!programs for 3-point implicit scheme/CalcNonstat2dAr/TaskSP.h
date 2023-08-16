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
 *  This file contains the headers for calculating the radial source part of primary field for nonstationary task
 *  
 *  
 *  Written by Prof. Yuri G. Soloveichik, Prof. Marina G. Persova, Ph.D. Denis V. Vagin 
 *  Novosibirsk State Technical University,                                             
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                   
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                    
 *  Version 2                                                                           
*/

#pragma once
#include "ArrayOf.h"

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

bool CompareFiles(const char* _f1, const char* _f2);

#define DIAGNULL 1e-60
#define DENOMNULL 1e-60

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
	VectorForSparseSLAE U, U1, U2; // U(j), U(j-1), U(j-2)
	VectorForSparseSLAE *pU, *pU1, *pU2; // pointers to U(j), U(j-1), U(j-2)

	SparseSLAE()
	{
		n=jsize=0;
		iptr=NULL;
		jptr=NULL;
		pU=&U;
		pU1=&U1;
		pU2=&U2;

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

	void DisplaceUPointers(int nn)
	{
		VectorForSparseSLAE *pU2_copy;
		pU2_copy=pU2;
		pU2=pU1;
		pU1=pU;
		pU=pU2_copy;

		for (int i=0; i<nn; i++)
			pU->v[i]=pU1->v[i];
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
/*! @} */

// The class contains a grid for 2D scalar problems
struct Mesh
{
	const char* taskpath;	//!<     	

	int meshoptions;		//!<  	

	int kpnt, krect, kt1, nc;	//!< - , ,   1,  
	PointRZ* pnt;				//!<  
	Rect* rect;					//!<  
	int* l1;					//!<    1
	short* nvk1;				//!<    1
	
	int ntimes;						//!< - 
	double times[MAXNUMBEROFTIMES];	//!< 
	double currc[MAXNUMBEROFTIMES];	//!<  
	double fromT, toT;				//!<     

	double mu[MAXNUMBEROFMATERIALS];	//!<  
	double sigma[MAXNUMBEROFMATERIALS];	//!<  
	int mtr3d2d[MAXNUMBEROFMATERIALS];	//!<   3d 
	bool bmtr[MAXNUMBEROFMATERIALS];

	double *q_all, *q_u;	//!<      r

	DeltaFunction dFunc; //!< -

	bool TryToReadV;	//!<    

	int nreg, qr, qz;
	int *reg;
	double *rm, *zm;
	
	Mesh(int mOptions=0):
		taskpath("."),
		meshoptions(mOptions)
	{
		pnt=NULL;
		q_all=NULL;
		q_u=NULL;
		rect=NULL;
		l1=NULL;
		nvk1=NULL;
		ntimes=0;
		TryToReadV=false;
		
		reg=NULL;
		rm=zm=NULL;
		nreg=qr=qz=0;
	}
	~Mesh()
	{
		if (pnt) delete [] pnt;
		if (q_all) delete [] q_all;
		if (q_u) delete [] q_u;
		if (rect) delete [] rect;
		if (l1) delete [] l1;
		if (nvk1) delete [] nvk1;

		if (rm) delete [] rm;
		if (zm) delete [] zm;
	}

	int Read(const char* path)
	{
		int i, j, retc, qrptotal=0;
		string buf;
		ifstream inf, rzf, nvtrf, nvkatf, l1f, nvk1f;

#ifdef USEREGINFO
		ifstream infreg;
#endif

		taskpath=path;

		strcpy__m(buf, path);
		inf.open(strcat__m(buf, "\\inf2tr.dat"));
		if (!inf) return RETCODE_NOFILE;
		inf.ignore(1000, '\n');
		inf.ignore(1000, '='); inf>>kpnt;
		inf.ignore(1000, '='); inf>>krect;
		inf.ignore(1000, '='); inf>>kt1;
		inf.close();

		strcpy__m(buf, path);
		rzf.open(strcat__m(buf, "\\rz.dat"), ios::binary);
		if (!rzf) return RETCODE_NOFILE;
		pnt=new PointRZ[kpnt];
		if (!pnt) return RETCODE_NOMEM;
		for (i=0; i<kpnt; i++)
			rzf > pnt[i];
		rzf.close();

		q_all=new double[kpnt];
		if (q_all==NULL) return RETCODE_NOMEM;
		q_u=new double[kpnt];
		if (q_u==NULL) return RETCODE_NOMEM;

		strcpy__m(buf, path);
		nvtrf.open(strcat__m(buf, "\\nvtr.dat"), ios::binary);
		if (!nvtrf) return RETCODE_NOFILE;
		rect=new Rect[krect];
		if (!rect) return RETCODE_NOMEM;
		for (i=0; i<krect; i++)
			nvtrf > rect[i];
		nvtrf.close();

		for (i=0; i<MAXNUMBEROFMATERIALS; i++)
			bmtr[i]=false;

		strcpy__m(buf, path);
		nvkatf.open(strcat__m(buf, "\\nvkat2d.dat"), ios::binary);
		if (!nvkatf) return RETCODE_NOFILE;
		for (i=0; i<krect; i++)
		{
			Rect& ri=rect[i];
			nvkatf > ri.mtr;
			bmtr[ri.mtr]=true;
		}
		nvkatf.close();

		strcpy__m(buf, path);
		l1f.open(strcat__m(buf, "\\l1.dat"), ios::binary);
		if (!l1f) return RETCODE_NOFILE;
		l1=new int[kt1];
		if (!l1) return RETCODE_NOMEM;
		for (i=0; i<kt1; i++)
			l1f > l1[i];
		l1f.close();


		nc=kpnt;


		strcpy__m(buf, path);
		if ((retc=ReadTimesFile(strcat__m(buf, "\\infite.0")))!=0)
			return retc;

		strcpy__m(buf, path);
		if ((retc=ReadCurrentFunctionFile(strcat__m(buf, "\\CurrentFunction")))!=0)
			return retc;

		strcpy__m(buf, path);
		if ((retc=ReadTimesLimitsFile(strcat__m(buf, "\\TimeIntervalForPrint")))!=0)
			return retc;


		if (!(meshoptions&moForDipole)&&!(meshoptions&moForDisk))
		{
			if ((retc=dFunc.Read(path))!=0)
				return retc;
			dFunc.dfNode=FindNearestNode(dFunc.dfPoint);
		}


		strcpy__m(buf, path);
		if ((retc=ReadMtrCatalog(strcat__m(buf, "\\mu"), mu, MU0))!=0)
			return retc;

		strcpy__m(buf, path);
		if ((retc=ReadMtrCatalog(strcat__m(buf, "\\sigma"), sigma))!=0)
			return retc;

		if (meshoptions&moUsePntMtrInfo)
		{
			strcpy__m(buf, path);
			if ((retc=ReadMtrCatalog(strcat__m(buf, "\\mtr3d2d"), mtr3d2d))!=0)
				return retc;
		}

		buf=path;
		TryToReadV=true;

		strcpy__m(buf, path);
		infreg.open(strcat__m(buf, "\\rz.txt"));
		if (!infreg) return RETCODE_NOFILE;
		infreg>>qr>>qz;
		if ((rm=new double[qr])==NULL) return RETCODE_NOMEM;
		if ((zm=new double[qz])==NULL) return RETCODE_NOMEM;
		infreg.close();
		infreg.clear();

		strcpy__m(buf, path);
		infreg.open(strcat__m(buf, "\\r.dat"), ios::binary);
		for (i=0; i<qr; i++)
			infreg>rm[i];
		infreg.close();
		infreg.clear();

		strcpy__m(buf, path);
		infreg.open(strcat__m(buf, "\\z.dat"), ios::binary);
		for (i=0; i<qz; i++)
			infreg>zm[i];
		infreg.close();
		infreg.clear();
		
		return RETCODE_OK;
	}

	/*!   xyzVector* */
	int ReadResCalcFile(const char* fname, vector<PointForResultCalculation>& _listP, ResCalcPointBounds& _boundsP,
						PointRZ& gp1, PointRZ& gp2, // gabarit points
						const bool& first=false)
	{
		int i, sz;
		_boundsP.beginP=static_cast<int>(_listP.size());
		if (meshoptions&moDontReadResPoints)
		{
			_boundsP.endP=static_cast<int>(_listP.size());
			return RETCODE_OK;
		}
		ifstream xyzf(fname);
		if (xyzf)
		{
			xyzf>>sz;
			for (i=0; i<sz; i++)
			{
				PointForResultCalculation Pnti;
				if (!xyzf.good()) 
					return RETCODE_BADFILE;
				xyzf>>Pnti.p;
				Pnti.r=sqrt(Pnti.p.x*Pnti.p.x+Pnti.p.y*Pnti.p.y);
				if (Pnti.r<1e-3) Pnti.r=1e-3;
				Pnti.sinfi=Pnti.p.y/Pnti.r;
				Pnti.cosfi=Pnti.p.x/Pnti.r;
#ifndef USEREGINFO
				if (first&&i==0)
				{
					gp1.r=gp2.r=Pnti.r;
					gp1.z=gp2.z=Pnti.p.z;
				}
				else
				{
					gp1.r=min(gp1.r, Pnti.r);
					gp2.r=max(gp2.r, Pnti.r);
					gp1.z=min(gp1.z, Pnti.p.z);
					gp2.z=max(gp2.z, Pnti.p.z);
				}
#endif
				_listP.push_back(Pnti); 
			}
			xyzf.close();
		}
		_boundsP.endP=static_cast<int>(_listP.size());
		return RETCODE_OK;
	}

	/*!      */ 
	int GetWeights(const double* q, const int& tnum=0, const bool& unload=true)
	{
		int i, j;
		double sum;

		for (i=0; i<nc; i++)
			q_all[i]=q[i];
		
		if (meshoptions&moUnloadSolutions&&unload)
		{
			char buf[256];
			ofstream outf;
			CHECKSTRING0(buf, strlen(taskpath)+32);
			sprintf_s(buf, "%s\\v2.%d", taskpath, tnum);
			outf.open(buf, ios::binary);
			for (i=0; i<kpnt; i++)
				outf < q_all[i];
			outf.close();
		}
		
		return RETCODE_OK;
	}

	/*!   */
	int ReadSolution(const int &tnum, double* q_sol)
	{
		int i;
		char buf[256];
		ifstream inf;
		CHECKSTRING0(buf, strlen(taskpath)+32);
		sprintf_s(buf, "%s\\v2.%d", taskpath, tnum);
		inf.open(buf, ios::binary);
		if (!inf) return RETCODE_NOFILE;
		for (i=0; i<nc; i++)
		{
			if (!inf.good())
				return RETCODE_BADFILE;
			inf>q_sol[i];
			q_all[i]=q_sol[i];
		}
		for (i=nc; i<kpnt; i++)
		{
			if (!inf.good())
				return RETCODE_BADFILE;
			inf>q_all[i];
		}
		inf.close();
		inf.clear();
		return RETCODE_OK;
	}

	/*!      +    */
	int GetWeights(const double* q, const double* q1, const double* q2, const int& tnum)
	{

		int i, j, jj;
		double dt, dt0, dt1, mt, mt1, mt2, sum;
		
		if (tnum==0)
		{
			dt=times[2]-times[0];
			dt1=times[1]-times[0];
			dt0=times[2]-times[1];
			mt=-dt1/(dt*dt0);
			mt1=dt/(dt0*dt1);
			mt2=-(dt+dt1)/(dt*dt1);
		}
		else
		{
			if (tnum==1)
			{
				dt=times[2]-times[0];
				dt1=times[1]-times[0];
				dt0=times[2]-times[1];
				mt=dt1/(dt*dt0);
				mt1=(dt0-dt1)/(dt0*dt1);
				mt2=-dt0/(dt*dt1);
			}
			else
			{
				dt=times[tnum]-times[tnum-2];
				dt1=times[tnum-1]-times[tnum-2];
				dt0=times[tnum]-times[tnum-1];
				mt=(dt+dt0)/(dt*dt0);
				mt1=-dt/(dt0*dt1);
				mt2=dt0/(dt*dt1);
			}
		}

		for (i=0; i<nc; i++)
			q_all[i]=mt*q[i]+mt1*q1[i]+mt2*q2[i];

		return RETCODE_OK;
	}
	
	/*!    */
	int ReadTimesFile(const char* fname)
	{
		ifstream inf(fname);
		int i, sz;
		if (!inf) return RETCODE_NOFILE;
		inf.ignore(1000, '=');
		inf >> sz;
		if (sz>MAXNUMBEROFTIMES) 
			return RETCODE_ERROR;
		ntimes=sz;
		inf.ignore(1000, '\n');
		inf.ignore(1000, '\n');
		inf.ignore(1000, '\n');
		inf.ignore(1000, '\n');
		for (i=0; i<sz; i++)
		{
			if (!inf.good()) return RETCODE_BADFILE;
			inf >> times[i];
			inf.ignore(1000, ';');
		}
		inf.close();
		return RETCODE_OK;
	}

	/*!     */
	int ReadCurrentFunctionFile(const char* fname)
	{
		int i;
		double cft;
		ifstream inf(fname);
		if (!inf) return RETCODE_NOFILE;
		for (i=0; i<ntimes; i++)
		{
			if (!inf.good()) return RETCODE_BADFILE;
			inf >> cft >> currc[i];
		}
		inf.close();
		return RETCODE_OK;
	}

	/*!      */
	int ReadTimesLimitsFile(const char* fname)
	{
		ifstream inf(fname);
		if (!inf) return RETCODE_NOFILE;
		inf >> fromT >> toT;
		inf.close();
		return RETCODE_OK;
	}

	/*!     */
	int ReadMtrCatalog(const char* fname, double* cat, double vcoeff=1)
	{
		ifstream inf(fname);
		if (!inf) return RETCODE_NOFILE;
		int i;
		double v;
		while (inf&&!inf.eof())
		{
			if (!inf.good()) return RETCODE_BADFILE;
			inf>>i>>v;
			inf.ignore(1000,'\n');
			if (!inf.eof()) cat[i]=v*vcoeff;
		}
		inf.close();
		return RETCODE_OK;
	}

	/*!     */
	int ReadMtrCatalog(const char* fname, int* cat)
	{
		ifstream inf(fname);
		if (!inf) return RETCODE_NOFILE;
		int i, v;
		while (inf&&!inf.eof())
		{
			if (!inf.good()) return RETCODE_BADFILE;
			inf>>i>>v;
			inf.ignore(1000,'\n');
			if (!inf.eof()) cat[i]=v;
		}
		inf.close();
		return RETCODE_OK;
	}

	/*!    */
	int FindNearestNode(const PointRZ& forP)
	{
		int i, NearestNode;
		double dmin, d;
		for (i=0; i<kpnt; i++)
		{
			const PointRZ& pnti=pnt[i];
			d=sqrt((pnti.r-forP.r)*(pnti.r-forP.r)+(pnti.z-forP.z)*(pnti.z-forP.z));
			if (i==0) 
			{ 
				dmin=d; 
				NearestNode=i;
				continue;
			}
			if (d<dmin)
			{
				dmin=d; 
				NearestNode=i;
			}
		}
		return NearestNode;
	}
};

/*!    */
struct MeshForDisk
{

	Mesh BilinearRZMesh;		//!<  ,  
	int InitializationRetCode;	//!<     


	bool fsdiff;


	MeshForDisk(const char* path):
		BilinearRZMesh(moForDisk|moUsePntMtrInfo|moUnloadSolutions)
	{
		if ((InitializationRetCode=BilinearRZMesh.Read(path))!=0)
			return;
		BilinearRZMesh.taskpath=path;
		if (BilinearRZMesh.ntimes<2)
		{
			InitializationRetCode=RETCODE_ERROR;
			return;
		}
	}
	~MeshForDisk()
	{
	}
};
