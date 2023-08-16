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
 *  This file contains the headers for calculating the stationary task for primary field
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

extern bool isFileExists(const char *fname);

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

	void DisplaceUPointers()
	{
		VectorForSparseSLAE *pU2_copy;
		pU2_copy=pU2;
		pU2=pU1;
		pU1=pU;
		pU=pU2_copy;

		for (int i=0; i<n; i++)
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

// The class contains a grid for 2D scalar problems
struct MeshForVP
{
	const char* taskpath;	//!<     

	int meshoptions;		//!<  	

	PointXYZ A, B;			//!<   

	int kpnt, krect, kt1, nc;	//!< - , ,   1,  
	PointRZ* pnt;				//!<  
	Rect* rect;					//!<  
	int* l1;					//!<    1

	vector<PointRZ> dfPointA,dfPointB;	//!< 
	vector<int> dfNodeA,dfNodeB;		//!<     ,  0

	int c1, c2, c1b, c2b;		//!<     2

	double Rbak;				

	int ntimes;						//!< - 
	double times[MAXNUMBEROFTIMES];	//!< 

	double sigma[MAXNUMBEROFMATERIALS];	//!<  
	double sigmaZ[MAXNUMBEROFMATERIALS];	//!<  
	double currc;	//!<  

	int qMN;		//!< -  

	double *qv;		//!< 
	double *q_all;

	MeshForVP(int mOptions=0):
		taskpath("."),
		meshoptions(mOptions)
	{
		pnt=NULL;
		q_all=NULL;
		qv=NULL;
		rect=NULL;
		l1=NULL;
		ntimes=0;
		Rbak=0;
	}
	~MeshForVP()
	{
		if (pnt) delete [] pnt;
		if (q_all) delete [] q_all;
		if (qv) delete [] qv;
		if (rect) delete [] rect;
		if (l1) delete [] l1;
	}

	int Read(const char* path, int pntVsize, const PointXYZ* pntV, bool read_times=false)
	{
		int i, retc, qrptotal=0;
		string buf;
		ifstream inf, rzf, nvtrf, nvkatf, l1f, curf, abf, wkf;

#ifdef USEREGINFO
		ifstream infreg;
		int nreg, qr, qz;
		int *reg;
		double *rm, *zm;
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
		inf.clear();

		strcpy__m(buf, path);
		rzf.open(strcat__m(buf, "\\rz.dat"), ios::binary);
		if (!rzf) return RETCODE_NOFILE;
		pnt=new PointRZ[kpnt];
		if (!pnt) return RETCODE_NOMEM;
		for (i=0; i<kpnt; i++)
			rzf > pnt[i];
		rzf.close();
		rzf.clear();

		q_all=new double[kpnt];
		if (q_all==NULL) return RETCODE_NOMEM;

		strcpy__m(buf, path);
		nvtrf.open(strcat__m(buf, "\\nvtr.dat"), ios::binary);
		if (!nvtrf) return RETCODE_NOFILE;
		rect=new Rect[krect];
		if (!rect) return RETCODE_NOMEM;
		for (i=0; i<krect; i++)
			nvtrf > rect[i];
		nvtrf.close();
		nvtrf.clear();

		strcpy__m(buf, path);
		nvkatf.open(strcat__m(buf, "\\nvkat2d.dat"), ios::binary);
		if (!nvkatf) return RETCODE_NOFILE;
		for (i=0; i<krect; i++)
			nvkatf > rect[i].mtr;
		nvkatf.close();
		nvkatf.clear();

		strcpy__m(buf, path);
		l1f.open(strcat__m(buf, "\\l1.dat"), ios::binary);
		if (!l1f) return RETCODE_NOFILE;
		l1=new int[kt1];
		if (!l1) return RETCODE_NOMEM;
		for (i=0; i<kt1; i++)
			l1f > l1[i];
		l1f.close();
		l1f.clear();

		nc=kpnt;

		qv=new double[nc];
		if (qv==NULL) return RETCODE_NOMEM;

		for(i=0;i<MAXNUMBEROFMATERIALS;i++){sigma[i]=0.0;}

		strcpy__m(buf, path);
		if ((retc=ReadMtrCatalog(strcat__m(buf, "\\sigma"), sigma))!=0)
			return retc;

		buf=string(path) + string("\\sigmaZ");
		if(isFileExists(buf.c_str()))
		{
			strcpy__m(buf, path);
			if ((retc=ReadMtrCatalog(strcat__m(buf, "\\sigmaZ"), sigmaZ))!=0)
				return retc;
		}
		else
		{
			for(i=0;i<MAXNUMBEROFMATERIALS;i++){sigmaZ[i]=sigma[i];}
		}


		strcpy__m(buf, path);
		infreg.open(strcat__m(buf, "\\rz.txt"));
		if(!infreg) return RETCODE_NOFILE;
		infreg>>qr>>qz;
		infreg.close();
		infreg.clear();

		nreg=krect;

		if(!(reg=new int[nreg])) return RETCODE_NOMEM;
		if(!(rm=new double[qr])) return RETCODE_NOMEM;
		if(!(zm=new double[qz])) return RETCODE_NOMEM;

		strcpy__m(buf, path);
		infreg.open(strcat__m(buf, "\\r.dat"), ios::binary);
		if(!infreg) return RETCODE_NOFILE;
		for(i=0;i<qr;i++){infreg>rm[i];}
		infreg.close();
		infreg.clear();

		strcpy__m(buf, path);
		infreg.open(strcat__m(buf, "\\z.dat"), ios::binary);
		if(!infreg) return RETCODE_NOFILE;
		for(i=0;i<qz;i++){infreg>zm[i];}
		infreg.close();
		infreg.clear();

			for(i=0;i<nreg;i++){reg[i]=i+1;}

			Rbak=rm[qr-1];

			return RETCODE_OK;
		}

	/*!      */
	int GetWeights(const double* q, bool storeq=false)
	{
		int i, j;
		double sum;
		if (storeq)
		{
			for (i=0; i<nc; i++)
				qv[i]=q[i];
		}
		for (i=0; i<nc; i++)
			q_all[i]=q[i];

		if (meshoptions&moUnloadSolutions)
		{
			char buf[256];
			ofstream outf;
			CHECKSTRING0(buf, strlen(taskpath)+32);
			sprintf_s(buf, "%s\\v2.dat", taskpath);
			outf.open(buf, ios::binary);
			for (i=0; i<kpnt; i++)
				outf < q_all[i];
			outf.close();
		}
		return RETCODE_OK;
	}

	/*!     */
	int ReadMtrCatalog(const char* fname, double* cat)
	{
		ifstream inf(fname);
		if (!inf) return RETCODE_NOFILE;
		int i;
		double v;
		i=0;
		while (inf&&!inf.eof())
		{
			if (!inf.good()) return RETCODE_BADFILE;
			inf>>i;
			if (!inf.eof()) {
				inf>>v;
				cat[i]=v;
			}
		}
		inf.close();
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
