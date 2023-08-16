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
 *  This file contains the code for calculating the eleptic source part of primary field for nonstationary two-layer scheme task
 *  
 *  
 *  Written by Prof. Yuri G. Soloveichik, Prof. Marina G. Persova, Ph.D. Denis V. Vagin 
 *  Novosibirsk State Technical University,                                             
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                   
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                    
 *  Version 2                                                                           
*/

#include "stdafx.h"

#include "TaskSP.h"

#include "pardiso.h"

#include "utils_direct_solver.h"

using namespace std;

extern int CalculateAllDecades;

struct SqLoop
{
	double Ax,Ay,Az;
	double Bx,By,Bz;
};

/*!       */
ofstream logfile;


extern int GetStiffnessMatrix(const double &r0, const double &z0, // coordinates of local node 2
							  const double &r3, const double &z3, // coordinates of local node 1
							  double m[4][4]);

extern int GetStiffnessMatrix0(const double &r0, const double &z0, // coordinates of local node 2
							   const double &r3, const double &z3, // coordinates of local node 1
							   double m[4][4]);

extern int GetMassMatrix(const double &r0, const double &z0, // coordinates of local node 2
						 const double &r3, const double &z3, // coordinates of local node 1
						 double m[4][4]);

extern int GetMatrixDr(const double &r0, const double &z0, // coordinates of local node 2
					 const double &r3, const double &z3, // coordinates of local node 1
					 double m[4][4]);

extern int GetMatrixDz(const double &r0, const double &z0, // coordinates of local node 2
					 const double &r3, const double &z3, // coordinates of local node 1
					 double m[4][4]);

extern int GetStiffnessMatrixForRect3D(
	const double &x1, const double &y1, const double &z1, // coordinates of local node 0
	const double &x2, const double &y2, const double &z2, // coordinates of local node 7
	double m[8][8]);

extern void GradMatrix_X(double Dy[8][8], double hx, double hy, double hz);
extern void GradMatrix_Y(double Dy[8][8], double hx, double hy, double hz);
extern void GradMatrix_Z(double Dy[8][8], double hx, double hy, double hz);

extern void Unpuk(double *uq,int npls,int slae_n,int tmap_n);

double Scal(double *v1, double *v2, int n)
{
	double s=0;
	int i;
	for (i=0; i<n; i++)
		s+=v1[i]*v2[i];
	return s;
}

/*!    */
#define GetErrorAndReturn(mes, retc) { GetError(mes, retc); return retc; }

extern bool CheckStop(void);

/*!    */
void GetError(const char* msg, int retc)
{
	printf("%s failed. Return code = %d.\n", msg, retc);
}

ifstream& operator>(ifstream& inf, PointRZ& p)
{
	inf > p.r > p.z;
	return inf;
}

PointXYZ operator*(const PointXYZ& p, const double& a)
{
	return PointXYZ(p.x*a, p.y*a, p.z*a);
}

PointXYZ operator+(const PointXYZ& p1, const PointXYZ& p2)
{
	return PointXYZ(p1.x+p2.x, p1.y+p2.y, p1.z+p2.z);
}

ifstream& operator>>(ifstream& inf, PointXYZ& p)
{
	inf >> p.x >> p.y >> p.z;
	return inf;
}

ifstream& operator>(ifstream& inf, PointXYZ& p)
{
	inf > p.x > p.y > p.z;
	return inf;
}

ofstream& operator<<(ofstream& outf, const PointXYZ& p)
{
	outf << p.x << " " << p.y << " " << p.z;
	return outf;
}

ofstream& operator<(ofstream& outf, const PointXYZ& p)
{
	outf < p.x < p.y < p.z;
	return outf;
}

ifstream& operator>(ifstream& inf, Rect& r)
{
	inf > r.nodes[2] > r.nodes[3] > r.nodes[0] > r.nodes[1] > r.nodes[4] > r.rtype;
	return inf;
}

void CopyV(const int &n, const double *from, double *to)
{
	int i;
	for (i=0; i<n; i++)
		to[i]=from[i];
}

double Scal(const int &n, const double *v1, const double *v2)
{
	double s=0;
	int i;
	for (i=0; i<n; i++)
		s+=v1[i]*v2[i];
	return s;
}

int MatrixOnVector(
	const double *x,	// z=Ax		// IN 
	const double *adiag,// matrix	//
	const double *altr,	//			//
	const double *autr,	//			//
	const int *iptr,	// portrait	//
	const int *jptr,	//			//
	const int &n,		// size		//
	double *z)			//			// OUT
{
	int i, j, jj;
	for (i=0; i<n; i++) 
		z[i]=x[i]*adiag[i];
	for (i=0; i<n; i++)  
		for (j=iptr[i]-1; j<iptr[i+1]-1; j++)
		{
			jj=jptr[j]-1;
			z[i]+=x[jj]*altr[j];
			z[jj]+=x[i]*autr[j];
		}
	return RETCODE_OK;
}

/*!   */
struct ListEl 
{
	int v;
	ListEl *next;
	ListEl() { next=NULL; };
};

/*!      */
class PortraitAL 
{
private:
	bool mem; //!<   
	int n;		//!<  
	ListEl **ListArray;	//!<  
	void DeleteList(ListEl* pListEl)
	{
		ListEl *cur=pListEl;
		ListEl *tmp;
		while (cur)
		{
			tmp=cur->next;
			delete cur;
			cur=tmp;
		}
	}
public:
	PortraitAL(int size)
	{
		ListArray=NULL;
		mem=((ListArray=new ListEl*[n=size])!=NULL);
		if (!mem) return;
		for (int i=0; i<n; i++) 
		{
			mem=((ListArray[i]=new ListEl)!=NULL);
			if (!mem) return;
			ListArray[i]->v=0;
		}
	}
	~PortraitAL()
	{
		if (ListArray)
		{
			for (int i=0; i<n; i++)
				DeleteList(ListArray[i]);
			delete [] ListArray;
		}
	}
	bool isMemoryAllocated()
	{
		return mem;
	}
	/*!      */
	int AddToPortrait(const int& isize, const int *ind)
	{
		bool iter;
		int el, i, j;
		ListEl *p, *t, *eln, *s;
		
		for (i=0; i<isize; i++)
			for (j=i+1; j<isize; j++) 
			{
				if (ind[i]>ind[j])
				{
					el=ind[j]; 
					p=ListArray[ind[i]-1]; 
				}
				else
				{
					el=ind[i]; 
					p=ListArray[ind[j]-1]; 
				}
				s=p;
				t=p->next;
				if (!t) 
				{ 
					p->next=new ListEl;
					if (p->next==NULL)
						return RETCODE_NOMEM;
					p->next->v=el;
					s->v++;
					continue;
				}
				iter=true;
				do 
				{
					if (el==t->v) 
						iter=false;
					if (el<t->v) 
					{
						eln=new ListEl;
						if (eln==NULL)
							return RETCODE_NOMEM;
						eln->v=el;
						eln->next=t;
						p->next=eln;
						s->v++;
						iter=false;
					}
					p=p->next;
					t=p->next;
				}
				while ((t)&&(iter));
				if (iter) 
				{
					p->next=new ListEl;
					if (p->next==NULL)
						return RETCODE_NOMEM;
					p->next->v=el;
					s->v++;
				}
			}
		return RETCODE_OK;
	}
	/*!      */
	int AddToPortrait(const int& bn, const int& isize, const int *ind)
	{
		bool iter;
		int el, j;
		ListEl *p, *t, *eln, *s;

		for (j=0; j<isize; j++) 
		{
			if (bn<=ind[j]) continue;
			el=ind[j];
			p=ListArray[bn-1]; 
			s=p;
			t=p->next;
			if (!t) 
			{ 
				p->next=new ListEl;
				if (p->next==NULL)
					return RETCODE_NOMEM;
				p->next->v=el;
				s->v++;
				continue;
			}
			iter=true;
			do 
			{
				if (el==t->v) 
					iter=false;
				if (el<t->v) 
				{
					eln=new ListEl;
					if (eln==NULL)
						return RETCODE_NOMEM;
					eln->v=el;
					eln->next=t;
					p->next=eln;
					s->v++;
					iter=false;
				}
				p=p->next;
				t=p->next;
			}
			while ((t)&&(iter));
			if (iter) 
			{
				p->next=new ListEl;
				if (p->next==NULL)
					return RETCODE_NOMEM;
				p->next->v=el;
				s->v++;
			}
		}
			return RETCODE_OK;
	}
	/*!   */
	int BuildPortrait(SparseSLAE* slae)
	{
		ListEl *p;
		int j=0, k;
		if ((slae->iptr=new int[n+1])==NULL)
			return RETCODE_NOMEM;
		slae->iptr[0]=1;
		for (k=0; k<n; k++)
			slae->iptr[k+1]=slae->iptr[k]+ListArray[k]->v;
		if ((slae->jptr=new int[(slae->jsize=slae->iptr[n]-slae->iptr[0])])==NULL)
			return RETCODE_NOMEM;
		for (k=0; k<n; k++) 
		{
			p=ListArray[k]->next;
			while (p) 
			{
				slae->jptr[j++]=p->v;
				p=p->next;
			}
		}
		if (j!=slae->jsize)
			return RETCODE_ERROR;
/*		
		ofstream out_iptr("iptr.txt");
		for (k=0; k<=n; k++)
			out_iptr<<slae->iptr[k]<<endl;
		out_iptr.close();

		ofstream out_jptr("jptr.txt");
		for (k=0; k<slae->jsize; k++)
			out_jptr<<slae->jptr[k]<<endl;
		out_jptr.close();
*/
		return RETCODE_OK;
	}
};
/*! @} */


inline int indX(const int& i)
{
	return div(i, 2).rem;
}

inline int indY(const int& i)
{
	return div(div(i, 2).quot, 2).rem;
}

inline int indZ(const int& i)
{
	return div(i, 4).quot;
}


extern double G1d[2][2];
extern double M1d[2][2];

void Mult_Plot(double *a, double *x, double *y, int n)
{
	int i, j, temp;
	double sum;

	for(i=0;i<n;i++)
	{
		sum = 0.0;
		temp = i*n;
		for(j=0;j<n;j++)
			sum += a[temp+j]*x[j];

		y[i] = sum;
	}
}

void AddV(const int &n, const double *from, double *to, const double& mlt=1)
{
	for (int i=0; i<n; i++)
		to[i]+=mlt*from[i];
}

/*!   2d  */
int BuildPortrait(const Mesh* ptrMesh, SparseSLAE* ptrSLAE, int npls)
{
	int i, j, k, retc, node, nextnode;
	PortraitAL P(ptrMesh->nc);
	ResizableArrayOf<int> indlist(32);
	for (i=0; i<ptrMesh->krect; i++)
	{
		const Rect& r=ptrMesh->rect[i];
		indlist.SetSize(0);
		for (j=0; j<4; j++)
		{
			node=r.nodes[j];
			if (!indlist.Find(node)) 
			{
				if (indlist.Add(node)==-1)
					return RETCODE_OUTOFRANGE;
			}
		}
		P.AddToPortrait(indlist.GetSize(), indlist.val);
	}
	ptrSLAE->n=ptrMesh->nc;
	if ((retc=P.BuildPortrait(ptrSLAE))!=0) 
		return retc;

	if ((retc=ptrSLAE->A.Init(ptrSLAE->n, ptrSLAE->jsize))!=0)
		return retc;
	if ((retc=ptrSLAE->B.Init(ptrSLAE->n, ptrSLAE->jsize))!=0)
		return retc;
	if ((retc=ptrSLAE->C.Init(ptrSLAE->n, ptrSLAE->jsize))!=0)
		return retc;
	if ((retc=ptrSLAE->C0.Init(ptrSLAE->n, ptrSLAE->jsize))!=0)
		return retc;
	if ((retc=ptrSLAE->Dr.Init(ptrSLAE->n, ptrSLAE->jsize))!=0)
		return retc;
	if ((retc=ptrSLAE->Dz.Init(ptrSLAE->n, ptrSLAE->jsize))!=0)
		return retc;
	if ((retc=ptrSLAE->F.Init(ptrSLAE->n*npls))!=0)
		return retc;
	if ((retc=ptrSLAE->G.Init(ptrSLAE->n*npls))!=0)
		return retc;
	if ((retc=ptrSLAE->tmpV.Init(ptrSLAE->n*npls))!=0)
		return retc;
	if ((retc=ptrSLAE->U.Init(ptrSLAE->n*npls))!=0)
		return retc;
	if ((retc=ptrSLAE->U1.Init(ptrSLAE->n*npls))!=0)
		return retc;
	if ((retc=ptrSLAE->U2.Init(ptrSLAE->n*npls))!=0)
		return retc;
	
	return RETCODE_OK;
}

/*!      */
void LocalMatrixOnVector(const double m[4][4], const double v[4], double res[4])
{
	int i, j;
	for (i=0; i<4; i++)
	{
		res[i]=0;
		for (j=0; j<4; j++)
			res[i]+=m[i][j]*v[j];
	}
}

/*!    2d  */
int GetLocalContributions(const Mesh* ptrMesh, SparseSLAE* ptrSLAE, bool withD=false, bool b0=false)
{
	int i, j, k, _i, _j, _ti, _tj, c, d;
	double lm[4][4];
	double lmm[4][4];
	double lmdr[4][4];
	double lmdz[4][4];
	double lmij, lmmij, el_mu, el_sigma, lmdrij, lmdzij;

	for (k=0; k<ptrMesh->krect; k++)
	{
		const Rect& r=ptrMesh->rect[k];
		const PointRZ& p2=ptrMesh->pnt[r.nodes[2]-1];
		const PointRZ& p1=ptrMesh->pnt[r.nodes[1]-1];
		if (b0)
			GetStiffnessMatrix0(p2.r, p2.z, p1.r, p1.z, lm);
		else
		GetStiffnessMatrix(p2.r, p2.z, p1.r, p1.z, lm);
		GetMassMatrix(p2.r, p2.z, p1.r, p1.z, lmm);
		if (withD)
		{
			GetMatrixDr(p2.r, p2.z, p1.r, p1.z, lmdr);
			GetMatrixDz(p2.r, p2.z, p1.r, p1.z, lmdz);
		}
		el_mu=1/(ptrMesh->mu[r.mtr]);
		el_sigma=ptrMesh->sigma[r.mtr];
		for (i=0; i<4; i++)
			for (j=0; j<4; j++)
			{
				_i=r.nodes[i];
				_j=r.nodes[j];
				lmij=lm[i][j];
				lmmij=lmm[i][j];
				if (withD)
				{
					lmdrij=lmdr[i][j];
					lmdzij=lmdz[i][j];
				}
               	ptrSLAE->AddToMatrix(_i-1, _j-1, el_mu*lmij, ptrSLAE->B);
				ptrSLAE->AddToMatrix(_i-1, _j-1, el_sigma*lmmij, ptrSLAE->C);
				if (withD)
				{
					ptrSLAE->AddToMatrix(_i-1, _j-1, lmmij, ptrSLAE->C0);
					ptrSLAE->AddToMatrix(_i-1, _j-1, lmdrij, ptrSLAE->Dr);
					ptrSLAE->AddToMatrix(_i-1, _j-1, lmdzij, ptrSLAE->Dz);
				}
			}
		
	}

	return RETCODE_OK;
}

/*!    2d  */
int ProcessBoundaryConditions(Mesh* ptrMesh, SparseSLAE* ptrSLAE, int npls)
{
	int ipls;
	for (int i=0; i<ptrMesh->kt1; i++)
	{
		const int& l1i=ptrMesh->l1[i];
		if (l1i>ptrMesh->nc) continue;
		const PointRZ& p=ptrMesh->pnt[l1i-1];
        ptrSLAE->SetMatrixDiagElement(l1i-1, 1e30, ptrSLAE->A);
		for(ipls=0;ipls<npls;ipls++)
		{
			ptrSLAE->SetVectorElement(ipls*ptrSLAE->n+(l1i-1), 0, ptrSLAE->F);
		}
	}
	return RETCODE_OK;
}

int ReadField2d(double *u,int n,int ind)
{
	FILE *fp;
	char buff[256];
	if(ind>=0)
	{
		sprintf(buff,"Ax\\v2.%d",ind);
		cout<<"Reading file "<<buff<<endl;
		fp=fopen(buff,"rb");
		if(!fp){return 1;}
		fread(u,sizeof(double),n,fp);
		fclose(fp);
		fflush(fp);
	}
	return 0;
}

int ReadField2d(double *u,int n,char *fname)
{
	FILE *fp;
	char buff[256];
	sprintf(buff,"Ax\\%s",fname);
	fp=fopen(buff,"rb");
	if(!fp){return 1;}
	cout<<"Reading file "<<fname<<endl;
	fread(u,sizeof(double),n,fp);
	fclose(fp);
	fflush(fp);
	return 0;
}

int WriteField2d(double *u,int n,int ind)
{
	FILE *fp;
	char buff[256];
	if(ind>=0)
	{
		sprintf(buff,"Ax\\v2.%d",ind);
		fp=fopen(buff,"wb");
		if(!fp){return 1;}
		fwrite(u,sizeof(double),n,fp);
		fclose(fp);
		fflush(fp);
	}
	return 0;
}

int WriteField2d(double *u,int n,char *fname)
{
	FILE *fp;
	char buff[256];
	cout<<"Writing file "<<fname<<endl;
	sprintf(buff,"Ax\\%s",fname);
	fp=fopen(buff,"wb");
	if(!fp){return 1;}
	fwrite(u,sizeof(double),n,fp);
	fclose(fp);
	fflush(fp);
	return 0;
}

/* 2d      () */
int SolveNonStationarProblemForSpiderWithoutResCalcPoints_2(const char* path)
{
	int i, k, indt, retc;

	SparseSLAE SLAE;

	double dt, dt1, dt0;

	ifstream inf;
	ofstream ofp;

	pardiso_solver prds;

	int p1,p2,npls,nprof;
	
	int ipls;
	vector<int> nodeDeltaFuncA,nodeDeltaFuncB;
	vector<SqLoop> GenSq;
	vector<double> currentDeltaFunc;
	
	vector<int> DecIg;

	int iDec,nDec,npre1,npre2;

	int *ig,*jg;

	int ntime;
	double *time;

	int tbeg,tend;

	bool fcalc;

	double ct1,ct2;

	ig=NULL;
	jg=NULL;
	time=NULL;

	GetNumberOfPlaces(npls);

	GenSq.resize(npls);
	nodeDeltaFuncA.resize(npls);
	nodeDeltaFuncB.resize(npls);
	currentDeltaFunc.resize(npls);



	inf.open("srsclca");
	if(!inf)
	{
		logfile<<"Error in open file "<<"srsclca"<<endl;
		cout<<"Error in open file "<<"srsclca"<<endl;
		return 1;
	}
	for(i=0;i<npls;i++){inf>>GenSq[i].Ax>>GenSq[i].Ay>>GenSq[i].Az>>GenSq[i].Bx>>GenSq[i].By>>GenSq[i].Bz;}
	inf.close();
	inf.clear();

	for(i=0;i<npls;i++){currentDeltaFunc[i]=1.0;}

	retc=read_time_mesh(ntime,time);
	if(retc)
	{
		logfile<<"Error in function "<<"read_time_mesh"<<endl;
		cout<<"Error in function "<<"read_time_mesh"<<endl;
		return retc;
	}

	MeshForSpiderWithoutResCalcPoints_2 SMesh(path);

	SMesh.BilinearRZMesh.dFunc.dfNodeA.resize(npls);
	SMesh.BilinearRZMesh.dFunc.dfNodeB.resize(npls);
	SMesh.BilinearRZMesh.dFunc.dfPointA.resize(npls);
	SMesh.BilinearRZMesh.dFunc.dfPointB.resize(npls);
	for(ipls=0;ipls<npls;ipls++)
	{
		SMesh.BilinearRZMesh.dFunc.dfPointA[ipls].r=0.0;
		SMesh.BilinearRZMesh.dFunc.dfPointA[ipls].z=GenSq[ipls].Az;

		SMesh.BilinearRZMesh.dFunc.dfPointB[ipls].r=0.0;
		SMesh.BilinearRZMesh.dFunc.dfPointB[ipls].z=GenSq[ipls].Bz;

		SMesh.BilinearRZMesh.dFunc.dfNodeA[ipls]=SMesh.BilinearRZMesh.FindNearestNode(SMesh.BilinearRZMesh.dFunc.dfPointA[ipls]);
		SMesh.BilinearRZMesh.dFunc.dfNodeB[ipls]=SMesh.BilinearRZMesh.FindNearestNode(SMesh.BilinearRZMesh.dFunc.dfPointB[ipls]);
	}

	if (SMesh.InitializationRetCode!=0)
		return SMesh.InitializationRetCode;


	if ((retc=BuildPortrait(&SMesh.BilinearRZMesh, &SLAE, npls))!=0)
		return retc;

	ig=new int[SLAE.n+1];
	jg=new int[SLAE.jsize];


	for(i=0;i<(SLAE.n+1);i++){ig[i]=SLAE.iptr[i]-1;}
	for(i=0;i<SLAE.jsize;i++){jg[i]=SLAE.jptr[i]-1;}

	if ((retc=GetLocalContributions(&SMesh.BilinearRZMesh, &SLAE, false, true))!=0)
		return retc;

	int giDec,gnDec;

	giDec=0;
	gnDec=0;
	
	do{
	_chdir(path);

	if(CalculateAllDecades)
	{
		if(!giDec)
		{
			GenerateDecInfo(ntime,time);
		}
		ofp.open("iDec");
		ofp<<giDec<<'\n';
		ofp.close();
		ofp.clear();
	}

	retc=ReadDecInfo(nDec,iDec,DecIg);
	if(retc)
	{
		logfile<<"Error in function "<<"ReadDecInfo"<<endl;
		cout<<"Error in function "<<"ReadDecInfo"<<endl;
		return retc;
	}

	_chdir("..");

	giDec=iDec;
	gnDec=nDec;

	npre1=-1;
	npre2=-1;
	if(iDec)
	{
		npre1=FindTimeLayer(ntime,time,time[DecIg[iDec]]-1*(time[DecIg[iDec]+1]-time[DecIg[iDec]]));
		npre2=FindTimeLayer(ntime,time,time[DecIg[iDec]]-2*(time[DecIg[iDec]+1]-time[DecIg[iDec]]));
		if(npre1<0)
		{
			cout<<"Error in npre1"<<endl;
			logfile<<"Error in npre1"<<endl;
			return 1;
		}
		if(npre2<0)
		{
			cout<<"Error in npre2"<<endl;
			logfile<<"Error in npre2"<<endl;
			return 1;
		}
	}

	cout<<"npre1= "<<npre1<<endl;
	cout<<"npre2= "<<npre2<<endl;

	tbeg=DecIg[iDec];
	tend=DecIg[iDec+1];

	cout<<"tbeg= "<<tbeg<<endl;
	cout<<"tend= "<<tend<<endl;

	SMesh.fsdiff=false;

	fcalc=false;

	if(0>=tbeg && 0<tend)
	{
	bool fstart;

	i=ReadField2d(SLAE.pU->v,SLAE.n*npls,"a2d.start");
	fstart=!i;


	SLAE.A.Clear(SLAE.n, SLAE.jsize);
	SLAE.F.Clear(SLAE.n*npls);

	SLAE.A.AddToMatrix(SLAE.n, SLAE.jsize, SLAE.B);

	if ((retc=ProcessBoundaryConditions(&SMesh.BilinearRZMesh, &SLAE, npls))!=0)
		return retc;

	for(ipls=0;ipls<npls;ipls++)
	{



	SLAE.F.v[SMesh.BilinearRZMesh.dFunc.dfNodeA[ipls]+SLAE.n*ipls]+=
			SMesh.BilinearRZMesh.currc[0]*currentDeltaFunc[ipls]/(2*_PI_);

	}

	if(!fstart)
	{
		logfile<<"factorize slae"<<'\n';
		prds.factorize(SLAE.n,ig,jg,SLAE.A.altr,SLAE.A.adiag,1);
		prds.solve_nrhs(npls,SLAE.F.v,SLAE.pU->v);
		prds.stop_solver();
	}
	

	WriteField2d(SLAE.pU->v,SLAE.n*npls,0);

	
	SLAE.DisplaceUPointers(SLAE.n*npls);

	cout<<"time layer "<<0<<" done"<<endl;

	fcalc=true;
	}

	bool SolverActive;

	SolverActive=false;

	bool FirstDecadeStep;

	FirstDecadeStep=true;

	for(indt=1; indt<SMesh.BilinearRZMesh.ntimes; indt++)
	{
		SMesh.fsdiff=false;

		if(indt>=tbeg && indt<tend)
		{
		if(!fcalc)
		{
			if((retc=ReadField2d(SLAE.pU1->v,SLAE.n*npls,npre1))!=0)return retc;
			SMesh.fsdiff=true;
		}

		SLAE.F.Clear(SLAE.n*npls);

		if(FirstDecadeStep)
		{
			if(npre1==-1)
			{
				npre1=indt-1;
			}

			SLAE.A.Clear(SLAE.n, SLAE.jsize);
			dt=SMesh.BilinearRZMesh.times[indt]-SMesh.BilinearRZMesh.times[npre1];
			SLAE.A.AddToMatrix(SLAE.n, SLAE.jsize, SLAE.B);
			SLAE.A.AddToMatrix(SLAE.n, SLAE.jsize, SLAE.C, 1.0/dt);
			if ((retc=ProcessBoundaryConditions(&SMesh.BilinearRZMesh, &SLAE, npls))!=0)
				return retc;
			prds.factorize(SLAE.n,ig,jg,SLAE.A.altr,SLAE.A.adiag,1);
			SolverActive=true;
			FirstDecadeStep=false;
		}

		for(ipls=0;ipls<npls;ipls++)
		{
			MatrixOnVector(SLAE.pU1->v+SLAE.n*ipls, SLAE.C.adiag, SLAE.C.altr, SLAE.C.autr, SLAE.iptr, SLAE.jptr, SLAE.n, SLAE.tmpV.v);
			for(i=0;i<SLAE.n;i++){SLAE.F.v[i+SLAE.n*ipls]+=SLAE.tmpV.v[i]/dt;}
		}

		for(ipls=0;ipls<npls;ipls++)
		{


		SLAE.F.v[SMesh.BilinearRZMesh.dFunc.dfNodeA[ipls]+SLAE.n*ipls]+=
				SMesh.BilinearRZMesh.currc[indt]*currentDeltaFunc[ipls]/(2*_PI_);

		}

		prds.solve_nrhs(npls,SLAE.F.v,SLAE.pU->v);


		WriteField2d(SLAE.pU->v,SLAE.n*npls,indt);

		if(indt==SMesh.BilinearRZMesh.ntimes-1)
		{
			WriteField2d(SLAE.pU->v,SLAE.n*npls,"a2d.last");
		}

		SLAE.DisplaceUPointers(SLAE.n*npls);

		cout<<"time layer "<<indt<<" done"<<endl;
		
		fcalc=true;
		}
	}

	if(SolverActive)
	{
		prds.stop_solver();
		SolverActive=false;
	}

		giDec++;
	}while(CalculateAllDecades && giDec<gnDec);


	if(ig){delete [] ig;ig=NULL;}
	if(jg){delete [] jg;jg=NULL;}
	if(time){delete [] time;time=NULL;}

	return RETCODE_OK;
}

int ProcTask2DLine2_2(string &tdir,char *tp=NULL)
{
	int retc;
	if ((retc=SolveNonStationarProblemForSpiderWithoutResCalcPoints_2(tdir.c_str()))!=0)
	{
		char str[1024];
		sprintf_s(str, "Solving nonstationar problem for spider failed. Error code = %d.", retc);
		throw logic_error(str);
	}
	return 0;
}

int CalcSP()
{
	string tdir;

	if(CheckStop())return 1;

	tdir="Ax";

	ProcTask2DLine2_2(tdir,"A");

	return 0;
}
