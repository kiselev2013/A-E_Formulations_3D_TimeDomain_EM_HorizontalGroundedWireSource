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
 *  This file contains the code for calculating the stationary task for primary field
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

using namespace std;

extern int FindIntervalInDoubleMassWithEps(double *vec,double elem,int size,double eps);

/*!       */
ofstream logfile;
bool fstop;

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

extern int GetStiffnessMatrix0Dz(const double &r0, const double &z0, // coordinates of local node 2
							   const double &r3, const double &z3, // coordinates of local node 1
							   double m[4][4]);

extern int GetStiffnessMatrix0Dr(const double &r0, const double &z0, // coordinates of local node 2
							   const double &r3, const double &z3, // coordinates of local node 1
							   double m[4][4]);

bool CheckStop(void)
{
	ifstream ifstop;
	ifstop.open("..\\stop");
	if(ifstop){
		fstop=true;
		ifstop.close();
	}
	ifstop.clear();
	return fstop;
}

int GetNumberOfPlaces(int &npls)
{
	int i,k,p1,p2,nprof;
	ifstream inf;

	npls=0;
	inf.open("clcnplsh");
	if(!inf)
	{
		logfile<<"Error in open file "<<"clcnplsh"<<endl;
		cout<<"Error in open file "<<"clcnplsh"<<endl;
		return 1;
	}
	inf>>npls;
	inf.close();
	inf.clear();

	return 0;
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

double Scal(double *x, double *y, int n)
{
	double sum = 0.0;
	long i;
	for(i=0;i<n;i++)
		sum += x[i]*y[i];
	return sum;
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

		return RETCODE_OK;
	}
};
/*! @} */


/*!   2d   */
int BuildPortrait(const MeshForVP* ptrMesh, SparseSLAE* ptrSLAE,int npls)
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
				if (indlist.Add(node)==-1)
					return RETCODE_OUTOFRANGE;

		}
		P.AddToPortrait(indlist.GetSize(), indlist.val);
	}
	ptrSLAE->n=ptrMesh->nc;
	if ((retc=P.BuildPortrait(ptrSLAE))!=0) 
		return retc;
	
	if ((retc=ptrSLAE->A.Init(ptrSLAE->n, ptrSLAE->jsize))!=0)
		return retc;
	if ((retc=ptrSLAE->B.Init(ptrSLAE->n, ptrSLAE->jsize))!=0) // Bij=sigma0*gradFIi*gradFIj
		return retc;
	if ((retc=ptrSLAE->C.Init(ptrSLAE->n, ptrSLAE->jsize))!=0) // Cij=gradFIi*gradFIj
		return retc;
	if ((retc=ptrSLAE->F.Init(ptrSLAE->n*npls))!=0)
		return retc;
	if ((retc=ptrSLAE->U.Init(ptrSLAE->n*npls))!=0)
		return retc;

	return RETCODE_OK;
}


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

/*!    2d   */
int GetLocalContributions(const MeshForVP* ptrMesh, SparseSLAE* ptrSLAE, int mtr) 
{
	int i, j, k, _i, _j, _ti, c, d;
	double lm[4][4];
	double lmij, el_sigma;
	for (k=0; k<ptrMesh->krect; k++)
	{
		const Rect& r=ptrMesh->rect[k];
		const PointRZ& p2=ptrMesh->pnt[r.nodes[2]-1];
		const PointRZ& p1=ptrMesh->pnt[r.nodes[1]-1];
		if (r.mtr!=mtr) continue;
		GetStiffnessMatrix0(p2.r, p2.z, p1.r, p1.z, lm);
		el_sigma=ptrMesh->sigma[r.mtr];
		for (i=0; i<4; i++)
			for (j=0; j<4; j++)
			{
				_i=r.nodes[i];
				_j=r.nodes[j];
				lmij=lm[i][j];
				ptrSLAE->AddToMatrix(_i-1, _j-1, el_sigma*lmij, ptrSLAE->C);
			}
	}
	return RETCODE_OK;
}

/*!    2d   */
int GetLocalContributions(const MeshForVP* ptrMesh, SparseSLAE* ptrSLAE)
{
	int i, j, k, _i, _j, _ti, c, d;
	double lmdr[4][4],lmdz[4][4];
	double lmij;
	for (k=0; k<ptrMesh->krect; k++)
	{
		const Rect& r=ptrMesh->rect[k];
		const PointRZ& p2=ptrMesh->pnt[r.nodes[2]-1];
		const PointRZ& p1=ptrMesh->pnt[r.nodes[1]-1];
		GetStiffnessMatrix0Dr(p2.r, p2.z, p1.r, p1.z, lmdr);
		GetStiffnessMatrix0Dz(p2.r, p2.z, p1.r, p1.z, lmdz);
		for (i=0; i<4; i++)
			for (j=0; j<4; j++)
			{
				_i=r.nodes[i];
				_j=r.nodes[j];
				lmij=lmdr[i][j]*(ptrMesh->sigma[r.mtr])+lmdz[i][j]*(ptrMesh->sigmaZ[r.mtr]);
				ptrSLAE->AddToMatrix(_i-1, _j-1, lmij, ptrSLAE->B);
			}
	}
	return RETCODE_OK;
}

/*!    2d   */
int ProcessBoundaryConditions(const MeshForVP* ptrMesh, SparseSLAE* ptrSLAE,int npls)
{
	int ipls;

	for (int i=0; i<ptrMesh->kt1; i++)
	{
		const int& l1i=ptrMesh->l1[i];
		if (l1i>ptrMesh->nc) continue;
		ptrSLAE->SetMatrixDiagElement(l1i-1, 1e30, ptrSLAE->A);
	}

	for(ipls=0;ipls<npls;ipls++)
	{
		for (int i=0; i<ptrMesh->kt1; i++)
		{
			const int& l1i=ptrMesh->l1[i];
			if (l1i>ptrMesh->nc) continue;
			ptrSLAE->SetVectorElement(l1i-1+ptrSLAE->n*ipls, 0, ptrSLAE->F);
		}
	}

	return RETCODE_OK;
}

/*!      2d   */
int ProcessCurrentBoundaryConditions(const MeshForVP* ptrMesh, SparseSLAE* ptrSLAE)
{
	const PointRZ& pc1=ptrMesh->pnt[ptrMesh->c1-1];
	const PointRZ& pc2=ptrMesh->pnt[ptrMesh->c2-1];
	if (ptrMesh->meshoptions&moForVEL)
	{
		ptrSLAE->AddToVector(ptrMesh->c1-1, (ptrMesh->currc)/(6*_PI_), ptrSLAE->F);
		ptrSLAE->AddToVector(ptrMesh->c2-1, (ptrMesh->currc)/(3*_PI_), ptrSLAE->F);
		ptrSLAE->AddToVector(ptrMesh->c1b-1, -(ptrMesh->currc)/(6*_PI_), ptrSLAE->F);
		ptrSLAE->AddToVector(ptrMesh->c2b-1, -(ptrMesh->currc)/(3*_PI_), ptrSLAE->F);
	}
	else
	{
		if (ptrMesh->meshoptions&moForCED)
		{
			const double& rk1=ptrMesh->pnt[ptrMesh->c1b-1].r;
			const double& rk2=ptrMesh->pnt[ptrMesh->c2b-1].r;
			const double hk=rk2-rk1;

			ptrSLAE->AddToVector(ptrMesh->c1-1, -(ptrMesh->currc)/(6*_PI_), ptrSLAE->F);
			ptrSLAE->AddToVector(ptrMesh->c2-1, -(ptrMesh->currc)/(3*_PI_), ptrSLAE->F);

			ptrSLAE->AddToVector(ptrMesh->c1b-1, (ptrMesh->currc)*hk*(3*rk1+  hk)/(6*_PI_*(rk2*rk2-rk1*rk1)), ptrSLAE->F);
			ptrSLAE->AddToVector(ptrMesh->c2b-1, (ptrMesh->currc)*hk*(3*rk1+2*hk)/(6*_PI_*(rk2*rk2-rk1*rk1)), ptrSLAE->F);
		}
		else
		{
			ptrSLAE->AddToVector(ptrMesh->c1-1, (ptrMesh->currc)/(6*_PI_), ptrSLAE->F);
			ptrSLAE->AddToVector(ptrMesh->c2-1, (ptrMesh->currc)/(3*_PI_), ptrSLAE->F);
		}
	}
	return RETCODE_OK;
}

void AddV(const int &n, const double *from, double *to, const double& mlt=1)
{
	for (int i=0; i<n; i++)
		to[i]+=mlt*from[i];
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

struct LineAB
{
	double Ax,Ay,Az;
	double Bx,By,Bz;
};

/*!   */
int SolveStationarProblem2d(const char* path,bool f3d, bool forVEL=false, bool forCED=false)
{
	ifstream inf;
	ofstream ofp;
	int retc,i,ipls,npls;
	pardiso_solver prds;
	vector<int> DecIg;
	int *ig,*jg;
	vector<LineAB> GenAB;
	vector<double> currentDeltaFunc;
	double eps_vel,ImpStatCur;

	MeshForVP meshvp(moUnloadSolutions|(forVEL?moForVEL:0)|(forCED?moForCED:0));
	SparseSLAE SLAE;
	
	GetNumberOfPlaces(npls);

	GenAB.resize(npls);
	currentDeltaFunc.resize(npls);

	inf.open("lc.txt");
	if(!inf)
	{
		logfile<<"Error in open file "<<"lc.txt"<<endl;
		cout<<"Error in open file "<<"lc.txt"<<endl;
		return 1;
	}
	inf>>i;
	inf>>eps_vel;
	inf.close();
	inf.clear();

	inf.open("srsclch");
	if(!inf)
	{
		logfile<<"Error in open file "<<"srsclch"<<endl;
		cout<<"Error in open file "<<"srsclch"<<endl;
		return 1;
	}
	for(ipls=0;ipls<npls;ipls++){inf>>GenAB[ipls].Ax>>GenAB[ipls].Ay>>GenAB[ipls].Az>>GenAB[ipls].Bx>>GenAB[ipls].By>>GenAB[ipls].Bz;}
	inf.close();
	inf.clear();

	for(ipls=0;ipls<npls;ipls++){currentDeltaFunc[ipls]=1.0;}

	inf.open("ImpStatCur");
	if(inf)
	{
		inf>>ImpStatCur;
		inf.close();
		for(ipls=0;ipls<npls;ipls++){currentDeltaFunc[ipls]*=ImpStatCur;}
	}
	inf.clear();

	if ((retc=meshvp.Read(path, 0, NULL))!=0)
		return retc;
	
	meshvp.dfNodeA.resize(npls);
	meshvp.dfNodeB.resize(npls);
	meshvp.dfPointA.resize(npls);
	meshvp.dfPointB.resize(npls);
	for(ipls=0;ipls<npls;ipls++)
	{
		meshvp.dfPointA[ipls].r=0.0;
		meshvp.dfPointA[ipls].z=GenAB[ipls].Az;

		meshvp.dfPointB[ipls].r=0.0;
		meshvp.dfPointB[ipls].z=GenAB[ipls].Bz;

		meshvp.dfNodeA[ipls]=meshvp.FindNearestNode(meshvp.dfPointA[ipls]);
		meshvp.dfNodeB[ipls]=meshvp.FindNearestNode(meshvp.dfPointB[ipls]);
	}

	if ((retc=BuildPortrait(&meshvp, &SLAE, npls))!=0)
		return retc;

	ig=new int[SLAE.n+1];
	jg=new int[SLAE.jsize];

	for(i=0;i<(SLAE.n+1);i++){ig[i]=SLAE.iptr[i]-1;}
	for(i=0;i<SLAE.jsize;i++){jg[i]=SLAE.jptr[i]-1;}


	if ((retc=GetLocalContributions(&meshvp, &SLAE))!=0)
		return retc;
	
	SLAE.A.Clear(SLAE.n, SLAE.jsize);
	SLAE.F.Clear(SLAE.n*npls);

	SLAE.A.AddToMatrix(SLAE.n, SLAE.jsize, SLAE.B);

	ofp.open("gens");
	ofp<<scientific<<setprecision(16);
	ofp<<npls<<'\n';
	for(ipls=0;ipls<npls;ipls++)
	{
		SLAE.F.v[meshvp.dfNodeA[ipls]+SLAE.n*ipls]+=currentDeltaFunc[ipls]/(2.0*_PI_);
		if(fabs(GenAB[ipls].Az-GenAB[ipls].Bz)>eps_vel)
		{
			SLAE.F.v[meshvp.dfNodeB[ipls]+SLAE.n*ipls]-=currentDeltaFunc[ipls]/(2.0*_PI_);
			ofp<<2<<' ';
		}
		else
		{
			ofp<<1<<' ';
		}
		ofp<<GenAB[ipls].Ax<<' '<<GenAB[ipls].Ay<<' '<<GenAB[ipls].Az<<' '<<GenAB[ipls].Bx<<' '<<GenAB[ipls].By<<' '<<GenAB[ipls].Bz<<'\n';
	}
	ofp.close();
	ofp.clear();

	if ((retc=ProcessBoundaryConditions(&meshvp, &SLAE, npls))!=0)
		return retc;

	prds.factorize(SLAE.n,ig,jg,SLAE.A.altr,SLAE.A.adiag,1);
	prds.solve_nrhs(npls,SLAE.F.v,SLAE.pU->v);
	prds.stop_solver();


	FILE *fp;
	fp=fopen("v2.dat","wb");
	fwrite(SLAE.pU->v,sizeof(double),meshvp.nc*npls,fp);
	fclose(fp);

	return RETCODE_OK;
}
int CalcSP()
{
	if(!CheckStop())
	{
		SolveStationarProblem2d(".", true, false, false);
	}
	return 0;
}
