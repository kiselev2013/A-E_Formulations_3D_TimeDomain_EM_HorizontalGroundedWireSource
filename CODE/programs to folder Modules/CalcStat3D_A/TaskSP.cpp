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
 *  This file contains the code for calculation of secondary part of field A for stationary task
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

#include "CheckInHex.h"


extern int FindIntervalInDoubleMassWithEps(double *vec,double elem,int size,double eps);

/*!       */
ofstream logfile;
bool fstop;


extern int GetStiffnessMatrixForRect3D(
	const double &x1, const double &y1, const double &z1, // coordinates of local node 0
	const double &x2, const double &y2, const double &z2, // coordinates of local node 7
	double m[8][8],double sigx=1,double sigy=1,double sigz=1);

extern void GradMatrix_X(double Dy[8][8], double hx, double hy, double hz);
extern void GradMatrix_Y(double Dy[8][8], double hx, double hy, double hz);
extern void GradMatrix_Z(double Dy[8][8], double hx, double hy, double hz);

void mult_matrix(double a[][3],double b[][3],double c[][3])
{
	int i,j,k;
	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			c[i][j]=0.0;
			for(k=0;k<3;k++)
			{
				c[i][j]+=a[i][k]*b[k][j];
			}
		}
	}
}

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

void loc_correct(double &ksi,double &eta,double &dzeta)
{
	ksi= (ksi<0.0)? 0.0 : (ksi>1.0)? 1.0 : ksi;
	eta= (eta<0.0)? 0.0 : (eta>1.0)? 1.0 : eta;
	dzeta= (dzeta<0.0)? 0.0 : (dzeta>1.0)? 1.0 : dzeta;
}

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

ifstream& operator>(ifstream& inf, Rect3D& r)
{
	for (int i=0; i<13; i++)
		inf > r.nodes[i];
	inf > r.rtype;
	return inf;
}

ifstream& operator>(ifstream& inf, Rect& r)
{
	inf > r.nodes[2] > r.nodes[3] > r.nodes[0] > r.nodes[1] > r.nodes[4] > r.rtype;
	return inf;
}

void GetNodeSolutionByLocalcCoords(double *lc,const double *q,double &Val)
{
	Val=q[0]*(1-lc[0])*(1-lc[1])*(1-lc[2])+q[1]*(lc[0])*(1-lc[1])*(1-lc[2])+
		q[2]*(1-lc[0])*(lc[1])*(1-lc[2])+q[3]*(lc[0])*(lc[1])*(1-lc[2])+
		q[4]*(1-lc[0])*(1-lc[1])*(lc[2])+q[5]*(lc[0])*(1-lc[1])*(lc[2])+
		q[6]*(1-lc[0])*(lc[1])*(lc[2])+q[7]*(lc[0])*(lc[1])*(lc[2]);
}

void GetNodeCoefficientsByLocalcCoords(double *lc,double *cff)
{
	cff[0]=(1-lc[0])*(1-lc[1])*(1-lc[2]);
	cff[1]=(lc[0])*(1-lc[1])*(1-lc[2]);
	cff[2]=(1-lc[0])*(lc[1])*(1-lc[2]);
	cff[3]=(lc[0])*(lc[1])*(1-lc[2]);
	cff[4]=(1-lc[0])*(1-lc[1])*(lc[2]);
	cff[5]=(lc[0])*(1-lc[1])*(lc[2]);
	cff[6]=(1-lc[0])*(lc[1])*(lc[2]);
	cff[7]=(lc[0])*(lc[1])*(lc[2]);
}

void GetGlobalCoordinates(PointXYZ *HexPnt,double *lc,double *gc)
{
	double ksi,eta,phi;

	ksi=lc[0];
	eta=lc[1];
	phi=lc[2];

	gc[0]=	HexPnt[0].x*(1-ksi)*(1-eta)*(1-phi)+HexPnt[1].x*(ksi)*(1-eta)*(1-phi)+
			HexPnt[2].x*(1-ksi)*(eta)*(1-phi)+HexPnt[3].x*(ksi)*(eta)*(1-phi)+
			HexPnt[4].x*(1-ksi)*(1-eta)*(phi)+HexPnt[5].x*(ksi)*(1-eta)*(phi)+
			HexPnt[6].x*(1-ksi)*(eta)*(phi)+HexPnt[7].x*(ksi)*(eta)*(phi);
	gc[1]=	HexPnt[0].y*(1-ksi)*(1-eta)*(1-phi)+HexPnt[1].y*(ksi)*(1-eta)*(1-phi)+
			HexPnt[2].y*(1-ksi)*(eta)*(1-phi)+HexPnt[3].y*(ksi)*(eta)*(1-phi)+
			HexPnt[4].y*(1-ksi)*(1-eta)*(phi)+HexPnt[5].y*(ksi)*(1-eta)*(phi)+
			HexPnt[6].y*(1-ksi)*(eta)*(phi)+HexPnt[7].y*(ksi)*(eta)*(phi);
	gc[2]=	HexPnt[0].z*(1-ksi)*(1-eta)*(1-phi)+HexPnt[1].z*(ksi)*(1-eta)*(1-phi)+
			HexPnt[2].z*(1-ksi)*(eta)*(1-phi)+HexPnt[3].z*(ksi)*(eta)*(1-phi)+
			HexPnt[4].z*(1-ksi)*(1-eta)*(phi)+HexPnt[5].z*(ksi)*(1-eta)*(phi)+
			HexPnt[6].z*(1-ksi)*(eta)*(phi)+HexPnt[7].z*(ksi)*(eta)*(phi);
}

void GetSolutionOnHex(PointXYZ &R,PointXYZ *CntHex,const double *q,double &Val)
{
	double lc[3];
	FindLocalCoordinates(R,CntHex,lc);
	GetNodeSolutionByLocalcCoords(lc,q,Val);
}

void FindLocalCoordinates(PointXYZ &R,PointXYZ *HexPnt,double *lc)
{
	double EpsForFindLocalCoord = 1e-3;
	int MaxDeep = 20;


	int p,t,m,deep,crd[3][2]={{1,2},{0,2},{0,1}};
	double LocalCoord[2][3],CentGlob[3],CentLoc[3],DopLoc[3];
	double dist,disc,h[3];
	const double ods=0.16666666666666666;
	const double hds=0.33333333333333333;
	deep=0;
	LocalCoord[0][0]=0.0;LocalCoord[0][1]=0.0;LocalCoord[0][2]=0.0;
	LocalCoord[1][0]=1.0;LocalCoord[1][1]=1.0;LocalCoord[1][2]=1.0;
	CentLoc[0]=0.5;CentLoc[1]=0.5;CentLoc[2]=0.5;
	GetGlobalCoordinates(HexPnt,CentLoc,CentGlob);
	disc=sqrt(
		(R.x-CentGlob[0])*(R.x-CentGlob[0])+
		(R.y-CentGlob[1])*(R.y-CentGlob[1])+
		(R.z-CentGlob[2])*(R.z-CentGlob[2])
		);
	h[0]=LocalCoord[1][0]-LocalCoord[0][0];
	h[1]=LocalCoord[1][1]-LocalCoord[0][1];
	h[2]=LocalCoord[1][2]-LocalCoord[0][2];
	do{
		if(disc<EpsForFindLocalCoord)break;

		for(m=0;m<3;m++){

			CentLoc[m]=LocalCoord[0][m]+ods*h[m];
			
			GetGlobalCoordinates(HexPnt,CentLoc,CentGlob);
			dist=sqrt(
					(R.x-CentGlob[0])*(R.x-CentGlob[0])+
					(R.y-CentGlob[1])*(R.y-CentGlob[1])+
					(R.z-CentGlob[2])*(R.z-CentGlob[2])
					);

			if(disc<dist){
				CentLoc[m]=LocalCoord[1][m]-ods*h[m];
			
				GetGlobalCoordinates(HexPnt,CentLoc,CentGlob);
				dist=sqrt(
						(R.x-CentGlob[0])*(R.x-CentGlob[0])+
						(R.y-CentGlob[1])*(R.y-CentGlob[1])+
						(R.z-CentGlob[2])*(R.z-CentGlob[2])
						);
				if(dist<disc){
					disc=dist;
					LocalCoord[0][m]=LocalCoord[1][m]-hds*h[m];
				}
				else{
					LocalCoord[0][m]=LocalCoord[0][m]+hds*h[m];
					LocalCoord[1][m]=LocalCoord[1][m]-hds*h[m];
				}
			}
			else{
				disc=dist;
				LocalCoord[1][m]=LocalCoord[0][m]+hds*h[m];				
			}
			h[m]=LocalCoord[1][m]-LocalCoord[0][m];
			CentLoc[m]=0.5*(LocalCoord[0][m]+LocalCoord[1][m]);
		}
		deep++;
	}while(deep<MaxDeep);
	
	if(deep==MaxDeep){
		DopLoc[0]=CentLoc[0];
		DopLoc[1]=CentLoc[1];
		DopLoc[2]=CentLoc[2];
		do{
			t=0;
			for(m=0;m<3;m++){
				do{
					p=0;
					DopLoc[m]=CentLoc[m]-h[m];				
					if(DopLoc[m]>0){
						GetGlobalCoordinates(HexPnt,DopLoc,CentGlob);
						dist=sqrt(
								(R.x-CentGlob[0])*(R.x-CentGlob[0])+
								(R.y-CentGlob[1])*(R.y-CentGlob[1])+
								(R.z-CentGlob[2])*(R.z-CentGlob[2])
								);
						if(dist<disc){
							disc=dist;
							CentLoc[m]=DopLoc[m];
							t=p=1;
						}
					}
				}while(p);
				do{
					p=0;
					DopLoc[m]=CentLoc[m]+h[m];				
					if(DopLoc[m]<1){
						GetGlobalCoordinates(HexPnt,DopLoc,CentGlob);
						dist=sqrt(
								(R.x-CentGlob[0])*(R.x-CentGlob[0])+
								(R.y-CentGlob[1])*(R.y-CentGlob[1])+
								(R.z-CentGlob[2])*(R.z-CentGlob[2])
								);
						if(dist<disc){
							disc=dist;
							CentLoc[m]=DopLoc[m];
							t=p=1;
						}
					}
				}while(p);
			}
		}while(t);
	}

	GetGlobalCoordinates(HexPnt,CentLoc,CentGlob);

	if(disc>EpsForFindLocalCoord)
	{
		logfile<<"Reciver Eps Fail:"<<endl;
		logfile<<R.x<<'\t'<<R.y<<'\t'<<R.z<<endl;
		logfile<<CentGlob[0]<<'\t'<<CentGlob[1]<<'\t'<<CentGlob[2]<<endl;
		logfile<<sqrt(R.x*R.x+R.y*R.y+R.z*R.z)<<endl;
		logfile<<sqrt((R.x-CentGlob[0])*(R.x-CentGlob[0])+
				   (R.y-CentGlob[1])*(R.y-CentGlob[1])+
				   (R.z-CentGlob[2])*(R.z-CentGlob[2]))<<endl;
		logfile<<endl;
	}
	

	lc[0]=CentLoc[0];
	lc[1]=CentLoc[1];
	lc[2]=CentLoc[2];
}

bool CheckPointInHex(PointXYZ &R,PointXYZ *HexPnt,double &mdist)
{
	bool mnf;
	double mnv;
	double EpsForFindLocalCoord = 1e-2;
	int MaxDeep = 10;

	int m,deep,crd[3][2]={{1,2},{0,2},{0,1}};
	double LocalCoord[2][3],CentGlob[3],CentLoc[3];
	double dist,disc,h[3],mind,ddg[4];
	const double ods=0.16666666666666666;
	const double hds=0.33333333333333333;
	deep=0;
	mnf=true;
	mnv=1e+30;
	LocalCoord[0][0]=0.0;LocalCoord[0][1]=0.0;LocalCoord[0][2]=0.0;
	LocalCoord[1][0]=1.0;LocalCoord[1][1]=1.0;LocalCoord[1][2]=1.0;
	CentLoc[0]=0.5;CentLoc[1]=0.5;CentLoc[2]=0.5;
	GetGlobalCoordinates(HexPnt,CentLoc,CentGlob);
	disc=sqrt(
		(R.x-CentGlob[0])*(R.x-CentGlob[0])+
		(R.y-CentGlob[1])*(R.y-CentGlob[1])+
		(R.z-CentGlob[2])*(R.z-CentGlob[2])
		);
	h[0]=LocalCoord[1][0]-LocalCoord[0][0];
	h[1]=LocalCoord[1][1]-LocalCoord[0][1];
	h[2]=LocalCoord[1][2]-LocalCoord[0][2];

	ddg[0]=sqrt((HexPnt[7].x-HexPnt[0].x)*(HexPnt[7].x-HexPnt[0].x)+(HexPnt[7].y-HexPnt[0].y)*(HexPnt[7].y-HexPnt[0].y)+(HexPnt[7].z-HexPnt[0].z)*(HexPnt[7].z-HexPnt[0].z));
	ddg[1]=sqrt((HexPnt[6].x-HexPnt[1].x)*(HexPnt[6].x-HexPnt[1].x)+(HexPnt[6].y-HexPnt[1].y)*(HexPnt[6].y-HexPnt[1].y)+(HexPnt[6].z-HexPnt[1].z)*(HexPnt[6].z-HexPnt[1].z));
	ddg[2]=sqrt((HexPnt[5].x-HexPnt[2].x)*(HexPnt[5].x-HexPnt[2].x)+(HexPnt[5].y-HexPnt[2].y)*(HexPnt[5].y-HexPnt[2].y)+(HexPnt[5].z-HexPnt[2].z)*(HexPnt[5].z-HexPnt[2].z));
	ddg[3]=sqrt((HexPnt[4].x-HexPnt[3].x)*(HexPnt[4].x-HexPnt[3].x)+(HexPnt[4].y-HexPnt[3].y)*(HexPnt[4].y-HexPnt[3].y)+(HexPnt[4].z-HexPnt[3].z)*(HexPnt[4].z-HexPnt[3].z));
	
	mind=ddg[0];

	for(m=1;m<4;m++)
	{
		if(ddg[m]<mind)mind=ddg[m];
	}

	do{
		if(disc<EpsForFindLocalCoord*mind)break;
		

		for(m=0;m<3;m++)
		{
			CentLoc[m]=LocalCoord[0][m]+ods*h[m];
			
			GetGlobalCoordinates(HexPnt,CentLoc,CentGlob);
			dist=sqrt(
					(R.x-CentGlob[0])*(R.x-CentGlob[0])+
					(R.y-CentGlob[1])*(R.y-CentGlob[1])+
					(R.z-CentGlob[2])*(R.z-CentGlob[2])
					);

			if(disc<dist)
			{
				CentLoc[m]=LocalCoord[1][m]-ods*h[m];
			
				GetGlobalCoordinates(HexPnt,CentLoc,CentGlob);
				dist=sqrt(
						(R.x-CentGlob[0])*(R.x-CentGlob[0])+
						(R.y-CentGlob[1])*(R.y-CentGlob[1])+
						(R.z-CentGlob[2])*(R.z-CentGlob[2])
						);
				if(dist<disc)
				{
					disc=dist;
					LocalCoord[0][m]=LocalCoord[1][m]-hds*h[m];
				}
				else
				{
					LocalCoord[0][m]=LocalCoord[0][m]+hds*h[m];
					LocalCoord[1][m]=LocalCoord[1][m]-hds*h[m];
				}
			}
			else
			{
				disc=dist;
				LocalCoord[1][m]=LocalCoord[0][m]+hds*h[m];				
			}
			h[m]=LocalCoord[1][m]-LocalCoord[0][m];
			CentLoc[m]=0.5*(LocalCoord[0][m]+LocalCoord[1][m]);
		}
		deep++;
	}while(deep<MaxDeep);


	mdist=dist;

	return (deep<MaxDeep);
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


void Calc_J_optimized(double *x, double *y, double *z, 
					  double (*J_1_T)[3], double *det_J_abs,  int n_of_point)
{
	int i, j;
	double J[3][3];
	double det_J;

	for(i=0; i<3; i++)
		for(j=0; j<3; j++)
			J[i][j] = 0.0;

	for(i=0; i<8; i++)
	{
		J[0][0] += x[i]*gauss_3_d_phi[n_of_point][i][0]; //  d_xi[i];
		J[0][1] += x[i]*gauss_3_d_phi[n_of_point][i][1]; //  d_eta[i];
		J[0][2] += x[i]*gauss_3_d_phi[n_of_point][i][2]; //  d_zeta[i];

		J[1][0] += y[i]*gauss_3_d_phi[n_of_point][i][0]; //  d_xi[i];
		J[1][1] += y[i]*gauss_3_d_phi[n_of_point][i][1]; //  d_eta[i];
		J[1][2] += y[i]*gauss_3_d_phi[n_of_point][i][2]; //  d_zeta[i];

		J[2][0] += z[i]*gauss_3_d_phi[n_of_point][i][0]; //  d_xi[i];
		J[2][1] += z[i]*gauss_3_d_phi[n_of_point][i][1]; //  d_eta[i];
		J[2][2] += z[i]*gauss_3_d_phi[n_of_point][i][2]; //  d_zeta[i];
	}
	det_J = J[0][0]*J[1][1]*J[2][2] - J[0][0]*J[1][2]*J[2][1] + J[1][0]*J[2][1]*J[0][2]
	- J[1][0]*J[0][1]*J[2][2] + J[2][0]*J[0][1]*J[1][2] - J[2][0]*J[1][1]*J[0][2];

	*det_J_abs = fabs(det_J);

	J_1_T[0][0] = (J[1][1]*J[2][2]-J[2][1]*J[1][2])/det_J;
	J_1_T[1][0] = (J[2][1]*J[0][2]-J[0][1]*J[2][2])/det_J;
	J_1_T[2][0] = (J[0][1]*J[1][2]-J[1][1]*J[0][2])/det_J;
	J_1_T[0][1] = (-J[1][0]*J[2][2]+J[2][0]*J[1][2])/det_J;
	J_1_T[1][1] = (J[0][0]*J[2][2]-J[2][0]*J[0][2])/det_J;
	J_1_T[2][1] = (-J[0][0]*J[1][2]+J[1][0]*J[0][2])/det_J;
	J_1_T[0][2] = (J[1][0]*J[2][1]-J[2][0]*J[1][1])/det_J;
	J_1_T[1][2] = (-J[0][0]*J[2][1]+J[2][0]*J[0][1])/det_J;
	J_1_T[2][2] = (J[0][0]*J[1][1]-J[1][0]*J[0][1])/det_J;
}


void Calc_Hex_Local_Matrix_Dx(double *x, double *y, double *z, double (*dx)[8])
{
	int i, j, i1, j1;
	double phi1;
	double gauss_3_mult;
	double J_1_T[3][3];
	double det_J_abs;
	double grad_all[8]; //   .    

	for(i=0; i<8; i++) //  
		for(j=0; j<8; j++)
		{
			dx[i][j] = 0.0;
		}

		for(i=0; i<27; i++) //    
		{
			Calc_J_optimized(x, y, z, J_1_T, &det_J_abs, i); //  
			gauss_3_mult = gauss_3_A_all[i]*det_J_abs; // A_i*A_j*A_k*|J|

			for(j=0; j<8; j++) //    -   ...
				grad_all[j]    = J_1_T[0][0]*gauss_3_d_phi[i][j][0]+  
								 J_1_T[0][1]*gauss_3_d_phi[i][j][1]+  
								 J_1_T[0][2]*gauss_3_d_phi[i][j][2];

			for(i1=0; i1<8; i1++)
			{
				phi1 = gauss_3_phi[i][i1]*gauss_3_mult;	
				for(j1=0; j1<8; j1++)
				{
					dx[i1][j1] += grad_all[j1]*phi1;
				}//j1
			}//i1
		}// i

	for(i=0; i<8; i++) // transpose
		for(j=0; j<i; j++)
			swap(dx[i][j], dx[j][i]);			
}
void Calc_Hex_Local_Matrix_Dy(double *x, double *y, double *z, double (*dy)[8])
{
	int i, j, i1, j1;
	double phi1;
	double gauss_3_mult;
	double J_1_T[3][3];
	double det_J_abs;
	double grad_all[8]; //   .    

	for(i=0; i<8; i++) //  
		for(j=0; j<8; j++)
		{
			dy[i][j] = 0.0;
		}

		for(i=0; i<27; i++) //    
		{
			Calc_J_optimized(x, y, z, J_1_T, &det_J_abs, i); //  
			gauss_3_mult = gauss_3_A_all[i]*det_J_abs; // A_i*A_j*A_k*|J|

			for(j=0; j<8; j++) //    -   ...
				grad_all[j]   =		J_1_T[1][0]*gauss_3_d_phi[i][j][0]+  
									J_1_T[1][1]*gauss_3_d_phi[i][j][1]+  
									J_1_T[1][2]*gauss_3_d_phi[i][j][2];

			for(i1=0; i1<8; i1++)
			{
				phi1 = gauss_3_phi[i][i1]*gauss_3_mult;	
				for(j1=0; j1<8; j1++)
				{
					dy[i1][j1] += grad_all[j1]*phi1;
				}//j1
			}//i1
		}// i

	for(i=0; i<8; i++) // transpose
		for(j=0; j<i; j++)
			swap(dy[i][j], dy[j][i]);
}
void Calc_Hex_Local_Matrix_Dz(double *x, double *y, double *z, double (*dz)[8])
{
	int i, j, i1, j1;
	double phi1;
	double gauss_3_mult;
	double J_1_T[3][3];
	double det_J_abs;
	double grad_all[8]; //   .    

	for(i=0; i<8; i++) //  
		for(j=0; j<8; j++)
		{
			dz[i][j] = 0.0;
		}

		for(i=0; i<27; i++) //    
		{
			Calc_J_optimized(x, y, z, J_1_T, &det_J_abs, i); //  
			gauss_3_mult = gauss_3_A_all[i]*det_J_abs; // A_i*A_j*A_k*|J|

			for(j=0; j<8; j++) //    -   ...
				grad_all[j] =		J_1_T[2][0]*gauss_3_d_phi[i][j][0]+  
									J_1_T[2][1]*gauss_3_d_phi[i][j][1]+  
									J_1_T[2][2]*gauss_3_d_phi[i][j][2];

			for(i1=0; i1<8; i1++)
			{
				phi1 = gauss_3_phi[i][i1]*gauss_3_mult;	
				for(j1=0; j1<8; j1++)
				{
					dz[i1][j1] += grad_all[j1]*phi1;
				}//j1
			}//i1
		}// i

	for(i=0; i<8; i++) // transpose
		for(j=0; j<i; j++)
			swap(dz[i][j], dz[j][i]);
}

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

/*!     */
double GetACoord(const double qc[8], 
				 const PointForResultCalculationInMesh3D& cpnt,
				 const RectForResultCalculationInMesh3D& crect)
{
	if(!cpnt.ishex)
	{
		return (qc[0]*cpnt.xb*cpnt.yb*cpnt.zb+ 
				qc[1]*cpnt.xa*cpnt.yb*cpnt.zb+
				qc[2]*cpnt.xb*cpnt.ya*cpnt.zb+
				qc[3]*cpnt.xa*cpnt.ya*cpnt.zb+
				qc[4]*cpnt.xb*cpnt.yb*cpnt.za+ 
				qc[5]*cpnt.xa*cpnt.yb*cpnt.za+
				qc[6]*cpnt.xb*cpnt.ya*cpnt.za+
				qc[7]*cpnt.xa*cpnt.ya*cpnt.za)/(crect.hx*crect.hy*crect.hz);
	}
	else
	{
		double val,lc[3];
		lc[0]=cpnt.ksi;
		lc[1]=cpnt.eta;
		lc[2]=cpnt.zet;
		GetNodeSolutionByLocalcCoords(lc,qc,val);
		return val;
	}
}

/*!       3d  */
int ProcessResCalcPoints(Mesh3DForVP* ptrMesh, const double* _v)
{
	int i, j;

	if (ptrMesh->boundsV.beginP==ptrMesh->boundsV.endP&&
		ptrMesh->boundsVforE.beginP==ptrMesh->boundsVforE.endP)
		return RETCODE_OK;

	for (i=0; i<ptrMesh->ResCalcRectNum; i++)
	{
		RectForResultCalculationInMesh3D& crect=ptrMesh->ResCalcRectArray[i];
		const Rect3D& elrect=ptrMesh->rect3d[crect.rnum];
		for (j=0; j<8; j++)
			crect.qv[j]=_v[elrect.nodes[j]-1];
	}
	for (i=ptrMesh->boundsV.beginP; i<ptrMesh->boundsV.endP; i++)
	{
		PointForResultCalculationInMesh3D& cpnt=ptrMesh->ResCalcPointArray[i];
		PointXYZ& cpntres=cpnt.res[0]; 
		const RectForResultCalculationInMesh3D& crect=ptrMesh->ResCalcRectArray[cpnt.rnum];
		cpntres.x=GetACoord(crect.qv, cpnt, crect);
	}

	return RETCODE_OK;
}
int MatrixOnVector(
	const double *x1,	// z=A(x1+x2)// IN 
	const double *x2,	// 			 //
	const double *adiag,// matrix	 //
	const double *altr,	//			 //
	const double *autr,	//			 //
	const int *iptr,	// portrait	 //
	const int *jptr,	//			 //
	const int &n,		// size		 //
	double *z)			//			 // OUT
{
	int i, j, jj;
	for (i=0; i<n; i++) 
		z[i]=(x1[i]+x2[i])*adiag[i];
	for (i=0; i<n; i++)  
		for (j=iptr[i]-1; j<iptr[i+1]-1; j++)
		{
			jj=jptr[j]-1;
			z[i]+=(x1[jj]+x2[jj])*altr[j];
			z[jj]+=(x1[i]+x2[i])*autr[j];
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

int GetLocalContributionsForRect3DInitMatixes(Mesh3DForVP* ptrMesh, SparseSLAE* ptrSLAE,
								   void (*GradF)(double dc[8][8], double hx, double hy, double hz),
								   void (*GradFHex)(double *x, double *y, double *z, double (*_d)[8]),
								   int _mtr, int coord)
{
	int i, j, k, _i, _j, _ti, c, d;
	double lm[8][8];
	double x_c[8], y_c[8], z_c[8];
	vector<LocalMatrix3D> *pvG;

	pvG=(!coord)? &(ptrMesh->vGx) : (coord==1)? &(ptrMesh->vGy) : &(ptrMesh->vGz);

	for (k=0; k<ptrMesh->krect3d; k++)
	{
		const Rect3D& r=ptrMesh->rect3d[k];
		const PointXYZ& p0=ptrMesh->pnt[r.nodes[0]-1];
		const PointXYZ& p7=ptrMesh->pnt[r.nodes[7]-1];
		if (_mtr!=-1&&_mtr!=r.mtr)
			continue;

		GradF(lm, p7.x-p0.x, p7.y-p0.y, p7.z-p0.z);

		for(i=0;i<8;i++)
		{
			for(j=0;j<8;j++)
			{
				(*pvG)[k].m[i][j]=lm[i][j];
			}
		}
	}

	return RETCODE_OK;
}

int GetLocalContributionsForRect3DWithMatrixes(const Mesh3DForVP* ptrMesh, SparseSLAE* ptrSLAE, 
											   bool sig, int _mtr, int coord)
{
	int i, j, k, _i, _j, _ti, c, d;
	double lm[8][8];
	double x_c[8], y_c[8], z_c[8];
	double lmij, el_sig, tsig, tdsig, vsig[3];

	for (k=0; k<ptrMesh->krect3d; k++)
	{
		const Rect3D& r=ptrMesh->rect3d[k];
		const PointXYZ& p0=ptrMesh->pnt[r.nodes[0]-1];
		const PointXYZ& p7=ptrMesh->pnt[r.nodes[7]-1];
		if (_mtr!=-1&&_mtr!=r.mtr)
			continue;

		for (i=0; i<8; i++)
		{
			const PointXYZ& phex=ptrMesh->pnt[r.nodes[i]-1];
			x_c[i]=phex.x;
			y_c[i]=phex.y;
			z_c[i]=phex.z;
		}

		tsig=ptrMesh->sig[r.mtr];
		tdsig=ptrMesh->sig[r.mtr]-ptrMesh->sig0[r.mtr];

		for (j=0; j<3; j++)
		{
			if(coord==j)
				vsig[j]=(sig) ? tsig : tdsig;
			else
				vsig[j]=0.0;
		}

		const LocalMatrix3D &lmx=ptrMesh->vGx[k];
		const LocalMatrix3D &lmy=ptrMesh->vGy[k];
		const LocalMatrix3D &lmz=ptrMesh->vGz[k];

		el_sig=1.0;

		for (i=0; i<8; i++)
		{
			for (j=0; j<8; j++)
			{
				lm[i][j]=vsig[0]*lmx.m[i][j]+vsig[1]*lmy.m[i][j]+vsig[2]*lmz.m[i][j];
			}
		}

		for (i=0; i<8; i++)
			for (j=0; j<8; j++)
			{
				_i=r.nodes[i];
				_j=r.nodes[j];
				lmij=lm[j][i];
				ptrSLAE->AddToMatrix(_i-1, _j-1, el_sig*lmij, ptrSLAE->C);
			}
	}
	return RETCODE_OK;
}

/*!    3d   */
int GetLocalContributionsForRect3D(const Mesh3DForVP* ptrMesh, SparseSLAE* ptrSLAE, bool withmu=false)
{
	int i, j, k, _i, _j, _ti, _tj, c, d;
	double lm[8][8], lm0[8][8];
	double x_c[8], y_c[8], z_c[8];
	double lmij, lmij0, el_sig, el_dsig, lsig[3], lsig0[3], tsig[3][3], tsig0[3][3], tdsig[3][3];

	for (k=0; k<ptrMesh->krect3d; k++)
	{
		const Rect3D& r=ptrMesh->rect3d[k];
		const PointXYZ& p0=ptrMesh->pnt[r.nodes[0]-1];
		const PointXYZ& p7=ptrMesh->pnt[r.nodes[7]-1];
		el_sig=ptrMesh->sig[r.mtr];
		el_dsig=ptrMesh->sig0[r.mtr]-el_sig;
		if (withmu)
		{
			el_sig=1/(MU0);
		}
		{
			for (i=0; i<8; i++)
			{
				const PointXYZ& phex=ptrMesh->pnt[r.nodes[i]-1];
				x_c[i]=phex.x;
				y_c[i]=phex.y;
				z_c[i]=phex.z;
			}
			GetStiffnessMatrixForRect3D(p0.x, p0.y, p0.z, p7.x, p7.y, p7.z, lm);
			for (i=0; i<8; i++)
			{
				for (j=0; j<8; j++)
				{
					lm0[i][j]=lm[i][j];
				}
			}
		}

		for (i=0; i<8; i++)
		{
			for (j=0; j<8; j++)
			{
				_i=r.nodes[i];
				_j=r.nodes[j];
				lmij=lm[i][j];
				lmij0=lm0[i][j];
				ptrSLAE->AddToMatrix(_i-1, _j-1, el_sig*lmij, ptrSLAE->B);
				ptrSLAE->AddToMatrix(_i-1, _j-1, el_dsig*lmij0, ptrSLAE->C);
			}
	}
	}
	return RETCODE_OK;
}

/*!    3d   */
int ProcessBoundaryConditionsForRect3DForVa(const Mesh3DForVP* ptrMesh, SparseSLAE* ptrSLAE, int npls)
{
	int ipls;
	for (int i=0; i<ptrMesh->nc; i++)
	{
		const PointXYZ& p=ptrMesh->pnt[i];
		if(fabs(ptrSLAE->A.adiag[i])<1e-30)
		{
			ptrSLAE->SetMatrixDiagElement(i, 1e30, ptrSLAE->A);
			for(ipls=0;ipls<npls;ipls++)
			{
				ptrSLAE->SetVectorElement(i+ptrSLAE->n*ipls, 0, ptrSLAE->F);
			}
		}
	}
	return RETCODE_OK;
}

/*!    3d   */
int ProcessBoundaryConditionsForRect3DMatrix(const Mesh3DForVP* ptrMesh, SparseSLAE* ptrSLAE)
{
	for (int i=0; i<ptrMesh->kt1; i++)
	{
		const int& l1i=ptrMesh->l1[i];
		if (l1i>ptrMesh->nc) continue;
		const PointXYZ& p=ptrMesh->pnt[l1i-1];
		ptrSLAE->SetMatrixDiagElement(l1i-1, 1e30, ptrSLAE->A);
	}
	return RETCODE_OK;
}

/*!    3d   */
int ProcessBoundaryConditionsForRect3D(const Mesh3DForVP* ptrMesh, SparseSLAE* ptrSLAE,int npls)
{
	int ipls;
	for (int i=0; i<ptrMesh->kt1; i++)
	{
		const int& l1i=ptrMesh->l1[i];
		if (l1i>ptrMesh->nc) continue;
		const PointXYZ& p=ptrMesh->pnt[l1i-1];
		for(ipls=0;ipls<npls;ipls++)
		{
			ptrSLAE->SetVectorElement(l1i-1+ptrSLAE->n*ipls, 0, ptrSLAE->F);
		}
	}
	return RETCODE_OK;
}

/*!    3d   */
int ProcessBoundaryConditionsForRect3DVector(const Mesh3DForVP* ptrMesh, SparseSLAE* ptrSLAE)
{
	for (int i=0; i<ptrMesh->kt1; i++)
	{
		const int& l1i=ptrMesh->l1[i];
		if (l1i>ptrMesh->nc) continue;
		const PointXYZ& p=ptrMesh->pnt[l1i-1];
		ptrSLAE->SetVectorElement(l1i-1, 0, ptrSLAE->F);
	}
	return RETCODE_OK;
}

/*!   3d   */
int BuildPortrait(const Mesh3DForVP* ptrMesh, SparseSLAE* ptrSLAE,int npls)
{
	int i, j, k, retc, node, nextnode;
	PortraitAL P(ptrMesh->nc);
	ResizableArrayOf<int> indlist(64);
	for (i=0; i<ptrMesh->krect3d; i++)
	{
		const Rect3D& r=ptrMesh->rect3d[i];
		indlist.SetSize(0);
		for (j=0; j<8; j++)
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
	if ((retc=ptrSLAE->B.Init(ptrSLAE->n, ptrSLAE->jsize))!=0) // Bij=sigma*gradFIi*gradFIj
		return retc;
	if ((retc=ptrSLAE->C.Init(ptrSLAE->n, ptrSLAE->jsize))!=0) // Cij=(sigma0-sigma)*gradFIi*gradFIj
		return retc;
	if ((retc=ptrSLAE->F.Init(ptrSLAE->n*npls))!=0)
		return retc;
	if ((retc=ptrSLAE->U.Init(ptrSLAE->n*npls))!=0)
		return retc;
	if ((retc=ptrSLAE->tmpV.Init(ptrSLAE->n*npls))!=0)
		return retc;
	if ((retc=ptrSLAE->F3.Init(ptrSLAE->n*npls*3))!=0)
		return retc;

	return RETCODE_OK;
}

/*!       3d  */
int GetLocalStiffnessContributionsForRect3D(const Mesh3DForVP* ptrMesh, SparseSLAE* ptrSLAE)
{
	int i, j, k, _i, _j, _ti, c, d;
	double lm[8][8];
	double x_c[8], y_c[8], z_c[8];
	double lmij, el_mu=1./(MU0);

	for (k=0; k<ptrMesh->krect3d; k++)
	{
		const Rect3D& r=ptrMesh->rect3d[k];
		const PointXYZ& p0=ptrMesh->pnt[r.nodes[0]-1];
		const PointXYZ& p7=ptrMesh->pnt[r.nodes[7]-1];
		GetStiffnessMatrixForRect3D(p0.x, p0.y, p0.z, p7.x, p7.y, p7.z, lm);
		for (i=0; i<8; i++)
			for (j=0; j<8; j++)
			{
				_i=r.nodes[i];
				_j=r.nodes[j];
				lmij=lm[i][j];
				ptrSLAE->AddToMatrix(_i-1, _j-1, el_mu*lmij, ptrSLAE->B);
			}
	}
	return RETCODE_OK;
}

int GetNumberOfPlaces(int &npls)
{
	int i,k,p1,p2,nprof;
	ifstream inf;

	npls=0;
	inf.open("group");
	if(!inf)
	{
		logfile<<"Error in open file "<<"group"<<endl;
		cout<<"Error in open file "<<"group"<<endl;
		return 1;
	}
	inf>>nprof;
	for(i=0;i<nprof;i++)
	{
		inf>>k>>p1>>p2;
		npls+=p2-p1+1;
	}
	inf.close();
	inf.clear();

	return 0;
}

int Mesh3DForVP::GetNearestElement(double *gc,double *lc,int fuseair)
{
	double scurr,sbest,m_s_best[8];
	int i,elem,ib,jb,i_best,n_best,m_i_best[8];
	PointXYZ Pmin,Pmax,Pcur,Hex[8],p,loc;
	if(!krect3d){return -1;}
	p.x=gc[0];
	p.y=gc[1];
	p.z=gc[2];
	n_best=8;
	for(ib=0;ib<n_best;ib++)
	{
		m_s_best[ib]=1e+30;
		m_i_best[ib]=-1;
	}
	sbest=m_s_best[n_best-1];
	if(krect3d<n_best){n_best=krect3d;}

	for(elem=0;elem<krect3d;elem++)
	{
		if(fuseair || rect3d[elem].mtr!=1)	//   
		{
			Pcur=pnt[rect3d[elem].nodes[0]-1];
			Hex[0]=Pmax=Pmin=Pcur;
			for(i=1;i<8;i++)
			{
				Pcur=pnt[rect3d[elem].nodes[i]-1];
				Hex[i]=Pcur;
				if(Pcur.x<Pmin.x){Pmin.x=Pcur.x;}
				if(Pcur.y<Pmin.y){Pmin.y=Pcur.y;}
				if(Pcur.z<Pmin.z){Pmin.z=Pcur.z;}
				if(Pcur.x>Pmax.x){Pmax.x=Pcur.x;}
				if(Pcur.y>Pmax.y){Pmax.y=Pcur.y;}
				if(Pcur.z>Pmax.z){Pmax.z=Pcur.z;}
			}

			scurr=0.0;
			scurr += (p.x<Pmin.x)? (Pmin.x-p.x) : (p.x>Pmax.x)? (p.x-Pmax.x) : 0.0;
			scurr += (p.y<Pmin.y)? (Pmin.y-p.y) : (p.y>Pmax.y)? (p.y-Pmax.y) : 0.0;
			scurr += (p.z<Pmin.z)? (Pmin.z-p.z) : (p.z>Pmax.z)? (p.z-Pmax.z) : 0.0;

			if(scurr<sbest)
			{
				for(ib=0;ib<n_best;ib++)
				{
					if(m_i_best[ib]==-1 || scurr<m_s_best[ib])
					{
						break;
					}
				}
					if(m_i_best[ib]==-1)
					{
						m_i_best[ib]=elem;
						m_s_best[ib]=scurr;
					}
					else
					{
						for(jb=n_best-1;jb>ib;jb--)
						{
							m_i_best[jb]=m_i_best[jb-1];
							m_s_best[jb]=m_s_best[jb-1];
						}
						m_i_best[ib]=elem;
						m_s_best[ib]=scurr;
					}
				sbest=m_s_best[n_best-1];
			}
		}
	}

	i_best=-1;
	sbest=1e+30;
	for(ib=0;ib<n_best;ib++)
	{
		elem=m_i_best[ib];
		if(elem==-1){break;}
		for(i=0;i<8;i++){Hex[i]=pnt[rect3d[elem].nodes[i]-1];}
		if(!CheckInHex(Hex,p,loc))
		{
			scurr=0.0;
			scurr += (loc.x<0.0)? (-loc.x) : (loc.x>1.0)? (loc.x-1.0) : 0.0;
			scurr += (loc.y<0.0)? (-loc.y) : (loc.y>1.0)? (loc.y-1.0) : 0.0;
			scurr += (loc.z<0.0)? (-loc.z) : (loc.z>1.0)? (loc.z-1.0) : 0.0;
			if(scurr<sbest)
			{
				sbest=scurr;
				i_best=elem;
				lc[0]=loc.x;
				lc[1]=loc.y;
				lc[2]=loc.z;
			}
		}
	}

	return i_best;
}

// Calculation of a stationary electric field
int SolveStationarProblem3d(const char* path,bool f3d,bool forA0)
{
	ifstream inf;
	char buf[256];
	int  i,j,k,retc,ipls,npls;
	FILE *fp;
	int *ig,*jg;
	pardiso_solver prds;
	ofstream outfMN,outfMOON;
	ifstream infmn;
	ofstream outf;
	vector<int> vNrec;

	Mesh3DForVP meshvp3d(moUnloadSolutions);

	SparseSLAE SLAE3D;

	retc=GetNumberOfPlaces(npls);
	if(retc)
	{
		cout << "Function GetNumberOfPlaces retrurned " << retc << '\n';
		logfile << "Function GetNumberOfPlaces retrurned " << retc << '\n';
		return 1;
	}


	meshvp3d.npls=npls;

	if ((retc=meshvp3d.Read(path))!=0)
		return retc;

	if ((retc=BuildPortrait(&meshvp3d, &SLAE3D, npls))!=0)
		return retc;

	ig=new int[SLAE3D.n+1];
	jg=new int[SLAE3D.jsize];

	for(i=0;i<(SLAE3D.n+1);i++){ig[i]=SLAE3D.iptr[i]-1;}
	for(i=0;i<SLAE3D.jsize;i++){jg[i]=SLAE3D.jptr[i]-1;}

	if ((retc=GetLocalContributionsForRect3D(&meshvp3d, &SLAE3D))!=0)
		return retc;
	
	SLAE3D.A.Clear(SLAE3D.n,SLAE3D.jsize);
	SLAE3D.F.Clear(SLAE3D.n*npls);

	SLAE3D.A.AddToMatrix(SLAE3D.n, SLAE3D.jsize, SLAE3D.B);

	for(i=0;i<meshvp3d.kpnt*npls;i++)
	{
		meshvp3d.qva[i]=meshvp3d.qvb[i]=0.0;
	}

#ifndef FDLINE
	{
		fp=fopen("va.dat","rb");
		fread(meshvp3d.qva,sizeof(double),meshvp3d.kpnt*npls,fp);
		fclose(fp);

		fp=fopen("vb.dat","rb");
		fread(meshvp3d.qvb,sizeof(double),meshvp3d.kpnt*npls,fp);
		fclose(fp);
	}
#endif

	ofstream ofp;

	ofp.open("v3n.dat",ios::binary);
	for(ipls=0;ipls<npls;ipls++)
	{
		for(i=0;i<meshvp3d.kpnt;i++)
		{
			ofp<(meshvp3d.qva[i+meshvp3d.kpnt*ipls]+meshvp3d.qvb[i+meshvp3d.kpnt*ipls]);
		}
	}
	ofp.close();
	ofp.clear();


	for(ipls=0;ipls<npls;ipls++)
	{
		MatrixOnVector(meshvp3d.qva+meshvp3d.kpnt*ipls, meshvp3d.qvb+meshvp3d.kpnt*ipls, 
			SLAE3D.C.adiag, SLAE3D.C.altr, SLAE3D.C.autr, 
			SLAE3D.iptr, SLAE3D.jptr, SLAE3D.n, SLAE3D.F.v+meshvp3d.nc*ipls);
	}


	for(i=0;i<meshvp3d.nc;i++)
	{
		if(!(SLAE3D.A.adiag[i]))
		{
			SLAE3D.A.adiag[i]=1.0;
			SLAE3D.F.v[i]=0.0;
		}
	}

	if ((retc=ProcessBoundaryConditionsForRect3DMatrix(&meshvp3d, &SLAE3D))!=0)
		return retc;

	if ((retc=ProcessBoundaryConditionsForRect3D(&meshvp3d, &SLAE3D, npls))!=0)
		return retc;





	prds.factorize(SLAE3D.n,ig,jg,SLAE3D.A.altr,SLAE3D.A.adiag,1);
	prds.solve_nrhs(npls,SLAE3D.F.v,SLAE3D.pU->v);
	prds.stop_solver();


	if(CheckStop())return 1;

	ofp.open("v3a.dat",ios::binary);
	for(ipls=0;ipls<npls;ipls++)
	{
		if ((retc=meshvp3d.GetWeights(SLAE3D.pU->v+ipls*SLAE3D.n, false))!=0)
			return retc;

		if ((retc=ProcessResCalcPoints(&meshvp3d, meshvp3d.q_all))!=0)
			return retc;

		ofp.write((char *)meshvp3d.q_all,meshvp3d.kpnt*sizeof(double));
	


	}
	ofp.close();
	ofp.clear();

	if(forA0)
	{
		meshvp3d.vGx.resize(meshvp3d.krect3d);
		meshvp3d.vGy.resize(meshvp3d.krect3d);
		meshvp3d.vGz.resize(meshvp3d.krect3d);

		if ((retc=GetLocalContributionsForRect3DInitMatixes(&meshvp3d, &SLAE3D, GradMatrix_X, Calc_Hex_Local_Matrix_Dx,-1,0))!=0)
			return retc;

		if ((retc=GetLocalContributionsForRect3DInitMatixes(&meshvp3d, &SLAE3D, GradMatrix_Y, Calc_Hex_Local_Matrix_Dy,-1,1))!=0)
			return retc;

		if ((retc=GetLocalContributionsForRect3DInitMatixes(&meshvp3d, &SLAE3D, GradMatrix_Z, Calc_Hex_Local_Matrix_Dz,-1,2))!=0)
			return retc;

		double *ax,*ay,*az;

		if ((ax=new double[meshvp3d.kpnt])==NULL)
			return RETCODE_NOMEM;
		if ((ay=new double[meshvp3d.kpnt])==NULL)
			return RETCODE_NOMEM;
		if ((az=new double[meshvp3d.kpnt])==NULL)
			return RETCODE_NOMEM;

		SLAE3D.A.Clear(SLAE3D.n, SLAE3D.jsize);
		SLAE3D.B.Clear(SLAE3D.n, SLAE3D.jsize);
		if ((retc=GetLocalStiffnessContributionsForRect3D(&meshvp3d, &SLAE3D))!=0)
			return retc;
		SLAE3D.A.AddToMatrix(SLAE3D.n, SLAE3D.jsize, SLAE3D.B);

		if ((retc=ProcessBoundaryConditionsForRect3DMatrix(&meshvp3d, &SLAE3D))!=0)
			return retc;

		prds.factorize(SLAE3D.n,ig,jg,SLAE3D.A.altr,SLAE3D.A.adiag,1);


		double *u3=NULL;

		u3=new double[3*meshvp3d.nc*npls];

		if(CheckStop())return 1;

		for(ipls=0;ipls<npls;ipls++)
		{
			for(i=0;i<SLAE3D.n;i++){meshvp3d.qv[i]=SLAE3D.pU->v[i+ipls*SLAE3D.n];}

		SLAE3D.F.Clear(SLAE3D.n);
		SLAE3D.C.Clear(SLAE3D.n, SLAE3D.jsize);
		if ((retc=GetLocalContributionsForRect3DWithMatrixes(&meshvp3d,&SLAE3D,true,-1,0))!=0)
			return retc;
		MatrixOnVector(
	#ifndef DLINE
	meshvp3d.qv,
	#else
	meshvp3d.qv_sum,
	#endif
						SLAE3D.C.adiag, SLAE3D.C.altr, SLAE3D.C.autr, 
						SLAE3D.iptr, SLAE3D.jptr, SLAE3D.n, SLAE3D.tmpV.v);
		SLAE3D.F.AddToVector(SLAE3D.n, SLAE3D.tmpV, -1);
		
	#ifndef DLINE	
		SLAE3D.C.Clear(SLAE3D.n, SLAE3D.jsize);
		if ((retc=GetLocalContributionsForRect3DWithMatrixes(&meshvp3d,&SLAE3D,false,-1,0))!=0)
			return retc;
		MatrixOnVector(meshvp3d.qva+meshvp3d.kpnt*ipls, meshvp3d.qvb+meshvp3d.kpnt*ipls, 
						SLAE3D.C.adiag, SLAE3D.C.altr, SLAE3D.C.autr, 
						SLAE3D.iptr, SLAE3D.jptr, SLAE3D.n, SLAE3D.tmpV.v);
		SLAE3D.F.AddToVector(SLAE3D.n, SLAE3D.tmpV, -1);
	#endif

		if ((retc=ProcessBoundaryConditionsForRect3DVector(&meshvp3d, &SLAE3D))!=0)
			return retc;
		for(i=0;i<SLAE3D.n;i++){SLAE3D.F3.v[i+ipls*SLAE3D.n]=SLAE3D.F.v[i];}

		}
		if(CheckStop())return 1;

		for(ipls=0;ipls<npls;ipls++)
		{
			for(i=0;i<SLAE3D.n;i++){meshvp3d.qv[i]=SLAE3D.pU->v[i+ipls*SLAE3D.n];}

		SLAE3D.F.Clear(SLAE3D.n);
		SLAE3D.C.Clear(SLAE3D.n, SLAE3D.jsize);
		if ((retc=GetLocalContributionsForRect3DWithMatrixes(&meshvp3d,&SLAE3D,true,-1,1))!=0)
			return retc;
		MatrixOnVector(
	#ifndef DLINE
	meshvp3d.qv,
	#else
	meshvp3d.qv_sum,
	#endif
						SLAE3D.C.adiag, SLAE3D.C.altr, SLAE3D.C.autr, 
						SLAE3D.iptr, SLAE3D.jptr, SLAE3D.n, SLAE3D.tmpV.v);
		SLAE3D.F.AddToVector(SLAE3D.n, SLAE3D.tmpV, -1);

	#ifndef DLINE
		SLAE3D.C.Clear(SLAE3D.n, SLAE3D.jsize);
		if ((retc=GetLocalContributionsForRect3DWithMatrixes(&meshvp3d,&SLAE3D,false,-1,1))!=0)
			return retc;
		MatrixOnVector(meshvp3d.qva+meshvp3d.kpnt*ipls, meshvp3d.qvb+meshvp3d.kpnt*ipls, 
						SLAE3D.C.adiag, SLAE3D.C.altr, SLAE3D.C.autr, 
						SLAE3D.iptr, SLAE3D.jptr, SLAE3D.n, SLAE3D.tmpV.v);
		SLAE3D.F.AddToVector(SLAE3D.n, SLAE3D.tmpV, -1);
	#endif

		if ((retc=ProcessBoundaryConditionsForRect3DVector(&meshvp3d, &SLAE3D))!=0)
			return retc;
		for(i=0;i<SLAE3D.n;i++){SLAE3D.F3.v[i+ipls*SLAE3D.n+npls*SLAE3D.n]=SLAE3D.F.v[i];}

		}
		if(CheckStop())return 1;

		for(ipls=0;ipls<npls;ipls++)
		{
			for(i=0;i<SLAE3D.n;i++){meshvp3d.qv[i]=SLAE3D.pU->v[i+ipls*SLAE3D.n];}

		SLAE3D.F.Clear(SLAE3D.n);
		SLAE3D.C.Clear(SLAE3D.n, SLAE3D.jsize);
		if ((retc=GetLocalContributionsForRect3DWithMatrixes(&meshvp3d,&SLAE3D,true,-1,2))!=0)
			return retc;
		MatrixOnVector(
	#ifndef DLINE
	meshvp3d.qv,
	#else
	meshvp3d.qv_sum,
	#endif
						SLAE3D.C.adiag, SLAE3D.C.altr, SLAE3D.C.autr, 
						SLAE3D.iptr, SLAE3D.jptr, SLAE3D.n, SLAE3D.tmpV.v);
		SLAE3D.F.AddToVector(SLAE3D.n, SLAE3D.tmpV, -1);

	#ifndef DLINE
		SLAE3D.C.Clear(SLAE3D.n, SLAE3D.jsize);
		if ((retc=GetLocalContributionsForRect3DWithMatrixes(&meshvp3d,&SLAE3D,false,-1,2))!=0)
			return retc;
		MatrixOnVector(meshvp3d.qva+meshvp3d.kpnt*ipls, meshvp3d.qvb+meshvp3d.kpnt*ipls, 
						SLAE3D.C.adiag, SLAE3D.C.altr, SLAE3D.C.autr, 
						SLAE3D.iptr, SLAE3D.jptr, SLAE3D.n, SLAE3D.tmpV.v);
		SLAE3D.F.AddToVector(SLAE3D.n, SLAE3D.tmpV, -1);
	#endif

		if ((retc=ProcessBoundaryConditionsForRect3DVector(&meshvp3d, &SLAE3D))!=0)
			return retc;
		for(i=0;i<SLAE3D.n;i++){SLAE3D.F3.v[i+ipls*SLAE3D.n+2*npls*SLAE3D.n]=SLAE3D.F.v[i];}
		
		if(ipls==npls-1)
		{
			prds.solve_nrhs(3*npls,SLAE3D.F3.v,u3);
			prds.stop_solver();
		}

		}

		for(ipls=0;ipls<npls;ipls++)
		{
			for(i=0;i<SLAE3D.n;i++){SLAE3D.pU->v[i]=u3[i+ipls*SLAE3D.n];}
			if ((retc=meshvp3d.GetWeights(SLAE3D.pU->v))!=0)
				return retc;
			CopyV(meshvp3d.kpnt, meshvp3d.q_all, ax);

			for(i=0;i<SLAE3D.n;i++){SLAE3D.pU->v[i]=u3[i+ipls*SLAE3D.n+npls*SLAE3D.n];}
			if ((retc=meshvp3d.GetWeights(SLAE3D.pU->v))!=0)
				return retc;
			CopyV(meshvp3d.kpnt, meshvp3d.q_all, ay);

			for(i=0;i<SLAE3D.n;i++){SLAE3D.pU->v[i]=u3[i+ipls*SLAE3D.n+2*npls*SLAE3D.n];}
			if ((retc=meshvp3d.GetWeights(SLAE3D.pU->v))!=0)
				return retc;
			CopyV(meshvp3d.kpnt, meshvp3d.q_all, az);

			CHECKSTRING0(buf, strlen(path)+32);
			sprintf_s(buf, "%s\\a0.dat", path);
			if(!ipls)
				outf.open(buf, ios::binary);
			else
				outf.open(buf, ios::binary|ios::app);
			for (i=0; i<meshvp3d.kpnt; i++)
				outf < ax[i] < ay[i] < az[i];
			outf.close();
			outf.clear();
		}
		
		if(u3){delete [] u3;u3=NULL;}

		if(CheckStop())return 1;

		delete [] ax;
		delete [] ay;
		delete [] az;
	}

	return RETCODE_OK;
}
int CalcSP()
{
	ifstream inf;
	bool forA0;

	forA0=false;
	inf.open("forA0");
	if(inf)
	{
		forA0=true;
		inf.close();
	}
	inf.clear();

	if(!CheckStop())
	{
		SolveStationarProblem3d(".", true, forA0);
	}

	return 0;
}
