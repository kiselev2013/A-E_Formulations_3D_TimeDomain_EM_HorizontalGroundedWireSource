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
 *  This file contains the code for calculation of field V for stationary task
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

#include "OutputResultant3d.h"

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

void normalize(PointXYZ &vec)
{
	double len=sqrt(vec.x*vec.x+vec.y*vec.y+vec.z*vec.z);
	vec.x/=len;
	vec.y/=len;
	vec.z/=len;
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

/*!      */
PointXYZ GetGradACoord(const double qc[8], 
			  		   const PointForResultCalculationInMesh3D& cpnt,
				       const RectForResultCalculationInMesh3D& crect)
{
	PointXYZ GradA;
	GradA.x=(qc[0]*(-1)*cpnt.yb*cpnt.zb+ 
			 qc[1]*cpnt.yb*cpnt.zb+
			 qc[2]*(-1)*cpnt.ya*cpnt.zb+
			 qc[3]*cpnt.ya*cpnt.zb+
			 qc[4]*(-1)*cpnt.yb*cpnt.za+ 
			 qc[5]*cpnt.yb*cpnt.za+
			 qc[6]*(-1)*cpnt.ya*cpnt.za+
			 qc[7]*cpnt.ya*cpnt.za)/(crect.hx*crect.hy*crect.hz);
	GradA.y=(qc[0]*cpnt.xb*(-1)*cpnt.zb+ 
 			 qc[1]*cpnt.xa*(-1)*cpnt.zb+
			 qc[2]*cpnt.xb*cpnt.zb+
			 qc[3]*cpnt.xa*cpnt.zb+
			 qc[4]*cpnt.xb*(-1)*cpnt.za+ 
			 qc[5]*cpnt.xa*(-1)*cpnt.za+
			 qc[6]*cpnt.xb*cpnt.za+
			 qc[7]*cpnt.xa*cpnt.za)/(crect.hx*crect.hy*crect.hz);
	GradA.z=(qc[0]*cpnt.xb*cpnt.yb*(-1)+ 
			 qc[1]*cpnt.xa*cpnt.yb*(-1)+
			 qc[2]*cpnt.xb*cpnt.ya*(-1)+
			 qc[3]*cpnt.xa*cpnt.ya*(-1)+
			 qc[4]*cpnt.xb*cpnt.yb+ 
			 qc[5]*cpnt.xa*cpnt.yb+
			 qc[6]*cpnt.xb*cpnt.ya+
			 qc[7]*cpnt.xa*cpnt.ya)/(crect.hx*crect.hy*crect.hz);
	return GradA;
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

	for (i=ptrMesh->boundsVforE.beginP; i<ptrMesh->boundsVforE.endP; i++)
	{
		PointForResultCalculationInMesh3D& cpnt=ptrMesh->ResCalcPointArray[i];
		PointXYZ& cpntres=cpnt.res[0]; 
		const RectForResultCalculationInMesh3D& crect=ptrMesh->ResCalcRectArray[cpnt.rnum];
		cpntres=GetGradACoord(crect.qv, cpnt, crect)*(-1.);
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

/*!    3d   */
int GetLocalContributionsForRect3D(const Mesh3DForVP* ptrMesh, SparseSLAE* ptrSLAE,
								   int m_0, //   
								   int m_size, //   
								   int* m_array, //  ,     
								   bool sig,	// true - sig*..., false - (sig0-sig)*..., 
								   bool sig3d)	// true - sig*..., false - sig0*..., 
{
	int i, j, k, _i, _j, _ti, c, d;
	double lm[8][8];
	double x_c[8], y_c[8], z_c[8];
	double lmij, scoeff;
	bool mtrfound;

	for (k=0; k<ptrMesh->krect3d; k++)
	{
		const Rect3D& r=ptrMesh->rect3d[k];
		const PointXYZ& p0=ptrMesh->pnt[r.nodes[0]-1];
		const PointXYZ& p7=ptrMesh->pnt[r.nodes[7]-1];
		mtrfound=(m_size==0);
		for (i=m_0; (i<m_size)&&(!mtrfound); i++)
			mtrfound=(r.mtr==m_array[i]);
		if (!mtrfound) continue;
		GetStiffnessMatrixForRect3D(p0.x, p0.y, p0.z, p7.x, p7.y, p7.z, lm);
		if (sig)
			scoeff=(sig3d?ptrMesh->sig[r.mtr]:ptrMesh->sig0[r.mtr]);
		else
			scoeff=ptrMesh->sig0[r.mtr]-ptrMesh->sig[r.mtr];
		for (i=0; i<8; i++)
			for (j=0; j<8; j++)
			{
				_i=r.nodes[i];
				_j=r.nodes[j];
				lmij=lm[i][j];
				ptrSLAE->AddToMatrix(_i-1, _j-1, scoeff*lmij, ptrSLAE->C);
			}
	}
	return RETCODE_OK;
}

/*!    3d   */
int GetLocalContributionsForRect3D(const Mesh3DForVP* ptrMesh, SparseSLAE* ptrSLAE, bool sig,
								   void (*GradF)(double dc[8][8], double hx, double hy, double hz),
								   void (*GradFHex)(double *x, double *y, double *z, double (*_d)[8]),
								   int _mtr, int coord)
{
	int i, j, k, _i, _j, _ti, c, d;
	double lm[8][8];
	double x_c[8], y_c[8], z_c[8];
	double lmij, el_sig;

	for (k=0; k<ptrMesh->krect3d; k++)
	{
		const Rect3D& r=ptrMesh->rect3d[k];
		const PointXYZ& p0=ptrMesh->pnt[r.nodes[0]-1];
		const PointXYZ& p7=ptrMesh->pnt[r.nodes[7]-1];
		if (_mtr!=-1&&_mtr!=r.mtr)
			continue;
		GradF(lm, p7.x-p0.x, p7.y-p0.y, p7.z-p0.z);
		el_sig = (sig ? ptrMesh->sig[r.mtr] : (ptrMesh->sig[r.mtr]-ptrMesh->sig0[r.mtr]));
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

struct loc_source
{
	int elem;					//  
	double x,y,z;				//     
	double dx,dy,dz;			//     
	double ksi,eta,dzeta;		//     
	double d_ksi,d_eta,d_dzeta;	//    
	double cur,len;				//      
	loc_source(){elem=-1;}
};

struct LineSource
{
	PointXYZ A,B;
};

void CheckAndGet(PointXYZ &p,PointXYZ &Pmin,PointXYZ &Pmax,PointXYZ *Hex,PointXYZ &loc,
				 loc_source &ls,const double bnd_loc_min,const double bnd_loc_max,int elem,int &nn)
{
	if(p.x>=Pmin.x && p.x<=Pmax.x && p.y>=Pmin.y && p.y<=Pmax.y && p.z>=Pmin.z && p.z<=Pmax.z)
	{
		nn=CheckInHex(Hex,p,loc);
		if(!nn &&
			(loc.x>bnd_loc_min && loc.x<bnd_loc_max &&
			loc.y>bnd_loc_min && loc.y<bnd_loc_max &&
			loc.z>bnd_loc_min && loc.z<bnd_loc_max ))
		{
			ls.elem=elem;
			ls.ksi=loc.x;
			ls.eta=loc.y;
			ls.dzeta=loc.z;
		}
	}
}

int GetNearestElement(loc_source &ls,Mesh3DForVP &meshvp3d)
{
	double scurr,sbest,m_s_best[8];
	int i,elem,ib,jb,i_best,n_best,m_i_best[8];
	PointXYZ Pmin,Pmax,Pcur,Hex[8],p,loc;
	if(!meshvp3d.krect3d){return -1;}
	p.x=ls.x;
	p.y=ls.y;
	p.z=ls.z;
	n_best=8;
	for(ib=0;ib<n_best;ib++)
	{
		m_s_best[ib]=1e+30;
		m_i_best[ib]=-1;
	}
	sbest=m_s_best[n_best-1];
	if(meshvp3d.krect3d<n_best){n_best=meshvp3d.krect3d;}

	for(elem=0;elem<meshvp3d.krect3d;elem++)
	{
		if(meshvp3d.rect3d[elem].mtr!=1)	//   
		{
			Pcur=meshvp3d.pnt[meshvp3d.rect3d[elem].nodes[0]-1];
			Hex[0]=Pmax=Pmin=Pcur;
			for(i=1;i<8;i++)
			{
				Pcur=meshvp3d.pnt[meshvp3d.rect3d[elem].nodes[i]-1];
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
		for(i=0;i<8;i++){Hex[i]=meshvp3d.pnt[meshvp3d.rect3d[elem].nodes[i]-1];}
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
				ls.elem=elem;
				ls.ksi=loc.x;
				ls.eta=loc.y;
				ls.dzeta=loc.z;
			}
		}
	}

	return i_best;
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

int GetRightPartForRect3D(const Mesh3DForVP* ptrMesh, SparseSLAE* ptrSLAE,loc_source &ls,int ipls,double crdcfflen)
{
	int i, _i, _ti, c;
	double value;


	const Rect3D& r=ptrMesh->rect3d[ls.elem];
	for(i=0;i<8;i++)
	{
		value=phi(i,ls.ksi,ls.eta,ls.dzeta)*crdcfflen;
		_i=r.nodes[i];
		ptrSLAE->AddToVector(_i-1+ipls*ptrSLAE->n,value, ptrSLAE->F);
	}
	
	return RETCODE_OK;
}

double CalcIntA(PointXYZ &A,PointXYZ &B,PointXYZ &P)
{ 
	const int d=20;
	int i;
	PointXYZ AB,pJ;
	AB.x=B.x-A.x;
	AB.y=B.y-A.y;
	AB.z=B.z-A.z;
	double a,dAB;
	a=0.0;
	dAB=sqrt(AB.x*AB.x+AB.y*AB.y+AB.z*AB.z)/d;
	normalize(AB);
	pJ.x=A.x+AB.x*0.5*dAB;
	pJ.y=A.y+AB.y*0.5*dAB;
	pJ.z=A.z+AB.z*0.5*dAB;
	for(i=0;i<d;i++)
	{
		a+=dAB/sqrt((pJ.x-P.x)*(pJ.x-P.x)+(pJ.y-P.y)*(pJ.y-P.y)+(pJ.z-P.z)*(pJ.z-P.z));
		pJ.x+=AB.x*dAB;
		pJ.y+=AB.y*dAB;
		pJ.z+=AB.z*dAB;
	}
	a*=MU0/(4.0*_PI_);
	return a;
}

int AddToAxFromPoissonEq(Mesh3DForVP& meshvp3d,double *ax,int ipls)
{
	for(int i=ipls*meshvp3d.kpnt;i<(ipls+1)*meshvp3d.kpnt;i++)
		ax[i]+=CalcIntA(meshvp3d.A,meshvp3d.B,meshvp3d.pnt[i]);
	return RETCODE_OK;
}

// Calculation of a stationary electric field
int SolveStationarProblem3d(const char* path,bool f3d,bool forA0)
{
	ifstream inf;
	char buf[256];
	int  i,j,k,l,m,retc,ipls,npls;
	FILE *fp;
	int *ig,*jg;
	pardiso_solver prds;
	ofstream outfMN,outfMOON;
	ifstream infmn;
	ofstream outf;
	vector<int> vNrec,vNrecBE;
	int nSrsRazb;
	double val,ImpStatCur;
	int nRecB;

	Mesh3DForVP meshvp3d(moUnloadSolutions);

	SparseSLAE SLAE3D;

	vector<LineSource> lsrs;
	vector<vector<loc_source>> srs;	//    
	vector<loc_source> srsV[2];		//  

	ImpStatCur=1.0;
	inf.open("ImpStatCur");
	if(inf)
	{
		inf>>ImpStatCur;
		inf.close();
	}
	inf.clear();

	retc=GetNumberOfPlaces(npls);
	if(retc)
	{
		cout << "Function GetNumberOfPlaces retrurned " << retc << '\n';
		logfile << "Function GetNumberOfPlaces retrurned " << retc << '\n';
		return 1;
	}

	vNrec.resize(npls);
	inf.open("lin");
	if(!inf)
	{
		logfile<<"Error in open file "<<"lin"<<endl;
		cout<<"Error in open file "<<"lin"<<endl;
		return 1;
	}
	for(i=0;i<npls;i++){inf>>vNrec[i];}
	inf.close();
	inf.clear();

	vNrecBE.resize(npls);
	for(i=0;i<npls;i++){vNrecBE[i]=0;}
	inf.open("recvsbe");
	if(inf)
	{
		for(i=0;i<npls;i++){inf>>vNrecBE[i];}
		inf.close();
	}
	else
	{
		inf.clear();
		inf.open("recvsb");
		if(inf)
		{
			for(i=0;i<npls;i++){inf>>vNrecBE[i];}
			inf.close();
		}
	}
	inf.clear();

	meshvp3d.npls=npls;

	if ((retc=meshvp3d.Read(path))!=0)
		return retc;

	if ((retc=BuildPortrait(&meshvp3d, &SLAE3D, npls))!=0)
		return retc;

	ig=new int[SLAE3D.n+1];
	jg=new int[SLAE3D.jsize];

	for(i=0;i<(SLAE3D.n+1);i++){ig[i]=SLAE3D.iptr[i]-1;}
	for(i=0;i<SLAE3D.jsize;i++){jg[i]=SLAE3D.jptr[i]-1;}

	nSrsRazb=500;
	inf.open("nSrsRazb");
	if(inf)
	{
		inf>>nSrsRazb;
		inf.close();
	}
	inf.clear();

	srs.resize(npls);
	srsV[0].resize(npls);
	srsV[1].resize(npls);

	inf.open("sours");
	if(!inf)
	{
		logfile<<"Error in open file "<<"sours"<<endl;
		cout<<"Error in open file "<<"sours"<<endl;
		return 1;
	}
	lsrs.resize(npls);
	for(i=0;i<npls;i++)
	{
		LineSource &ls=lsrs[i];
		inf>>ls.A.x>>ls.A.y>>ls.A.z;
		inf>>ls.B.x>>ls.B.y>>ls.B.z;
		srs[i].resize(nSrsRazb);
		srsV[0][i].x=ls.A.x;
		srsV[0][i].y=ls.A.y;
		srsV[0][i].z=ls.A.z;
		srsV[1][i].x=ls.B.x;
		srsV[1][i].y=ls.B.y;
		srsV[1][i].z=ls.B.z;
	}
	inf.close();
	inf.clear();

	int elem,nn;
	const double eps_loc_crd=1e-2;
	const double bnd_loc_min=0.0-eps_loc_crd;
	const double bnd_loc_max=1.0+eps_loc_crd;

	for(i=0;i<npls;i++)
	{
		LineSource &gen=lsrs[i];
		k=(int)srs[i].size();
		for(j=0;j<k;j++)
		{
			loc_source &ls=srs[i][j];
			PointXYZ dir,pp,pm;
			pm=gen.A;
			pp=gen.B;
			dir.x=pp.x-pm.x;
			dir.y=pp.y-pm.y;
			dir.z=pp.z-pm.z;
			dir.x/=nSrsRazb;
			dir.y/=nSrsRazb;
			dir.z/=nSrsRazb;
			ls.len=sqrt(dir.x*dir.x+dir.y*dir.y+dir.z*dir.z);
			ls.x=pm.x+(j+0.5)*dir.x;
			ls.y=pm.y+(j+0.5)*dir.y;
			ls.z=pm.z+(j+0.5)*dir.z;
			ls.dx=dir.x;
			ls.dy=dir.y;
			ls.dz=dir.z;
		}
	}

	for(elem=0;elem<meshvp3d.krect3d;elem++)
	{
		if(meshvp3d.rect3d[elem].mtr!=1)	//   
		{
			PointXYZ Pmin,Pmax,Pcur,Hex[8];	// Pmin,Pmax -  

			Pcur=meshvp3d.pnt[meshvp3d.rect3d[elem].nodes[0]-1];
			Hex[0]=Pmax=Pmin=Pcur;
			for(i=1;i<8;i++)
			{
				Pcur=meshvp3d.pnt[meshvp3d.rect3d[elem].nodes[i]-1];
				Hex[i]=Pcur;
				if(Pcur.x<Pmin.x){Pmin.x=Pcur.x;}
				if(Pcur.y<Pmin.y){Pmin.y=Pcur.y;}
				if(Pcur.z<Pmin.z){Pmin.z=Pcur.z;}
				if(Pcur.x>Pmax.x){Pmax.x=Pcur.x;}
				if(Pcur.y>Pmax.y){Pmax.y=Pcur.y;}
				if(Pcur.z>Pmax.z){Pmax.z=Pcur.z;}
			}

			Pmin.x-=1e-3;
			Pmin.y-=1e-3;
			Pmin.z-=1e-3;
			Pmax.x+=1e-3;
			Pmax.y+=1e-3;
			Pmax.z+=1e-3;

			for(i=0;i<npls;i++)
			{
				LineSource &gen=lsrs[i];

				PointXYZ loc;

				loc.x=loc.y=loc.z=0.5;

				k=(int)srs[i].size();
				for(j=0;j<k;j++)
				{
					loc_source &ls=srs[i][j];
					if(ls.elem==-1)
					{
						const double delta=0.1;
						PointXYZ p,dir,pp,pm;

						pm=gen.A;
						pp=gen.B;
						p.x=ls.x;
						p.y=ls.y;
						p.z=ls.z;
						dir.x=ls.dx;
						dir.y=ls.dy;
						dir.z=ls.dz;

						if(p.x>=Pmin.x && p.x<=Pmax.x && p.y>=Pmin.y && p.y<=Pmax.y && p.z>=Pmin.z && p.z<=Pmax.z)
						{
							normalize(dir);
							pm.x=p.x-delta*dir.x;
							pm.y=p.y-delta*dir.y;
							pm.z=p.z-delta*dir.z;
							pp.x=p.x+delta*dir.x;
							pp.y=p.y+delta*dir.y;
							pp.z=p.z+delta*dir.z;

							nn=CheckInHex(Hex,p,loc);
							if(!nn &&
								(loc.x>bnd_loc_min && loc.x<bnd_loc_max &&
								loc.y>bnd_loc_min && loc.y<bnd_loc_max &&
								loc.z>bnd_loc_min && loc.z<bnd_loc_max ))
							{
								ls.elem=elem;
								ls.ksi=loc.x;
								ls.eta=loc.y;
								ls.dzeta=loc.z;
								nn=CheckInHex(Hex,pp,dir);
								if(nn)
								{
									logfile<<"Can't find local coordinates for point ";
									logfile<<pp.x<<' '<<pp.y<<' '<<pp.z<<' '<<endl;
									exit(1);
								}
								ls.d_ksi=dir.x;
								ls.d_eta=dir.y;
								ls.d_dzeta=dir.z;
								nn=CheckInHex(Hex,pm,dir);
								if(nn)
								{
									logfile<<"Can't find local coordinates for point ";
									logfile<<pm.x<<' '<<pm.y<<' '<<pm.z<<' '<<endl;
									exit(1);
								}
								ls.d_ksi-=dir.x;
								ls.d_eta-=dir.y;
								ls.d_dzeta-=dir.z;
								dir.x=ls.d_ksi;
								dir.y=ls.d_eta;
								dir.z=ls.d_dzeta;
								normalize(dir);
								ls.d_ksi=dir.x;
								ls.d_eta=dir.y;
								ls.d_dzeta=dir.z;
							}
						}
					}
				}
				if(srsV[0][i].elem==-1)
				{
					CheckAndGet(lsrs[i].A,Pmin,Pmax,Hex,loc,srsV[0][i],bnd_loc_min,bnd_loc_max,elem,nn);
				}
				if(srsV[1][i].elem==-1)
				{
					CheckAndGet(lsrs[i].B,Pmin,Pmax,Hex,loc,srsV[1][i],bnd_loc_min,bnd_loc_max,elem,nn);
				}
			}
		}
	}


	for(i=0;i<npls;i++)
	{
		k=(int)srs[i].size();
		for(j=0;j<k;j++)
		{
			loc_source &ls=srs[i][j];
			if(ls.elem==-1)
			{
				elem=GetNearestElement(ls,meshvp3d);
				if(elem!=-1)
				{
					PointXYZ p,dir,pp,pm,Hex[8];
					const double delta=0.1;
					for(l=0;l<8;l++){Hex[l]=meshvp3d.pnt[meshvp3d.rect3d[elem].nodes[l]-1];}
					p.x=ls.x;
					p.y=ls.y;
					p.z=ls.z;
					dir.x=ls.dx;
					dir.y=ls.dy;
					dir.z=ls.dz;
					normalize(dir);
					pm.x=p.x-delta*dir.x;
					pm.y=p.y-delta*dir.y;
					pm.z=p.z-delta*dir.z;
					pp.x=p.x+delta*dir.x;
					pp.y=p.y+delta*dir.y;
					pp.z=p.z+delta*dir.z;
					nn=CheckInHex(Hex,pp,dir);
					if(nn)
					{
						logfile<<"Can't find local coordinates for point ";
						logfile<<pp.x<<' '<<pp.y<<' '<<pp.z<<' '<<endl;
						exit(1);
					}
					ls.d_ksi=dir.x;
					ls.d_eta=dir.y;
					ls.d_dzeta=dir.z;
					nn=CheckInHex(Hex,pm,dir);
					if(nn)
					{
						logfile<<"Can't find local coordinates for point ";
						logfile<<pm.x<<' '<<pm.y<<' '<<pm.z<<' '<<endl;
						exit(1);
					}
					ls.d_ksi-=dir.x;
					ls.d_eta-=dir.y;
					ls.d_dzeta-=dir.z;
					dir.x=ls.d_ksi;
					dir.y=ls.d_eta;
					dir.z=ls.d_dzeta;
					normalize(dir);
					ls.d_ksi=dir.x;
					ls.d_eta=dir.y;
					ls.d_dzeta=dir.z;
				}
				else
				{
					logfile<<"Error in find element for sours"<<endl;
					exit(1);
				}
			}
		}
		if(srsV[0][i].elem==-1)
		{
			elem=GetNearestElement(srsV[0][i],meshvp3d);
			if(elem==-1)
			{
				logfile<<"Error in find element for sours"<<endl;
				exit(1);
			}
		}
		if(srsV[1][i].elem==-1)
		{
			elem=GetNearestElement(srsV[1][i],meshvp3d);
			if(elem==-1)
			{
				logfile<<"Error in find element for sours"<<endl;
				exit(1);
			}
		}
	}

	outf.open("sours_dip_calc");
	for(i=0;i<npls;i++)
	{
		k=(int)srs[i].size();
		for(j=0;j<k;j++)
		{
			loc_source &ls=srs[i][j];
			outf<<' '<<ls.x<<' '<<ls.y<<' '<<ls.z<<'\n';
		}
	}
	outf.close();
	outf.clear();


	SLAE3D.A.Clear(SLAE3D.n,SLAE3D.jsize);
	SLAE3D.F.Clear(SLAE3D.n*npls);

	if ((retc=GetLocalContributionsForRect3D(&meshvp3d, &SLAE3D))!=0)
		return retc;

	for(i=0;i<npls;i++)
	{
		GetRightPartForRect3D(&meshvp3d,&SLAE3D,srsV[0][i],i,-1.0*ImpStatCur);
		GetRightPartForRect3D(&meshvp3d,&SLAE3D,srsV[1][i],i,1.0*ImpStatCur);
	}

	SLAE3D.A.AddToMatrix(SLAE3D.n, SLAE3D.jsize, SLAE3D.B);

	ofstream ofp;



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

	infmn.open("xyzmn");
	if(!infmn){
		logfile<<"Error in open file xyzmn"<<endl;
		return 1;
	}
	infmn>>meshvp3d.qMN;
	infmn.close();
	infmn.clear();

	int nRecBE;

	infmn.open("xyzVectorErec");
	if(infmn)
	{
		infmn>>nRecBE;
		infmn.close();
	}
	else
	{
		infmn.clear();
		infmn.open("xyzVectorB");
		if(infmn)
		{
			infmn>>nRecBE;
			infmn.close();
		}
	}
	infmn.clear();

	if(CheckStop())return 1;

	ofp.open("v3n.dat",ios::binary);
	for(i=0;i<meshvp3d.kpnt;i++){meshvp3d.q_all[i]=0.0;}
	for(ipls=0;ipls<npls;ipls++)
	{
		ofp.write((char *)meshvp3d.q_all,meshvp3d.kpnt*sizeof(double));
	}
	ofp.close();
	ofp.clear();

	ofp.open("v3a.dat",ios::binary);
	for(ipls=0;ipls<npls;ipls++)
	{
		if ((retc=meshvp3d.GetWeights(SLAE3D.pU->v+ipls*SLAE3D.n, false))!=0)
			return retc;

		if ((retc=ProcessResCalcPoints(&meshvp3d, meshvp3d.q_all))!=0)
			return retc;

		ofp.write((char *)meshvp3d.q_all,meshvp3d.kpnt*sizeof(double));


		sprintf(buf,"dV_anom.000.%d",ipls+1);
		outfMN.open(buf);
		sprintf(buf,"ddV_anom.000.%d",ipls+1);
		outfMOON.open(buf);
		k=0;
		for(j=ipls-1;j>=0;j--){k+=vNrec[j];}
		for(i=k;i<k+vNrec[ipls];i++)
		{
			const double& vM=meshvp3d.ResCalcPointArray[meshvp3d.boundsV.beginP+3*i  ].res[0].x;
			const double& vN=meshvp3d.ResCalcPointArray[meshvp3d.boundsV.beginP+3*i+1].res[0].x;
			const double& vO=meshvp3d.ResCalcPointArray[meshvp3d.boundsV.beginP+3*i+2].res[0].x;
			outfMN<<1000*(vM-vN)<<'\n';
			outfMOON<<1000*((vM-vO)-(vO-vN))<<'\n';	
		}
		outfMN.close();
		outfMN.clear();
		outfMOON.close();
		outfMOON.clear();

		sprintf(buf,"e3dstat_anom.%d",ipls+1);
		outfMN.open(buf);
		k=0;
		for(j=ipls-1;j>=0;j--){k+=vNrecBE[j];}
		for(i=k;i<k+vNrecBE[ipls];i++)
		{
			PointXYZ pp=meshvp3d.ResCalcPointArray[meshvp3d.boundsVforE.beginP+i].res[0];
			outfMN<<1000*pp.x<<'\t'<<1000*pp.y<<'\t'<<1000*pp.z<<'\n';
		}
		outfMN.close();
		outfMN.clear();
	}
	ofp.close();
	ofp.clear();

	nRecB=0;
	inf.open("xyzVectorB");
	if(inf)
	{
		inf>>nRecB;
		inf.close();
	}
	inf.clear();

	if(nRecB)
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

		vector<int> RecvPlsIgB;

		RecvPlsIgB.resize(npls+1);

		inf.open("recvsb");
		if(!inf)
		{
			cout<<"Error in open file "<<"recvsb"<<endl;
			return 1;
		}
		RecvPlsIgB[0]=0;
		for(i=0;i<npls;i++){inf>>RecvPlsIgB[i+1];}
		inf.close();
		inf.clear();

		for(i=0;i<npls;i++){RecvPlsIgB[i+1]+=RecvPlsIgB[i];}

		OutputResultant3d resultantA(&meshvp3d, vtWithoutDiscontinuity);

		meshvp3d.SetPointsTrue();

		resultantA.Prepare();

		meshvp3d.RestorePointsTrue();

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
			MatrixOnVector(meshvp3d.qv/*meshvp3d.qv_sum*/,SLAE3D.C.adiag, SLAE3D.C.altr, SLAE3D.C.autr, 
							SLAE3D.iptr, SLAE3D.jptr, SLAE3D.n, SLAE3D.tmpV.v);
			SLAE3D.F.AddToVector(SLAE3D.n, SLAE3D.tmpV, -1);
			

			m=(int)srs[ipls].size();
			for(l=0;l<m;l++)
			{
				loc_source &ls=srs[ipls][l];
				val=ls.d_ksi*ls.len*ImpStatCur;
				GetRightPartForRect3D(&meshvp3d,&SLAE3D,ls,0,val);
			}

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
			MatrixOnVector(meshvp3d.qv/*meshvp3d.qv_sum*/,SLAE3D.C.adiag, SLAE3D.C.altr, SLAE3D.C.autr, 
							SLAE3D.iptr, SLAE3D.jptr, SLAE3D.n, SLAE3D.tmpV.v);
			SLAE3D.F.AddToVector(SLAE3D.n, SLAE3D.tmpV, -1);


			m=(int)srs[ipls].size();
			for(l=0;l<m;l++)
			{
				loc_source &ls=srs[ipls][l];
				val=ls.d_eta*ls.len*ImpStatCur;
				GetRightPartForRect3D(&meshvp3d,&SLAE3D,ls,0,val);
			}

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
			MatrixOnVector(meshvp3d.qv/*meshvp3d.qv_sum*/,SLAE3D.C.adiag, SLAE3D.C.altr, SLAE3D.C.autr, 
							SLAE3D.iptr, SLAE3D.jptr, SLAE3D.n, SLAE3D.tmpV.v);
			SLAE3D.F.AddToVector(SLAE3D.n, SLAE3D.tmpV, -1);


			m=(int)srs[ipls].size();
			for(l=0;l<m;l++)
			{
				loc_source &ls=srs[ipls][l];
				val=ls.d_dzeta*ls.len*ImpStatCur;
				GetRightPartForRect3D(&meshvp3d,&SLAE3D,ls,0,val);
			}

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


			meshvp3d.q_res1=ay;
			meshvp3d.q_res2=az;
			resultantA.ValueType=Res3DValueType::vtRotxA;
			resultantA.Output(0);

			meshvp3d.q_res1=ax;
			meshvp3d.q_res2=az;
			resultantA.ValueType=Res3DValueType::vtRotyA;
			resultantA.Output(0);

			meshvp3d.q_res1=ax;
			meshvp3d.q_res2=ay;
			resultantA.ValueType=Res3DValueType::vtRotzA;
			resultantA.Output(0);

			sprintf(buf,"b3dstat_anom.%d",ipls+1);
			ofp.open(buf);
			for(i=RecvPlsIgB[ipls];i<RecvPlsIgB[ipls+1];i++)
			{
				PointForResultCalculationInMesh3D& cpnt=meshvp3d.ResCalcPointArray[meshvp3d.boundsB.beginP+i];
				ofp<<1000*cpnt.res[0].x<<' '<<1000*cpnt.res[0].y<<' '<<1000*cpnt.res[0].z<<'\n';
			}
			ofp.close();
			ofp.clear();
		}
		
		if(u3){delete [] u3;u3=NULL;}

		if(CheckStop())return 1;

		delete [] ax;
		delete [] ay;
		delete [] az;
	}
	else
	{
		for(ipls=0;ipls<npls;ipls++)
		{
			sprintf(buf,"b3dstat_anom.%d",ipls+1);
			ofp.open(buf);
			ofp.close();
			ofp.clear();
		}
	}

	return RETCODE_OK;
}
int CalcSP()
{
	if(!CheckStop())
	{
		SolveStationarProblem3d(".", true, true);
	}
	return 0;
}
