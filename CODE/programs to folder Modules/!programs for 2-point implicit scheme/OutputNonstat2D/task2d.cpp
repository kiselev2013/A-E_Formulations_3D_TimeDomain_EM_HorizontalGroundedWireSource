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
 *  This file contains the code for working with a 2D mesh and solution
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 January 16, 2023
*/

#include "stdafx.h"
#include "task2d.h"
#include "utils.h"
#include "MemFile.h"

extern ofstream logfile;

bool inTriangle(PointRZ *trg,const PointRZ &ptmp)
{
	int i,n1,n2;
	int edges[3][2]={{0,1},{1,2},{2,0}};
	bool fsig[3];
	for(i=0;i<3;i++)
	{
		n1=edges[i][0];
		n2=edges[i][1];
		fsig[i]=((trg[n1].z-trg[n2].z)*ptmp.r+(trg[n2].r-trg[n1].r)*ptmp.z+(trg[n1].r*trg[n2].z-trg[n2].r*trg[n1].z)<=0.0);
	}
	return (fsig[0]==fsig[1] && fsig[1]==fsig[2]);
}

bool CalcBaseFunctions(PointRZ *trg,double a[3][3])
{
	double det;

	a[0][0] = (trg[1].r*trg[2].z-trg[2].r*trg[1].z);
	a[0][1] = (trg[1].z - trg[2].z);
	a[0][2] = (trg[2].r - trg[1].r);


	a[1][0] = (trg[2].r*trg[0].z-trg[0].r*trg[2].z);
	a[1][1] = (trg[2].z - trg[0].z);
	a[1][2] = (trg[0].r - trg[2].r);

	a[2][0] = (trg[0].r*trg[1].z-trg[1].r*trg[0].z);
	a[2][1] = (trg[0].z - trg[1].z);
	a[2][2] = (trg[1].r - trg[0].r);

	det = a[2][2]*a[1][1]- a[1][2]*a[2][1];

	if(fabs(det)>1e-12)
	{
		a[0][0] /= det;	a[0][1] /= det;	a[0][2] /= det;
		a[1][0] /= det;	a[1][1] /= det;	a[1][2] /= det;
		a[2][0] /= det;	a[2][1] /= det;	a[2][2] /= det;
		return 1;
	}
	else
	{
		return 0;
	}
}

double scal(const PointRZ &p1,const PointRZ &p2)
{
	return p1.r*p2.r+p1.z*p2.z;
}

double dist(const PointRZ &p1,const PointRZ &p2)
{
	return sqrt((p2.r-p1.r)*(p2.r-p1.r)+(p2.z-p1.z)*(p2.z-p1.z));
}

task2d::task2d()
{
	rect=NULL;
	pnt=NULL;

	ntimes=0;

	nmat3d=0;
	times=NULL;
	mtr3d2d=NULL;

	nmat2d=0;
	sigma=NULL;
	sigmaZ=NULL;
	mu=NULL;

	reg=NULL;
	rm=NULL;
	zm=NULL;
}

void task2d::SetPath(char *_path)
{
	strcpy(path,_path);
}

task2d::~task2d()
{
	if(rect){delete [] rect; rect=NULL;}
	if(pnt){delete [] pnt; pnt=NULL;}
	if(times){delete [] times; times=NULL;}
	if(mtr3d2d){delete [] mtr3d2d; mtr3d2d=NULL;}
	if(sigma){delete [] sigma; sigma=NULL;}
	if(sigmaZ){delete [] sigmaZ; sigmaZ=NULL;}
	if(mu){delete [] mu; mu=NULL;}
	if(reg){delete [] reg; reg=NULL;}
	if(rm){delete [] rm; rm=NULL;}
	if(zm){delete [] zm; zm=NULL;}
}

/*!    */
int task2d::ReadTimesFile()
{
	sprintf(buf,"%s\\infite.0",path);
	inf.open(buf);
	int i;
	if (!inf) return RETCODE_NOFILE;
	inf.ignore(1000, '=');
	inf >> ntimes;
	times=new double[ntimes];
	inf.ignore(1000, '\n');
	inf.ignore(1000, '\n');
	inf.ignore(1000, '\n');
	inf.ignore(1000, '\n');
	for(i=0;i<ntimes;i++)
	{
		if (!inf.good()) return RETCODE_BADFILE;
		inf >> times[i];
		inf.ignore(1000, ';');
	}
	inf.close();
	inf.clear();
	return RETCODE_OK;
}

int task2d::ReadMtr3d2d(char *path3d)
{
	int j,k;

	sprintf(buf,"%s\\mtr3d2d",path3d);
	inf.open(buf);
	if(!inf)return RETCODE_NOFILE;
	nmat3d=0;
	while(!inf.eof())
	{
		inf>>k;
		if(inf.eof() || !inf.good())break;
		inf>>j;
		if(k>nmat3d){nmat3d=k;}
	}
	inf.close();
	inf.clear();

	if(nmat3d)
	{
		mtr3d2d=new int[nmat3d];
		if(!mtr3d2d)return RETCODE_NOMEM;
		sprintf(buf,"%s\\mtr3d2d",path3d);
		inf.open(buf);
		if(!inf)return RETCODE_NOFILE;
		nmat3d=0;
		while(!inf.eof())
		{
			inf>>k;
			if(inf.eof() || !inf.good())break;
			inf>>j;
			mtr3d2d[k-1]=j;
		}
		inf.close();
		inf.clear();
	}

	return 0;
}

int task2d::Read(int fSigmaZ)
{
	int i,k,retc;
	ifstream rzf, nvtrf, nvkatf, l1f, nvk1f;
	double sum;

	sprintf(buf,"%s\\inf2tr.dat",path);
	inf.open(buf);
	if (!inf) return RETCODE_NOFILE;
	inf.ignore(1000, '\n');
	inf.ignore(1000, '='); inf>>kpnt;
	inf.ignore(1000, '='); inf>>krect;
	inf.close();
	inf.clear();

	sprintf(buf,"%s\\rz.dat",path);
	rzf.open(buf, ios::binary);
	if (!rzf) return RETCODE_NOFILE;
	pnt=new PointRZ[kpnt];
	if (!pnt) return RETCODE_NOMEM;
	for (i=0; i<kpnt; i++)
	{
		rzf > pnt[i].r;
		rzf > pnt[i].z;
		if(fabs(pnt[i].z)<1e-6)
		{
			pnt[i].z=0.0;
		}
	}
	rzf.close();
	rzf.clear();

	sprintf(buf,"%s\\nvtr.dat",path);
	nvtrf.open(buf, ios::binary);
	if (!nvtrf) return RETCODE_NOFILE;
	rect=new Rect[krect];
	if (!rect) return RETCODE_NOMEM;
	for (i=0; i<krect; i++)
	{
		nvtrf > rect[i].nodes[2] > rect[i].nodes[3] > rect[i].nodes[0] > rect[i].nodes[1] > rect[i].nodes[4] > rect[i].rtype;
	}
	nvtrf.close();
	nvtrf.clear();

	sprintf(buf,"%s\\nvkat2d.dat",path);
	nvkatf.open(buf, ios::binary);
	if (!nvkatf) return RETCODE_NOFILE;
	for (i=0; i<krect; i++)
	{
		nvkatf > rect[i].mtr;
	}
	nvkatf.close();
	nvkatf.clear();

	nc=kpnt;

	retc=ReadTimesFile();
	if (retc)
	{
		cout<<"ReadTimesFile returned "<<retc<<endl;
		logfile<<"ReadTimesFile returned "<<retc<<endl;
		return retc;
	}

	retc=ReadDecadeInfoTimeAll(ntimes,times,nDec,DecIg,tbeg,tend,path);
	if(retc)
	{
		cout<<"ReadDecadeInfo returned "<<retc<<endl;
		logfile<<"ReadDecadeInfo returned "<<retc<<endl;
		return retc;
	}

	ntimesdec=tend-tbeg;

	sprintf(buf,"%s\\sigma",path);
	inf.open(buf);
	if(!inf) return RETCODE_NOFILE;
	nmat2d=0;
	while(!inf.eof())
	{
		inf>>k;
		if(inf.eof() || !inf.good())break;
		inf>>sum;
		if(k>nmat2d)nmat2d=k;
	}
	inf.close();
	inf.clear();

	sigma=new double[nmat2d];
	sigmaZ=new double[nmat2d];
	mu=new double[nmat2d];

	sprintf(buf,"%s\\sigma",path);
	inf.open(buf);
	if(!inf) return RETCODE_NOFILE;
	while(!inf.eof())
	{
		inf>>k;
		if(inf.eof() || !inf.good())break;
		inf>>sum;
		sigma[k-1]=sum;
	}
	inf.close();
	inf.clear();

	if(fSigmaZ)
	{
		sprintf(buf,"%s\\sigmaZ",path);
		inf.open(buf);
		if(!inf) return RETCODE_NOFILE;
		while(!inf.eof())
		{
			inf>>k;
			if(inf.eof() || !inf.good())break;
			inf>>sum;
			sigmaZ[k-1]=sum;
		}
		inf.close();
		inf.clear();
	}
	else
	{
		for(i=0;i<nmat2d;i++)
		{
			sigmaZ[i]=sigma[i];
		}
	}

	sprintf(buf,"%s\\mu",path);
	inf.open(buf);
	if(!inf) return RETCODE_NOFILE;
	while(!inf.eof())
	{
		inf>>k;
		if(inf.eof() || !inf.good())break;
		inf>>sum;
		mu[k-1]=sum;
	}
	inf.close();
	inf.clear();


	sprintf(buf,"%s\\rz.txt",path);
	inf.open(buf);
	if (!inf) return RETCODE_NOFILE;
	inf>>qr>>qz;
	inf.close();
	inf.clear();

	nreg=krect;

	if(!(reg=new int[nreg])) return RETCODE_NOMEM;
	if(!(rm=new double[qr])) return RETCODE_NOMEM;
	if(!(zm=new double[qz])) return RETCODE_NOMEM;

	sprintf(buf,"%s\\r.dat",path);
	inf.open(buf, ios::binary);
	if(!inf) return RETCODE_NOFILE;
	for(i=0;i<qr;i++){inf>rm[i];}
	inf.close();
	inf.clear();

	sprintf(buf,"%s\\z.dat",path);
	inf.open(buf, ios::binary);
	if(!inf) return RETCODE_NOFILE;
	for(i=0;i<qz;i++)
	{
		inf>zm[i];
		if(fabs(zm[i])<1e-6)
		{
			zm[i]=0.0;
		}
	}
	inf.close();
	inf.clear();

	for(i=0;i<nreg;i++){reg[i]=i+1;}

	return RETCODE_OK;
}

int	task2d::GetField(int it,double *A,char *str)
{
	FILE *fp;
	int ipls;
	sprintf(str,"%s\\v2.%d",path,it+tbeg);
	fp=fopen(str,"rb");
	if(!fp) return RETCODE_NOFILE;
	fread(A,sizeof(double),nc*npls,fp);
	fclose(fp);
	return 0;
}


void _diff_t2(double *du_j,double *u_j,double *u_j1,int n,double *time,int tnum)
{
	int i;
	double invdt;
	invdt=1.0/(time[tnum]-time[tnum-1]);
	for(i=0;i<n;i++)du_j[i]=-((u_j[i]-u_j1[i])*invdt);
}

void _diff_t3(double *du_j,double *u_j,double *u_j1,double *u_j2,int n,double *time,int tnum)
{
	int i;
	double dt,dt0,dt1;
	double mt0,mt1,mt2;

	dt =  time[tnum]  - time[tnum-2];
	dt0 = time[tnum]  - time[tnum-1];
	dt1 = time[tnum-1] - time[tnum-2];

	mt0 = (dt + dt0)/(dt*dt0);
	mt1 = -dt/(dt1*dt0);
	mt2 = dt0/(dt*dt1);
	
	for(i=0;i<n;i++)du_j[i]=-(mt0*u_j[i]+mt1*u_j1[i]+mt2*u_j2[i]);
}

void _diff_t3(double *du_j,double *u_j,double *u_j1,double *u_j2,int n,double t,double t1,double t2)
{
	int i;
	double dt,dt0,dt1;
	double mt0,mt1,mt2;

	dt =  t  - t2;
	dt0 = t  - t1;
	dt1 = t1 - t2;

	mt0 = (dt + dt0)/(dt*dt0);
	mt1 = -dt/(dt1*dt0);
	mt2 = dt0/(dt*dt1);
	
	for(i=0;i<n;i++)du_j[i]=-(mt0*u_j[i]+mt1*u_j1[i]+mt2*u_j2[i]);
}

double get_angle(double px,double py,double ox,double oy)
{		
	double x[2],r;
	x[0]=px-ox;
	x[1]=py-oy;
	r=sqrt(x[0]*x[0]+x[1]*x[1]);
	if(!r)return 0;
	x[0]/=r;
	x[1]/=r;
	if(!x[1]){
		if(x[0]>0)return 0;
		else return PI;
	}
	if(!x[0]){
		if(x[1]>0)return PI/2;
		else return 3*PI/2;
	}
	if(x[0]>0 && x[1]>0)return atan2(x[1],x[0]);
	if(x[0]<0 && x[1]>0)return PI-atan2(x[1],-x[0]);
	if(x[0]<0 && x[1]<0)return PI+atan2(-x[1],-x[0]);
	return 2*PI-atan2(-x[1],x[0]);
}

void AddReciver(PointXYZ &xyzVectorPi,PointXYZ &GenMin,PointRZ &recPi,double &sinfiPi,double &cosfiPi,double rc0)
{
	double x,y,r;
	x=xyzVectorPi.x-GenMin.x;
	y=xyzVectorPi.y-GenMin.y;
	r=sqrt(x*x+y*y);
	if(r<rc0)r=rc0;
	recPi.r=r;
	recPi.z=xyzVectorPi.z;
	sinfiPi=y/r;
	cosfiPi=x/r;
}

void RollFields(double *(&u),double *(&u_1),double *(&u_2))
{
	double *tmpu;
	tmpu=u_2;
	u_2=u_1;
	u_1=u;
	u=tmpu;
}

void ClearVector(float *(&v)){if(v){delete [] v; v=NULL;}}
void ClearVector(double *(&v)){if(v){delete [] v; v=NULL;}}
void ClearVector(PointRZ *(&v)){if(v){delete [] v; v=NULL;}}
void ClearVector(PointXYZ *(&v)){if(v){delete [] v; v=NULL;}}
void ClearVector2d(double **(&v),int n)
{
	int i;
	for(i=0;i<n;i++)
	{
		ClearVector(v[i]);
	}
}


void task2d_Ax::GetSource(PointXYZ &GenMin,PointXYZ &GenMax,PointXYZ &pA,PointXYZ &pB,PointXYZ &vAB,PointXYZ &Ic,double &len,double &cosphi,double &sinphi)
{
	double phi,cedr;
	PointXYZ pT;

	pA.x=GenMin.x;
	pA.y=GenMin.y;
	pA.z=GenMin.z;

	pB.x=GenMax.x;
	pB.y=GenMax.y;
	pB.z=GenMin.z;

	vAB.x=pB.x-pA.x;
	vAB.y=pB.y-pA.y;
	vAB.z=pB.z-pA.z;

	cedr=sqrt(vAB.x*vAB.x+vAB.y*vAB.y+vAB.z*vAB.z);
	phi=acos(vAB.x/cedr);

	if(vAB.y<0) phi*=-1.;

	cosphi=cos(-phi);
	sinphi=sin(-phi);
	
	pB.x=pB.x-pA.x;
	pB.y=pB.y-pA.y;

	pT.x=pB.x*cosphi-pB.y*sinphi;
	pT.y=pB.x*sinphi+pB.y*cosphi;

	pB.x=pT.x+pA.x;
	pB.y=pT.y+pA.y;

	vAB.x=pB.x-pA.x;
	vAB.y=pB.y-pA.y;
	vAB.z=pB.z-pA.z;

	Ic=vAB;

	Ic.x/=cedr;
	Ic.y/=cedr;
	Ic.z/=cedr;

	len=cedr/NUMBEROFABDIPOLES;

	vAB.x/=NUMBEROFABDIPOLES;
	vAB.y/=NUMBEROFABDIPOLES;
	vAB.z/=NUMBEROFABDIPOLES;
}

int task2d_Ax::init()
{
	NUMBEROFABDIPOLES=20;

	nthreads = omp_get_max_threads();
	nthreads = (nthreads>3)? 3 : nthreads;
	omp_set_num_threads(nthreads);

	PathTo2d[0]='\0';
	inf.open("PathTo2d");
	if(inf)
	{
		inf>>PathTo2d;
		logfile<<"Outputing from "<<PathTo2d<<'\n';
		inf.close();
	}
	else
	{
		strcpy(PathTo2d,".");
	}
	inf.clear();
	strcat(PathTo2d,"\\Ax");

	inf.open("numberofabdipoles");
	if(!inf)
	{
		logfile<<"Error in open file "<<"numberofabdipoles"<<endl;
		cout<<"Error in open file "<<"numberofabdipoles"<<endl;
		return 1;	
	}
	inf>>NUMBEROFABDIPOLES;
	inf.close();
	inf.clear();

	GetNumberOfPlaces(npls);

	nGenByPls.resize(npls);
	NrecB.resize(npls);
	NrecE.resize(npls);

	inf.open("srsclcgsz");
	if(!inf)
	{
		logfile<<"Error in open file "<<"srsclcgsz"<<endl;
		cout<<"Error in open file "<<"srsclcgsz"<<endl;
		return 1;
	}
	p1=0;
	for(i=0;i<npls;i++)
	{
		inf>>nGenByPls[i];
		p1+=nGenByPls[i];
	}
	inf.close();
	inf.clear();


	RecvPlsIgB.resize(p1+1);
	RecvPlsIgE.resize(p1+1);
	RecvPlsIgE0.resize(p1+1);

	inf.open("xyzVectorE0");
	if(inf)
	{
		inf>>npntE0;
		inf.close();
	}
	else
	{
		npntE0=0;
	}
	inf.clear();

	if(npntE0)
	{
		inf.open("TgCompE0");
		if(!inf)
		{
			logfile<<"Error in open file "<<"TgCompE0"<<endl;
			cout<<"Error in open file "<<"TgCompE0"<<endl;
			return 1;
		}
		TgCompE0.resize(npntE0);
		for(i=0;i<npntE0;i++)
		{
			inf>>TgCompE0[i].x>>TgCompE0[i].y>>TgCompE0[i].z;
		}
		inf.close();
		inf.clear();
	}

	inf.open("recvsb");
	if(!inf)
	{
		logfile<<"Error in open file "<<"recvsb"<<endl;
		cout<<"Error in open file "<<"recvsb"<<endl;
		return 1;
	}
	RecvPlsIgB[0]=0;
	m=0;
	for(i=0;i<npls;i++)
	{
		inf>>NrecB[i];
		if(nGenByPls[i])
		{
			RecvPlsIgB[m+1]=NrecB[i];
			m++;
		}
	}
	inf.close();
	inf.clear();

	inf.open("recvse");
	if(!inf)
	{
		logfile<<"Error in open file "<<"recvse"<<endl;
		cout<<"Error in open file "<<"recvse"<<endl;
		return 1;
	}
	RecvPlsIgE[0]=0;
	m=0;
	for(i=0;i<npls;i++)
	{
		inf>>NrecE[i];
		if(nGenByPls[i])
		{
			RecvPlsIgE[m+1]=NrecE[i];
			m++;
		}
	}
	inf.close();
	inf.clear();

	m=0;
	RecvPlsIgE0[0]=0;
	for(i=0;i<npls;i++)
	{
		if(nGenByPls[i])
		{
			RecvPlsIgE0[m+1]=npntE0;
			m++;
		}
	}

	p2=p1;
	j=npls;
	for(i=m;i>0;i--)
	{
		while(j>0 && !nGenByPls[j-1]){j--;}
		RecvPlsIgB[p2]=RecvPlsIgB[i];
		for(k=p2-1;k>p2-nGenByPls[j-1];k--){RecvPlsIgB[k]=RecvPlsIgB[p2];}
		p2-=nGenByPls[j-1];
		j--;
	}


	p2=p1;
	j=npls;
	for(i=m;i>0;i--)
	{
		while(j>0 && !nGenByPls[j-1]){j--;}
		RecvPlsIgE[p2]=RecvPlsIgE[i];
		for(k=p2-1;k>p2-nGenByPls[j-1];k--){RecvPlsIgE[k]=RecvPlsIgE[p2];}
		p2-=nGenByPls[j-1];
		j--;
	}

	p2=p1;
	j=npls;
	for(i=m;i>0;i--)
	{
		while(j>0 && !nGenByPls[j-1]){j--;}
		RecvPlsIgE0[p2]=RecvPlsIgE0[i];
		for(k=p2-1;k>p2-nGenByPls[j-1];k--){RecvPlsIgE0[k]=RecvPlsIgE0[p2];}
		p2-=nGenByPls[j-1];
		j--;
	}

	npls_main=npls;
	npls=p1;

	for(i=0;i<npls;i++)
	{
		RecvPlsIgB[i+1]+=RecvPlsIgB[i];
		RecvPlsIgE[i+1]+=RecvPlsIgE[i];
		RecvPlsIgE0[i+1]+=RecvPlsIgE0[i];
	}

	npntB=RecvPlsIgB[npls];
	npntE=RecvPlsIgE[npls];
	npntE0=RecvPlsIgE0[npls];

	RecToSourceB.resize(npntB);
	RecToSourceE.resize(npntE);
	RecToSourceE0.resize(npntE0);

	for(i=0;i<npls;i++)
	{
		for(k=RecvPlsIgB[i];k<RecvPlsIgB[i+1];k++)
		{
			RecToSourceB[k]=i;
		}
	}

	for(i=0;i<npls;i++)
	{
		for(k=RecvPlsIgE[i];k<RecvPlsIgE[i+1];k++)
		{
			RecToSourceE[k]=i;
		}
	}

	for(i=0;i<npls;i++)
	{
		for(k=RecvPlsIgE0[i];k<RecvPlsIgE0[i+1];k++)
		{
			RecToSourceE0[k]=i;
		}
	}

	GenSq.resize(npls);
	inf.open("srsclca");
	if(!inf)
	{
		logfile<<"Error in open file "<<"srsclca"<<endl;
		cout<<"Error in open file "<<"srsclca"<<endl;
		return 1;
	}
	for(i=0;i<npls;i++){inf>>GenSq[i].A.x>>GenSq[i].A.y>>GenSq[i].A.z>>GenSq[i].B.x>>GenSq[i].B.y>>GenSq[i].B.z;}
	inf.close();
	inf.clear();

	inf.open("timeintervalforprint");
	if(!inf)
	{
		logfile<<"No file "<<"timeintervalforprint"<<'\n';
		return 1;
	}
	inf>>ftime>>ltime;
	inf.close();
	inf.clear();

	time_gaps_for_loop=0;

	IcB.resize(npntB);
	IcE.resize(npntE);
	IcE0.resize(npntE0);
	cosphiB.resize(npntB);
	sinphiB.resize(npntB);
	cosphiE.resize(npntE);
	sinphiE.resize(npntE);
	cosphiE0.resize(npntE0);
	sinphiE0.resize(npntE0);
	lenB.resize(npntB);
	lenE.resize(npntE);
	lenE0.resize(npntE0);

	SetPath(PathTo2d);

	retp=Read(0);
	if(retp)
	{
		logfile<<"Function Read() for Ax returned "<<retp<<'\n';
		return 1;
	}

	retp=ReadMtr3d2d(".");
	if(retp)
	{
		logfile<<"Function ReadMtr3d2d() for . returned "<<retp<<'\n';
		return 1;
	}

	rc0=0.5*(rm[0]+rm[1]);

	DminSrc=rc0;
	inf.open("dminsrc");
	if(inf)
	{
		inf>>DminSrc;
		inf.close();
	}
	inf.clear();

	inf.open("xyzVectorB");
	if(!inf)
	{
		logfile<<"No file "<<"xyzVectorB"<<'\n';
		return 1;
	}
	inf>>npntBor;
	xyzVectorB=new PointXYZ[npntBor];
	recB=new PointRZ[npntB*NUMBEROFABDIPOLES];
	sinfiB=new double[npntB*NUMBEROFABDIPOLES];
	cosfiB=new double[npntB*NUMBEROFABDIPOLES];
	npc=0;
	l=0;
	for(ipm=0;ipm<npls_main;ipm++)
	{
		if(!nGenByPls[ipm])
		{
			for(i=0;i<NrecB[ipm];i++)
			{
				inf>>x>>y>>z;
			}
		}
		m=l;
		for(ipls=0;ipls<nGenByPls[ipm];ipls++)
		{
			for(i=RecvPlsIgB[npc+ipls];i<RecvPlsIgB[npc+ipls+1];i++)
			{
				GetSource(GenSq[RecToSourceB[i]].A,GenSq[RecToSourceB[i]].B,pA,pB,vAB,Ic,len,cosphi,sinphi);
				IcB[i]=Ic;
				cosphiB[i]=cosphi;
				sinphiB[i]=sinphi;
				lenB[i]=len;
				if(!ipls)
				{
					inf>>Tp.x>>Tp.y>>Tp.z;
					xyzVectorB[l]=Tp;
					l++;
				}
				else
				{
					Tp=xyzVectorB[m+i-RecvPlsIgB[npc+ipls]];
				}
				for(idp=0;idp<NUMBEROFABDIPOLES;idp++)
				{
					k=i*NUMBEROFABDIPOLES+idp;
					x=Tp.x;
					y=Tp.y;

					x=x-pA.x;
					y=y-pA.y;
					pT.x=x*cosphi-y*sinphi;
					pT.y=x*sinphi+y*cosphi;	
					x=pT.x+pA.x;
					y=pT.y+pA.y;

					pT.x=pA.x+vAB.x*(idp+0.5);
					pT.y=pA.y+vAB.y*(idp+0.5);
					x-=pT.x;
					y-=pT.y;
					r=sqrt(x*x+y*y);
					if(r<rc0)r=rc0;
					recB[k].r=r;
					recB[k].z=Tp.z;
					sinfiB[k]=y/r;
					cosfiB[k]=x/r;
				}
			}
		}
		npc+=nGenByPls[ipm];
	}
	inf.close();
	inf.clear();


	inf.open("xyzVectorE");
	if(!inf)
	{
		logfile<<"No file "<<"xyzVectorE"<<'\n';
		return 1;
	}
	inf>>npntEor;
	xyzVectorE=new PointXYZ[npntEor];
	recE=new PointRZ[npntE*NUMBEROFABDIPOLES];
	sinfiE=new double[npntE*NUMBEROFABDIPOLES];
	cosfiE=new double[npntE*NUMBEROFABDIPOLES];
	matrec=new int[npntE*NUMBEROFABDIPOLES];
	npc=0;
	l=0;
	for(ipm=0;ipm<npls_main;ipm++)
	{
		if(!nGenByPls[ipm])
		{
			for(i=0;i<NrecE[ipm];i++)
			{
				inf>>x>>y>>z;
			}
		}
		m=l;
		for(ipls=0;ipls<nGenByPls[ipm];ipls++)
		{
			for(i=RecvPlsIgE[npc+ipls];i<RecvPlsIgE[npc+ipls+1];i++)
			{
				GetSource(GenSq[RecToSourceE[i]].A,GenSq[RecToSourceE[i]].B,pA,pB,vAB,Ic,len,cosphi,sinphi);
				IcE[i]=Ic;
				cosphiE[i]=cosphi;
				sinphiE[i]=sinphi;
				lenE[i]=len;
				if(!ipls)
				{
					inf>>Tp.x>>Tp.y>>Tp.z;
					xyzVectorE[l]=Tp;
					l++;
				}
				else
				{
					Tp=xyzVectorE[m+i-RecvPlsIgE[npc+ipls]];
				}
				for(idp=0;idp<NUMBEROFABDIPOLES;idp++)
				{
					k=i*NUMBEROFABDIPOLES+idp;
					x=Tp.x;
					y=Tp.y;

					x=x-pA.x;
					y=y-pA.y;
					pT.x=x*cosphi-y*sinphi;
					pT.y=x*sinphi+y*cosphi;	
					x=pT.x+pA.x;
					y=pT.y+pA.y;

					pT.x=pA.x+vAB.x*(idp+0.5);
					pT.y=pA.y+vAB.y*(idp+0.5);
					x-=pT.x;
					y-=pT.y;
					r=sqrt(x*x+y*y);
					if(r<rc0)r=rc0;
					recE[k].r=r;
					recE[k].z=Tp.z;
					sinfiE[k]=y/r;
					cosfiE[k]=x/r;
					matrec[k]=0;
				}
			}
		}
		npc+=nGenByPls[ipm];
	}
	inf.close();
	inf.clear();

	npntE0or=0;
	if(npntE0)
	{
		inf.open("xyzVectorE0");
		if(!inf)
		{
			logfile<<"No file "<<"xyzVectorE0"<<'\n';
			return 1;
		}
		inf>>npntE0or;
		xyzVectorE0=new PointXYZ[npntE0or];
		for(i=0;i<npntE0or;i++)
		{
			inf>>Tp.x>>Tp.y>>Tp.z;
			xyzVectorE0[i]=Tp;
		}
		inf.close();
		inf.clear();
	}

	recE0=new PointRZ[npntE0*NUMBEROFABDIPOLES];
	sinfiE0=new double[npntE0*NUMBEROFABDIPOLES];
	cosfiE0=new double[npntE0*NUMBEROFABDIPOLES];
	matrec0=new int[npntE0*NUMBEROFABDIPOLES];
	npc=0;
	l=0;
	for(ipls=0;ipls<npls;ipls++)
	{
		l=0;
		for(i=RecvPlsIgE0[npc+ipls];i<RecvPlsIgE0[npc+ipls+1];i++)
		{
			GetSource(GenSq[RecToSourceE0[i]].A,GenSq[RecToSourceE0[i]].B,pA,pB,vAB,Ic,len,cosphi,sinphi);
			IcE0[i]=Ic;
			cosphiE0[i]=cosphi;
			sinphiE0[i]=sinphi;
			lenE0[i]=len;
			Tp=xyzVectorE0[l];
			l++;
			for(idp=0;idp<NUMBEROFABDIPOLES;idp++)
			{
				k=i*NUMBEROFABDIPOLES+idp;
				x=Tp.x;
				y=Tp.y;

				x=x-pA.x;
				y=y-pA.y;
				pT.x=x*cosphi-y*sinphi;
				pT.y=x*sinphi+y*cosphi;	
				x=pT.x+pA.x;
				y=pT.y+pA.y;

				pT.x=pA.x+vAB.x*(idp+0.5);
				pT.y=pA.y+vAB.y*(idp+0.5);
				x-=pT.x;
				y-=pT.y;
				r=sqrt(x*x+y*y);
				if(r<rc0)r=rc0;
				recE0[k].r=r;
				recE0[k].z=Tp.z;
				sinfiE0[k]=y/r;
				cosfiE0[k]=x/r;
				matrec0[k]=0;
			}
		}
	}

	inf.close();
	inf.clear();



	nt=ntimesdec;

	A=new double *[3];
	Bx=new double*[nt];
	By=new double*[nt];
	Bz=new double*[nt];
	Fx=new double*[nt];
	Fy=new double*[nt];
	Fz=new double*[nt];
	E=new double*[nt];
	H=new double[kpnt*npls];
	E0=new double[npntE0*NUMBEROFABDIPOLES];

	for(it=0;it<nt;it++)
	{
		Bx[it]=new double[npntB*NUMBEROFABDIPOLES];
		By[it]=new double[npntB*NUMBEROFABDIPOLES];
		Bz[it]=new double[npntB*NUMBEROFABDIPOLES];
		Fx[it]=new double[npntB*NUMBEROFABDIPOLES];
		Fy[it]=new double[npntB*NUMBEROFABDIPOLES];
		Fz[it]=new double[npntB*NUMBEROFABDIPOLES];
		E[it]=new double[npntE*NUMBEROFABDIPOLES];
	}

	retp=SobjB.Init(npntB*NUMBEROFABDIPOLES,kpnt,krect,recB,pnt,rect,0,NULL,sigma,
		sigmaZ,mu,nreg,qr,qz,reg,rm,zm,0);
	if(retp)
	{
		logfile<<"SobjB.Init() returned "<<retp<<'\n';
		return retp;
	}

	retp=SobjE.Init(npntE*NUMBEROFABDIPOLES,kpnt,krect,recE,pnt,rect,0,NULL,sigma,
		sigmaZ,mu,nreg,qr,qz,reg,rm,zm,0);
	if(retp)
	{
		logfile<<"SobjE.Init() returned "<<retp<<'\n';
		return retp;
	}
	
	vRecOutData = new RectOutData[npntE0*NUMBEROFABDIPOLES];
	if(!vRecOutData)return RETCODE_NOMEM;

	rmin = zmin = 1e+30;
	rmax = zmax = -1e+30;

	for(i=0;i<npntE0*NUMBEROFABDIPOLES;i++)
	{
		if(recE0[i].r<rmin){rmin=recE0[i].r; }
		if(recE0[i].r>rmax){rmax=recE0[i].r; }
		if(recE0[i].z<zmin){zmin=recE0[i].z; }
		if(recE0[i].z>zmax){zmax=recE0[i].z; }
	}

	retp=FindElementsForReceivers(npntE0*NUMBEROFABDIPOLES,qr,qz,rm,zm,reg,recE0,vRecOutData,pnt,rect,sigma,sigmaZ,rmin,zmin,rmax,zmax);
	if (retp)
	{
		logfile << "FindElementsForReceivers returned " << retp << endl;
		return retp;
	}
	

	for(it=0;it<3;it++)
	{
		A[2-it]=new double[kpnt*npls];
		if(!A[2-it])return RETCODE_NOMEM;
		GetField(it,A[2-it],buf);
	}

	return 0;
}

void task2d_Ax::output(int it,vector<float> &enor)
{
	if(tbeg+it>2)
	{
		RollFields(A[0],A[1],A[2]);
		GetField(it,A[0],buf);
	}

	double *A0=(tbeg+it>=2)? A[0] : (tbeg+it==1)? A[1] : A[2];

	SobjB.OutputB(A0,Bx[it],By[it],Bz[it],IcB,sinfiB,cosfiB,1,RecToSourceB);

	if(tbeg+it>=1)
	{
		_diff_t2(Fx[it],Bx[it],Bx[it-1],npntB*NUMBEROFABDIPOLES,times,it);
		_diff_t2(Fy[it],By[it],By[it-1],npntB*NUMBEROFABDIPOLES,times,it);
		_diff_t2(Fz[it],Bz[it],Bz[it-1],npntB*NUMBEROFABDIPOLES,times,it);
	}
	else
	{
		k=npntB*NUMBEROFABDIPOLES;
		for(i=0;i<k;i++){Fx[it][i]=Fy[it][i]=Fz[it][i]=0.0;}
	}

	if(tbeg+it>=1)
	{
		_diff_t2(H,A[0],A[1],kpnt*npls,times+tbeg,it);
	}
	else
	{
		k=npntE*NUMBEROFABDIPOLES;
		for(i=0;i<k;i++){E[it][i]=0.0;}
		k=npntE0*NUMBEROFABDIPOLES;
		for(i=0;i<k;i++){E0[i]=0.0;}
	}

	if(tbeg+it>=1)
	{
		SobjE.OutputE(H,E[it],1,RecToSourceE);
		Output_Field_Pr(npntE0*NUMBEROFABDIPOLES,H,E0,vRecOutData,RecToSourceE0,kpnt);
	}

	if(npntE0or)
	{
		kk=0;
		for(ipls=0;ipls<npls;ipls++)
		{
			l=0;
			for(i=RecvPlsIgE0[ipls];i<RecvPlsIgE0[ipls+1];i++)
			{
				Ic=IcE0[i];
				cosphi=cosphiE0[i];
				sinphi=-sinphiE0[i];
				len=lenE0[i];

				U[0]=0.0;
				for(idp=0;idp<NUMBEROFABDIPOLES;idp++)
				{
					k=i*NUMBEROFABDIPOLES+idp;
					U[0]+=E0[k];
				}
				U[0]*=len;

				ux=U[0]*Ic.x;
				uy=U[0]*Ic.y;
				uz=U[0]*Ic.z;

				x=ux*cosphi-uy*sinphi;
				y=ux*sinphi+uy*cosphi;
		
				ux=x;
				uy=y;

				if(recE0[i].r>DminSrc || (fabs(recE0[i].z-GenSq[ipls].A.z)>DminSrc && fabs(recE0[i].z-GenSq[ipls].B.z)>DminSrc))
				{
					enor[kk]=float(ux*TgCompE0[l].x+uy*TgCompE0[l].y+uz*TgCompE0[l].z);
				}
				else
				{
					enor[kk]=0.0;
				}
				l++;
				kk++;
			}
		}

	}
}

void task2d_Ax::output_rec()
{
	sprintf_s(buf, "ball.ax.stat");
	ofpstat.open(buf);
	sprintf_s(buf, "ball.ax");
	ofp.open(buf);
	ofp<<scientific<<setprecision(14);
	for(ipls=0;ipls<npls;ipls++)
	{
		m=0;
		for(l=0;l<npls_main;l++)
		{
			m+=nGenByPls[l];
			if(ipls<m)
			{
				break;
			}
		}
		npc=0;
		for(j=l-1;j>=0;j--)
		{
			if(nGenByPls[j])
			{
				npc+=NrecB[j];
			}
		}
		for(i=RecvPlsIgB[ipls];i<RecvPlsIgB[ipls+1];i++)
		{
			Ic=IcB[i];
			cosphi=cosphiB[i];
			sinphi=-sinphiB[i];
			len=lenB[i];

			Tp=xyzVectorB[npc+i-RecvPlsIgB[ipls]];
			ofp<<(int)Tp.x<<" "<<(int)Tp.y<<endl;
			ofp<<Tp.x<<" "<<Tp.y<<" "<<Tp.z<<endl;
			ofp<<"1"<<endl;
			ofp<<"1  1  1  1"<<endl;
			ofp<<"1"<<endl;
			ofp<<"      t (s)         Bx (T)     By (T)     Bz (T)"<<endl;
			for(it=0;it<ntimesdec;it++)
			{
				const double& tj=times[it+tbeg];
				if (tj>=ftime&&tj<=ltime)
				{
					U[0]=U[1]=U[2]=0.0;
					for(idp=0;idp<NUMBEROFABDIPOLES;idp++)
					{
						k=i*NUMBEROFABDIPOLES+idp;
						U[0]+=Bx[it][k];
						U[1]+=By[it][k];
						U[2]+=Bz[it][k];
					}
					U[0]*=len;
					U[1]*=len;
					U[2]*=len;

					ux=U[0];
					uy=U[1];
					uz=U[2];

					x=ux*cosphi-uy*sinphi;
					y=ux*sinphi+uy*cosphi;

					ux=x;
					uy=y;

					ofp<<tj*1000.0<<"   "<<ux*1000.0<<" "<<uy*1000.0<<" "<<uz*1000.0<<'\n';
				}

				if(!it /* && !iDec */)
				{
					U[0]=U[1]=U[2]=0.0;
					for(idp=0;idp<NUMBEROFABDIPOLES;idp++)
					{
						k=i*NUMBEROFABDIPOLES+idp;
						U[0]+=Bx[it][k];
						U[1]+=By[it][k];
						U[2]+=Bz[it][k];
					}
					U[0]*=len;
					U[1]*=len;
					U[2]*=len;

					ux=U[0];
					uy=U[1];
					uz=U[2];

					x=ux*cosphi-uy*sinphi;
					y=ux*sinphi+uy*cosphi;

					ux=x;
					uy=y;

					ofpstat<<ux*1000.0<<" "<<uy*1000.0<<" "<<uz*1000.0<<'\n';
				}
			}
		}
	}
	ofp.close();
	ofp.clear();
	ofpstat.close();
	ofpstat.clear();

	sprintf_s(buf, "edsall.ax");
	ofp.open(buf);
	ofp<<scientific<<setprecision(14);
	for(ipls=0;ipls<npls;ipls++)
	{
		m=0;
		for(l=0;l<npls_main;l++)
		{
			m+=nGenByPls[l];
			if(ipls<m)
			{
				break;
			}
		}
		npc=0;
		for(j=l-1;j>=0;j--)
		{
			if(nGenByPls[j])
			{
				npc+=NrecB[j];
			}
		}
		
		for(i=RecvPlsIgB[ipls];i<RecvPlsIgB[ipls+1];i++)
		{
			Ic=IcB[i];
			cosphi=cosphiB[i];
			sinphi=-sinphiB[i];
			len=lenB[i];

			Tp=xyzVectorB[npc+i-RecvPlsIgB[ipls]];
			ofp<<(int)Tp.x<<" "<<(int)Tp.y<<endl;
			ofp<<Tp.x<<" "<<Tp.y<<" "<<Tp.z<<endl;
			ofp<<"1"<<endl;
			ofp<<"1  1  1  1"<<endl;
			ofp<<"1"<<endl;
			ofp<<"      t (s)         Emfx (T)     Emfy (T)     Emfz (T)"<<endl;
			for(it=0;it<ntimesdec;it++)
			{
				const double& tj=times[it+tbeg];
				if (tj>=ftime&&tj<=ltime)
				{
					U[0]=U[1]=U[2]=0.0;
					for(idp=0;idp<NUMBEROFABDIPOLES;idp++)
					{
						k=i*NUMBEROFABDIPOLES+idp;
						U[0]+=Fx[it][k];
						U[1]+=Fy[it][k];
						U[2]+=Fz[it][k];
					}
					U[0]*=len;
					U[1]*=len;
					U[2]*=len;

					ux=U[0];
					uy=U[1];
					uz=U[2];

					x=ux*cosphi-uy*sinphi;
					y=ux*sinphi+uy*cosphi;

					ux=x;
					uy=y;

					ofp<<tj*1000.0<<"   "<<ux*1000.0<<" "<<uy*1000.0<<" "<<uz*1000.0<<'\n';
				}
			}
		}
	}
	ofp.close();
	ofp.clear();


	sprintf_s(buf, "eall.ax");
	ofp.open(buf);
	ofp<<scientific<<setprecision(14);
	for(ipls=0;ipls<npls;ipls++)
	{
		m=0;
		for(l=0;l<npls_main;l++)
		{
			m+=nGenByPls[l];
			if(ipls<m)
			{
				break;
			}
		}
		npc=0;
		for(j=l-1;j>=0;j--)
		{
			if(nGenByPls[j])
			{
				npc+=NrecE[j];
			}
		}

		for(i=RecvPlsIgE[ipls];i<RecvPlsIgE[ipls+1];i++)
		{
			Ic=IcE[i];
			cosphi=cosphiE[i];
			sinphi=-sinphiE[i];
			len=lenE[i];

			Tp=xyzVectorE[npc+i-RecvPlsIgE[ipls]];
			ofp<<(int)Tp.x<<" "<<(int)Tp.y<<endl;
			ofp<<Tp.x<<" "<<Tp.y<<" "<<Tp.z<<endl;
			ofp<<"1"<<endl;
			ofp<<"1  1  1  1"<<endl;
			ofp<<"1"<<endl;
			ofp<<"      t (s)         Ex (V/m)     Ey (V/m)     Ez (V/m)"<<endl;
			for(it=0;it<ntimesdec;it++)
			{
				const double& tj=times[it+tbeg];
				if (tj>=ftime&&tj<=ltime)
				{
					U[0]=0.0;
					for(idp=0;idp<NUMBEROFABDIPOLES;idp++)
					{
						k=i*NUMBEROFABDIPOLES+idp;
						U[0]+=E[it][k];
					}
					U[0]*=len;

					ux=U[0]*Ic.x;
					uy=U[0]*Ic.y;
					uz=U[0]*Ic.z;

					x=ux*cosphi-uy*sinphi;
					y=ux*sinphi+uy*cosphi;
				
					ux=x;
					uy=y;

					ofp<<tj*1000.0<<"   "<<ux*1000.0<<" "<<uy*1000<<" "<<uz*1000.0<<'\n';
				}
			}
		}
	}
	ofp.close();
	ofp.clear();
}

void task2d_Ax::clear()
{
	ClearVector2d(Bx,nt);
	ClearVector2d(By,nt);
	ClearVector2d(Bz,nt);
	ClearVector2d(Fx,nt);
	ClearVector2d(Fy,nt);
	ClearVector2d(Fz,nt);
	ClearVector2d(E,nt);

	ClearVector(H);
	ClearVector(E0);
	ClearVector(recE);
	ClearVector(recE0);
	ClearVector(recB);
	ClearVector(sinfiE);
	ClearVector(sinfiE0);
	ClearVector(cosfiE);
	ClearVector(cosfiE0);
	ClearVector(sinfiB);
	ClearVector(cosfiB);
}


int task2d_Ar::init()
{
	PathTo2d[0]='\0';
	inf.open("PathTo2d");
	if(inf)
	{
		inf>>PathTo2d;
		logfile<<"Outputing from "<<PathTo2d<<'\n';
		inf.close();
	}
	else
	{
		strcpy(PathTo2d,".");
	}
	inf.clear();
	strcat(PathTo2d,"\\Ar");

	GetNumberOfPlaces(npls);

	nGenByPls.resize(npls);
	NrecB.resize(npls);
	NrecE.resize(npls);

	inf.open("srsclcgsz");
	if(!inf)
	{
		logfile<<"Error in open file "<<"srsclcgsz"<<endl;
		cout<<"Error in open file "<<"srsclcgsz"<<endl;
		return 1;
	}
	p1=0;
	for(i=0;i<npls;i++)
	{
		inf>>nGenByPls[i];
		p1+=nGenByPls[i];
	}
	inf.close();
	inf.clear();


	RecvPlsIgB.resize(p1+1);
	RecvPlsIgE.resize(p1+1);
	RecvPlsIgE0.resize(p1+1);

	inf.open("xyzVectorE0");
	if(inf)
	{
		inf>>npntE0;
		inf.close();
	}
	else
	{
		npntE0=0;
	}
	inf.clear();

	if(npntE0)
	{
		inf.open("TgCompE0");
		if(!inf)
		{
			logfile<<"Error in open file "<<"TgCompE0"<<endl;
			cout<<"Error in open file "<<"TgCompE0"<<endl;
			return 1;
		}
		TgCompFlagX.resize(npntE0);
		TgCompFlagY.resize(npntE0);
		TgCompFlagZ.resize(npntE0);
		TgCompE0.resize(npntE0);
		for(i=0;i<npntE0;i++)
		{
			inf>>TgCompE0[i].x>>TgCompE0[i].y>>TgCompE0[i].z;
			TgCompFlagX[i]=(fabs(TgCompE0[i].x)>d_eps_3);
			TgCompFlagY[i]=(fabs(TgCompE0[i].y)>d_eps_3);
			TgCompFlagZ[i]=(fabs(TgCompE0[i].z)>d_eps_3);
		}
		inf.close();
		inf.clear();
	}

	inf.open("recvsb");
	if(!inf)
	{
		logfile<<"Error in open file "<<"recvsb"<<endl;
		cout<<"Error in open file "<<"recvsb"<<endl;
		return 1;
	}
	RecvPlsIgB[0]=0;
	m=0;
	for(i=0;i<npls;i++)
	{
		inf>>NrecB[i];
		if(nGenByPls[i])
		{
			RecvPlsIgB[m+1]=NrecB[i];
			m++;
		}
	}
	inf.close();
	inf.clear();

	inf.open("recvse");
	if(!inf)
	{
		logfile<<"Error in open file "<<"recvse"<<endl;
		cout<<"Error in open file "<<"recvse"<<endl;
		return 1;
	}
	RecvPlsIgE[0]=0;
	m=0;
	for(i=0;i<npls;i++)
	{
		inf>>NrecE[i];
		if(nGenByPls[i])
		{
			RecvPlsIgE[m+1]=NrecE[i];
			m++;
		}
	}
	inf.close();
	inf.clear();


	m=0;
	RecvPlsIgE0[0]=0;
	for(i=0;i<npls;i++)
	{
		if(nGenByPls[i])
		{
			RecvPlsIgE0[m+1]=npntE0;
			m++;
		}
	}


	p2=p1;
	j=npls;
	for(i=m;i>0;i--)
	{
		while(j>0 && !nGenByPls[j-1]){j--;}
		RecvPlsIgB[p2]=RecvPlsIgB[i];
		for(k=p2-1;k>p2-nGenByPls[j-1];k--){RecvPlsIgB[k]=RecvPlsIgB[p2];}
		p2-=nGenByPls[j-1];
		j--;
	}


	p2=p1;
	j=npls;
	for(i=m;i>0;i--)
	{
		while(j>0 && !nGenByPls[j-1]){j--;}
		RecvPlsIgE[p2]=RecvPlsIgE[i];
		for(k=p2-1;k>p2-nGenByPls[j-1];k--){RecvPlsIgE[k]=RecvPlsIgE[p2];}
		p2-=nGenByPls[j-1];
		j--;
	}

	p2=p1;
	j=npls;
	for(i=m;i>0;i--)
	{
		while(j>0 && !nGenByPls[j-1]){j--;}
		RecvPlsIgE0[p2]=RecvPlsIgE0[i];
		for(k=p2-1;k>p2-nGenByPls[j-1];k--){RecvPlsIgE0[k]=RecvPlsIgE0[p2];}
		p2-=nGenByPls[j-1];
		j--;
	}

	npls_main=npls;
	npls=p1;

	for(i=0;i<npls;i++)
	{
		RecvPlsIgB[i+1]+=RecvPlsIgB[i];
		RecvPlsIgE[i+1]+=RecvPlsIgE[i];
		RecvPlsIgE0[i+1]+=RecvPlsIgE0[i];
	}

	npntB=RecvPlsIgB[npls];
	npntE=RecvPlsIgE[npls];
	npntE0=RecvPlsIgE0[npls];

	RecToSourceB.resize(npntB);
	RecToSourceE.resize(npntE);
	RecToSourceE0.resize(npntE0);

	for(i=0;i<npls;i++)
	{
		for(k=RecvPlsIgB[i];k<RecvPlsIgB[i+1];k++)
		{
			RecToSourceB[k]=i;
		}
	}

	for(i=0;i<npls;i++)
	{
		for(k=RecvPlsIgE[i];k<RecvPlsIgE[i+1];k++)
		{
			RecToSourceE[k]=i;
		}
	}

	for(i=0;i<npls;i++)
	{
		for(k=RecvPlsIgE0[i];k<RecvPlsIgE0[i+1];k++)
		{
			RecToSourceE0[k]=i;
		}
	}

	TgCompFlagXAll.resize(npntE0);
	TgCompFlagYAll.resize(npntE0);
	TgCompFlagZAll.resize(npntE0);
	TgCompFlagRAll.resize(npntE0);
	for(i=0;i<npls;i++)
	{
		j=0;
		for(k=RecvPlsIgE0[i];k<RecvPlsIgE0[i+1];k++)
		{
			TgCompFlagXAll[k]=TgCompFlagX[j];
			TgCompFlagYAll[k]=TgCompFlagY[j];
			TgCompFlagZAll[k]=TgCompFlagZ[j];
			TgCompFlagRAll[k]=(TgCompFlagX[j] || TgCompFlagY[j]);
			j++;
		}
	}
	
	GenSq.resize(npls);
	inf.open("srsclca");
	if(!inf)
	{
		logfile<<"Error in open file "<<"srsclca"<<endl;
		cout<<"Error in open file "<<"srsclca"<<endl;
		return 1;
	}
	for(i=0;i<npls;i++){inf>>GenSq[i].A.x>>GenSq[i].A.y>>GenSq[i].A.z>>GenSq[i].B.x>>GenSq[i].B.y>>GenSq[i].B.z;}
	inf.close();
	inf.clear();

	inf.open("timeintervalforprint");
	if(!inf)
	{
		logfile<<"No file "<<"timeintervalforprint"<<'\n';
		return 1;
	}
	inf>>ftime>>ltime;
	inf.close();
	inf.clear();

	time_gaps_for_loop=0;




	SetPath(PathTo2d);

	retp=Read(0);
	if(retp)
	{
		logfile<<"Function Read() for Ar returned "<<retp<<'\n';
		return 1;
	}

	retp=ReadMtr3d2d(".");
	if(retp)
	{
		logfile<<"Function ReadMtr3d2d() for . returned "<<retp<<'\n';
		return 1;
	}

	rc0=0.5*(rm[0]+rm[1]);

	DminSrc=rc0;
	inf.open("dminsrc");
	if(inf)
	{
		inf>>DminSrc;
		inf.close();
	}
	inf.clear();

	inf.open("xyzVectorB");
	if(!inf)
	{
		logfile<<"No file "<<"xyzVectorB"<<'\n';
		return 1;
	}
	inf>>npntBor;
	xyzVectorB=new PointXYZ[npntBor];
	recB=new PointRZ[npntB*2];
	sinfiB=new double[npntB*2];
	cosfiB=new double[npntB*2];
	npc=0;
	l=0;
	for(ipm=0;ipm<npls_main;ipm++)
	{
		if(!nGenByPls[ipm])
		{
			for(i=0;i<NrecB[ipm];i++)
			{
				inf>>x>>y>>z;
			}
		}
		m=l;
		for(ipls=0;ipls<nGenByPls[ipm];ipls++)
		{
			for(i=RecvPlsIgB[npc+ipls];i<RecvPlsIgB[npc+ipls+1];i++)
			{
				Av=GenSq[RecToSourceB[i]].A;
				Bv=GenSq[RecToSourceB[i]].B;
				if(!ipls)
				{
					inf>>Tp.x>>Tp.y>>Tp.z;
					xyzVectorB[l]=Tp;
					l++;
				}
				else
				{
					Tp=xyzVectorB[m+i-RecvPlsIgB[npc+ipls]];
				}
				AddReciver(Tp,Av,recB[2*i],sinfiB[2*i],cosfiB[2*i],rc0);
				AddReciver(Tp,Bv,recB[2*i+1],sinfiB[2*i+1],cosfiB[2*i+1],rc0);
			}
		}
		npc+=nGenByPls[ipm];
	}
	inf.close();
	inf.clear();


	inf.open("xyzVectorE");
	if(!inf)
	{
		logfile<<"No file "<<"xyzVectorE"<<'\n';
		return 1;
	}
	inf>>npntEor;
	xyzVectorE=new PointXYZ[npntEor];
	recE=new PointRZ[npntE*2];
	sinfiE=new double[npntE*2];
	cosfiE=new double[npntE*2];
	matrec=new int[npntE*2];
	npc=0;
	l=0;
	for(ipm=0;ipm<npls_main;ipm++)
	{
		if(!nGenByPls[ipm])
		{
			for(i=0;i<NrecE[ipm];i++)
			{
				inf>>x>>y>>z;
			}
		}
		m=l;
		for(ipls=0;ipls<nGenByPls[ipm];ipls++)
		{
			for(i=RecvPlsIgE[npc+ipls];i<RecvPlsIgE[npc+ipls+1];i++)
			{
				Av=GenSq[RecToSourceE[i]].A;
				Bv=GenSq[RecToSourceE[i]].B;
				if(!ipls)
				{
					inf>>Tp.x>>Tp.y>>Tp.z;
					xyzVectorE[l]=Tp;
					l++;
				}
				else
				{
					Tp=xyzVectorE[m+i-RecvPlsIgE[npc+ipls]];
				}
				AddReciver(Tp,Av,recE[2*i],sinfiE[2*i],cosfiE[2*i],rc0);
				AddReciver(Tp,Bv,recE[2*i+1],sinfiE[2*i+1],cosfiE[2*i+1],rc0);
				matrec[2*i]=0;
				matrec[2*i+1]=0;
			}
		}
		npc+=nGenByPls[ipm];
	}
	inf.close();
	inf.clear();

	npntE0or=0;
	if(npntE0)
	{
		inf.open("xyzVectorE0");
		if(!inf)
		{
			logfile<<"No file "<<"xyzVectorE0"<<'\n';
			return 1;
		}
		inf>>npntE0or;
		xyzVectorE0=new PointXYZ[npntE0or];
		for(i=0;i<npntE0or;i++)
		{
			inf>>Tp.x>>Tp.y>>Tp.z;
			xyzVectorE0[i]=Tp;
		}
		inf.close();
		inf.clear();
	}

	recE0=new PointRZ[npntE0*2];
	sinfiE0=new double[npntE0*2];
	cosfiE0=new double[npntE0*2];
	matrec0=new int[npntE0*2];
	npc=0;
	for(ipls=0;ipls<npls;ipls++)
	{
		l=0;
		for(i=RecvPlsIgE0[ipls];i<RecvPlsIgE0[ipls+1];i++)
		{
			Av=GenSq[RecToSourceE0[i]].A;
			Bv=GenSq[RecToSourceE0[i]].B;
			Tp=xyzVectorE0[l];
			l++;
			AddReciver(Tp,Av,recE0[2*i],sinfiE0[2*i],cosfiE0[2*i],rc0);
			AddReciver(Tp,Bv,recE0[2*i+1],sinfiE0[2*i+1],cosfiE0[2*i+1],rc0);
			matrec0[2*i]=0;
			matrec0[2*i+1]=0;
		}
	}


	nt=ntimesdec;

	A=new double *[3];
	H=new double[kpnt*npls];
	B=new double*[nt];
	F=new double *[nt];
	E=new double*[nt];
	E0=new double[2*npntE0];

	for(it=0;it<nt;it++)
	{
		B[it]=new double[2*npntB];
		F[it]=new double[2*npntB];
		E[it]=new double[2*npntE];
	}

	retp=SobjB.Init(npntB*2,kpnt,krect,recB,pnt,rect,0,NULL,sigma,
		sigmaZ,mu,nreg,qr,qz,reg,rm,zm,0);
	if(retp)
	{
		logfile<<"SobjB.Init() returned "<<retp<<endl;
		return retp;
	}

	retp=SobjE.Init(npntE*2,kpnt,krect,recE,pnt,rect,0,NULL,sigma,
		sigmaZ,mu,nreg,qr,qz,reg,rm,zm,0);
	if(retp)
	{
		logfile<<"SobjE.Init() returned "<<retp<<endl;
		return retp;
	}


	double rmin,zmin,rmax,zmax;

	vRecOutData = new RectOutData[npntE0*2];
	if(!vRecOutData)return RETCODE_NOMEM;

	rmin = zmin = 1e+30;
	rmax = zmax = -1e+30;

	for(i=0;i<npntE0*2;i++)
	{
		if(recE0[i].r<rmin){rmin=recE0[i].r; }
		if(recE0[i].r>rmax){rmax=recE0[i].r; }
		if(recE0[i].z<zmin){zmin=recE0[i].z; }
		if(recE0[i].z>zmax){zmax=recE0[i].z; }
	}

	retp=FindElementsForReceivers(npntE0*2,qr,qz,rm,zm,reg,recE0,vRecOutData,pnt,rect,sigma,sigmaZ,rmin,zmin,rmax,zmax);
	if (retp)
	{
		logfile << "FindElementsForReceivers returned " << retp << endl;
		return retp;
	}


	for(it=0;it<3;it++)
	{
		A[2-it]=new double[kpnt*npls];
		if(!A[2-it])return RETCODE_NOMEM;
		GetField(it,A[2-it],str);
	}

	return 0;
}

void task2d_Ar::output(int it,vector<float> &enor)
{
	if(tbeg+it>2)
	{
		RollFields(A[0],A[1],A[2]);
		GetField(it,A[0],str);
	}

	double *A0=(tbeg+it>=2)? A[0] : (tbeg+it==1)? A[1] : A[2];

	SobjB.Output_dArdz(A0,B[it],1,RecToSourceB);

	if(tbeg+it>=1)
	{
		_diff_t2(F[it],B[it],B[it-1],2*npntB,times,it);
	}
	else
	{
		for(i=0;i<2*npntB;i++){F[it][i]=0.0;}
	}

	for(i=0;i<2*npntE0;i++){E0[i]=0.0;}

	if(tbeg+it>=1)
	{
		_diff_t2(H,A[0],A[1],kpnt*npls,times+tbeg,it);
	}
	else
	{
		for(i=0;i<2*npntE;i++){E[it][i]=0.0;}
		for(i=0;i<2*npntE0;i++){E0[i]=0.0;}
	}
	if(tbeg+it>=1)
	{
		SobjE.Output_Er(H,E[it],1,RecToSourceE);
		Output_FieldAr_Pr(npntE0*2,H,E0,vRecOutData,RecToSourceE0,kpnt);
	}

	if(npntE0or)
	{
		k=0;
		for(ipls=0;ipls<npls;ipls++)
		{
			j=0;
			for(i=RecvPlsIgE0[ipls];i<RecvPlsIgE0[ipls+1];i++)
			{
				if(recE0[i].r>DminSrc || (fabs(recE0[i].z-GenSq[ipls].A.z)>DminSrc && fabs(recE0[i].z-GenSq[ipls].B.z)>DminSrc))
				{
					x=(E0[2*i+1]*cosfiE0[2*i+1]-E0[2*i]*cosfiE0[2*i]);
					y=(E0[2*i+1]*sinfiE0[2*i+1]-E0[2*i]*sinfiE0[2*i]);
					enor[k]=float(x*TgCompE0[j].x+y*TgCompE0[j].y);
				}
				else
				{
					enor[k]=0.0;
				}
				j++;
				k++;
			}
		}

	}
}

void task2d_Ar::output_rec()
{
	sprintf_s(buf, "ball.ar.stat");
	ofpstat.open(buf);
	sprintf_s(buf, "ball.ar");
	ofp.open(buf);
	ofp<<scientific<<setprecision(14);
	for(ipls=0;ipls<npls;ipls++)
	{
		m=0;
		for(l=0;l<npls_main;l++)
		{
			m+=nGenByPls[l];
			if(ipls<m)
			{
				break;
			}
		}
		npc=0;
		for(j=l-1;j>=0;j--)
		{
			if(nGenByPls[j])
			{
				npc+=NrecB[j];
			}
		}
		for(i=RecvPlsIgB[ipls];i<RecvPlsIgB[ipls+1];i++)
		{
			Tp=xyzVectorB[npc+i-RecvPlsIgB[ipls]];
			ofp<<(int)Tp.x<<" "<<(int)Tp.y<<endl;
			ofp<<Tp.x<<" "<<Tp.y<<" "<<Tp.z<<endl;
			ofp<<"1"<<endl;
			ofp<<"1  1  1  1"<<endl;
			ofp<<"1"<<endl;
			ofp<<"      t (s)         Bx (T)     By (T)     Bz (T)"<<endl;
			for(it=0;it<ntimesdec;it++)
			{
				const double& tj=times[it+tbeg];
				if (tj>=ftime&&tj<=ltime)
				{
					x=(-B[it][2*i+1]*sinfiB[2*i+1]+B[it][2*i]*sinfiB[2*i]);
					y=(B[it][2*i+1]*cosfiB[2*i+1]-B[it][2*i]*cosfiB[2*i]);
					ofp<<tj*1000<<"   "<<x*1000<<" "<<y*1000<<" "<<0.0<<'\n';
				}

				if(!it /* && !iDec */)
				{
					x=(-B[it][2*i+1]*sinfiB[2*i+1]+B[it][2*i]*sinfiB[2*i]);
					y=(B[it][2*i+1]*cosfiB[2*i+1]-B[it][2*i]*cosfiB[2*i]);
					ofpstat<<x*1000<<" "<<y*1000<<" "<<0.0<<'\n';
				}
			}
		}
	}
	ofp.close();
	ofp.clear();
	ofpstat.close();
	ofpstat.clear();
	
	sprintf_s(buf, "edsall.ar");
	ofp.open(buf);
	ofp<<scientific<<setprecision(14);
	for(ipls=0;ipls<npls;ipls++)
	{
		m=0;
		for(l=0;l<npls_main;l++)
		{
			m+=nGenByPls[l];
			if(ipls<m)
			{
				break;
			}
		}
		npc=0;
		for(j=l-1;j>=0;j--)
		{
			if(nGenByPls[j])
			{
				npc+=NrecB[j];
			}
		}
		for(i=RecvPlsIgB[ipls];i<RecvPlsIgB[ipls+1];i++)
		{
			Tp=xyzVectorB[npc+i-RecvPlsIgB[ipls]];
			ofp<<(int)Tp.x<<" "<<(int)Tp.y<<endl;
			ofp<<Tp.x<<" "<<Tp.y<<" "<<Tp.z<<endl;
			ofp<<"1"<<endl;
			ofp<<"1  1  1  1"<<endl;
			ofp<<"1"<<endl;
			ofp<<"      t (s)         Emfx (T)     Emfy (T)     Emfz (T)"<<endl;
			for(it=0;it<ntimesdec;it++)
			{
				const double& tj=times[it+tbeg];
				if (tj>=ftime&&tj<=ltime)
				{
					x=(-F[it][2*i+1]*sinfiB[2*i+1]+F[it][2*i]*sinfiB[2*i]);
					y=(F[it][2*i+1]*cosfiB[2*i+1]-F[it][2*i]*cosfiB[2*i]);
					ofp<<tj*1000<<"   "<<x*1000<<" "<<y*1000<<" "<<0.0<<'\n';
				}
			}
		}
	}
	ofp.close();
	ofp.clear();


	sprintf_s(buf, "eall.ar");
	ofp.open(buf);
	ofp<<scientific<<setprecision(14);

	for(ipls=0;ipls<npls;ipls++)
	{
		m=0;
		for(l=0;l<npls_main;l++)
		{
			m+=nGenByPls[l];
			if(ipls<m)
			{
				break;
			}
		}
		npc=0;
		for(j=l-1;j>=0;j--)
		{
			if(nGenByPls[j])
			{
				npc+=NrecE[j];
			}
		}

		for(i=RecvPlsIgE[ipls];i<RecvPlsIgE[ipls+1];i++)
		{
			Tp=xyzVectorE[npc+i-RecvPlsIgE[ipls]];

			ofp<<(int)Tp.x<<" "<<(int)Tp.y<<endl;
			ofp<<Tp.x<<" "<<Tp.y<<" "<<Tp.z<<endl;
			ofp<<"1"<<endl;
			ofp<<"1  1  1  1"<<endl;
			ofp<<"1"<<endl;
			ofp<<"      t (s)         Ex (V/m)     Ey (V/m)     Ez (V/m)"<<endl;
			for(it=0;it<ntimesdec;it++)
			{
				const double& tj=times[it+tbeg];
				if (tj>=ftime&&tj<=ltime)
				{
					x=(E[it][2*i+1]*cosfiE[2*i+1]-E[it][2*i]*cosfiE[2*i]);
					y=(E[it][2*i+1]*sinfiE[2*i+1]-E[it][2*i]*sinfiE[2*i]);
					ofp<<tj*1000<<"   "<<x*1000<<" "<<y*1000<<" "<<0.0<<'\n';
				}
			}
		}
	}
	ofp.close();
	ofp.clear();
}

void task2d_Ar::clear()
{
	if(A)
	{
		for(it=0;it<3;it++)
		{
			if(A[it]){delete [] A[it]; A[it]=NULL;}
		}
		delete [] A; A=NULL;
	}
	if(B)
	{
		for(it=0;it<nt;it++)
		{
			if(B[it]){delete [] B[it]; B[it]=NULL;}
		}
		delete [] B; B=NULL;
	}
	if(E)
	{
		for(it=0;it<nt;it++)
		{
			if(E[it]){delete [] E[it]; E[it]=NULL;}
		}
		delete [] E; E=NULL;
	}

	if(H){delete [] H; H=NULL;}
	if(recE){delete [] recE; recE=NULL;}
	if(recB){delete [] recB; recB=NULL;}
	if(sinfiE){delete [] sinfiE; sinfiE=NULL;}
	if(cosfiE){delete [] cosfiE; cosfiE=NULL;}
	if(sinfiB){delete [] sinfiB; sinfiB=NULL;}
	if(cosfiB){delete [] cosfiB; cosfiB=NULL;}


	if(F)
	{
		for(it=0;it<nt;it++)
		{
			if(F[it]){delete [] F[it]; F[it]=NULL;}
		}
		delete [] F; F=NULL;
	}
}


int task2d_Hfi::ReadDiscontinues(int npls)
{
	int i,ipls;
	ifstream infn,infs;

	kt1dc.resize(npls);
	l1dc.resize(npls);

	sprintf(buf,"%s\\kt1dc",path);
	infs.open(buf);
	sprintf(buf,"%s\\l1dc.dat",path);
	infn.open(buf,ios::binary);

	if(!infs)
	{
		cout<<"Error in open file "<<"kt1dc"<<endl;
		logfile<<"Error in open file "<<"kt1dc"<<endl;
		return RETCODE_NOFILE;
	}
	if(!infn)
	{
		cout<<"Error in open file "<<"l1dc.dat"<<endl;
		logfile<<"Error in open file "<<"l1dc.dat"<<endl;
		return RETCODE_NOFILE;
	}

	for(ipls=0;ipls<npls;ipls++)
	{
		infs>>kt1dc[ipls];
		l1dc[ipls].resize(kt1dc[ipls]);
		for(i=0;i<kt1dc[ipls];i++){infn>l1dc[ipls][i];}
	}

	infs.close();
	infs.clear();
	infn.close();
	infn.clear();

	return RETCODE_OK;
}

void task2d_Hfi::RepairGeneratorOrder(vector<int> &GenToPls,vector<int> &nGels,vector<int> &nVels,vector<SqLoop> &GenSq,int n,int npls)
{
	int i,l,m,r,ipls;
	vector<SqLoop> GenSqT=GenSq;

	l=0;
	for(ipls=0;ipls<npls;ipls++)
	{
		m=0;
		for(r=ipls-1;r>=0;r--)
		{
			m+=nGels[r]+nVels[r];
		}
		for(i=0;i<nGels[ipls];i++)
		{
			GenSq[m+i]=GenSqT[l];
			l++;
		}
	}
	for(ipls=0;ipls<npls;ipls++)
	{
		m=0;
		for(r=ipls-1;r>=0;r--)
		{
			m+=nGels[r]+nVels[r];
		}
		for(i=0;i<nVels[ipls];i++)
		{
			GenSq[nGels[ipls]+m+i]=GenSqT[l];
			l++;
		}
	}
}

void task2d_Hfi::RepairFieldOrder(vector<int> &GenToPls,vector<int> &nGels,vector<int> &nVels,double *(&Hphi),double *(&THphi),int krect,int n,int npls)
{
	int i,j,k,l,m,r,ipls,t,p,tt,pp;
	double *tu;

	tu=THphi;
	THphi=Hphi;
	Hphi=tu;

	l=0;
	for(ipls=0;ipls<npls;ipls++)
	{
		m=0;
		for(r=ipls-1;r>=0;r--)
		{
			m+=nGels[r]+nVels[r];
		}
		for(i=0;i<nGels[ipls];i++)
		{
			tt=(m+i)*4*krect;
			pp=l*4*krect;
			for(k=0;k<krect;k++)
			{
				t=tt+4*k;
				p=pp+4*k;
				for(j=0;j<4;j++)
				{
					Hphi[t+j]=THphi[p+j];
				}
			}
			l++;
		}
	}
	for(ipls=0;ipls<npls;ipls++)
	{
		m=0;
		for(r=ipls-1;r>=0;r--)
		{
			m+=nGels[r]+nVels[r];
		}
		for(i=0;i<nVels[ipls];i++)
		{
			tt=(nGels[ipls]+m+i)*4*krect;
			pp=l*4*krect;
			for(k=0;k<krect;k++)
			{
				t=tt+4*k;
				p=pp+4*k;
				for(j=0;j<4;j++)
				{
					Hphi[t+j]=THphi[p+j];
				}
			}
			l++;
		}
	}
}



int task2d_Hfi::GetField(int it,double *v2,double *(&Hphi),char *str,int *DsFlgNode,
			 double *(&THphi),vector<int> &nGels,vector<int> &nVels,vector<int> &GenToPls,int npls_main)
{
	int i,j,k,ipls;
	FILE *fp;

	sprintf(str,"%s\\v2.%d",path,it+tbeg);
	fp=fopen(str,"rb");
	if(!fp) return RETCODE_NOFILE;
	fread(v2,sizeof(double),nc*npls,fp);
	fclose(fp);
	
	for(ipls=0;ipls<npls;ipls++)
	{
		for(i=0;i<kpnt;i++){DsFlgNode[i]=0;}
		for(k=0;k<kt1dc[ipls];k++){DsFlgNode[l1dc[ipls][k]-1]=1;}

		for(i=0;i<krect;i++)
		{
			for(j=0;j<4;j++)
			{
				Hphi[4*i+j+4*krect*ipls]=v2[rect[i].nodes[j]-1+kpnt*ipls];
			}
			

			if(DsFlgNode[rect[i].nodes[2]-1] && DsFlgNode[rect[i].nodes[3]-1])
			{
				for(j=2;j<=3;j++)
				{
					Hphi[4*i+j+4*krect*ipls] += curr[it+tbeg]/(2.0*PI*pnt[rect[i].nodes[j]-1].r);
				}
			}

			if(DsFlgNode[rect[i].nodes[0]-1] && DsFlgNode[rect[i].nodes[2]-1])
			{
				for(j=0;j<=2;j+=2)
				{
					Hphi[4*i+j+4*krect*ipls] += curr[it+tbeg]/(2.0*PI*pnt[rect[i].nodes[j]-1].r);
				}
			}
		}
	}


	RepairFieldOrder(GenToPls,nGels,nVels,Hphi,THphi,krect,npls,npls_main);

	return 0;
}

int task2d_Hfi::init()
{
	double tmpd;

	PathTo2d[0]='\0';
	inf.open("PathTo2d");
	if(inf)
	{
		inf>>PathTo2d;
		logfile<<"Outputing from "<<PathTo2d<<'\n';
		inf.close();
	}
	else
	{
		strcpy(PathTo2d,".");
	}
	inf.clear();
	strcat(PathTo2d,"\\Hfi");

	GetNumberOfPlaces(npls);

	nGenByPls.resize(npls);
	NrecB.resize(npls);
	NrecE.resize(npls);

	nGels.resize(npls);
	inf.open("srsclcgsz");
	if(!inf)
	{
		logfile<<"Error in open file "<<"srsclcgsz"<<endl;
		cout<<"Error in open file "<<"srsclcgsz"<<endl;
		return 1;
	}
	p1=0;
	for(i=0;i<npls;i++)
	{
		inf>>nGels[i];
		nGenByPls[i]=nGels[i];
	}
	inf.close();
	inf.clear();

	nGelsAll=0;
	for(i=0;i<npls;i++){nGelsAll+=nGels[i];}

	nVels.resize(npls);
	inf.open("srsclcvsz");
	if(!inf)
	{
		logfile<<"Error in open file "<<"srsclcvsz"<<endl;
		cout<<"Error in open file "<<"srsclcvsz"<<endl;
		return 1;
	}
	for(i=0;i<npls;i++)
	{
		inf>>nVels[i];
		nGenByPls[i]+=nVels[i];
	}
	inf.close();
	inf.clear();

	p1=0;
	for(i=0;i<npls;i++)
	{
		p1+=nGenByPls[i];
	}

	GenToPls.resize(p1);
	GenType.resize(p1);
	GenTypeTrue.resize(p1);

	k=0;
	for(i=0;i<npls;i++)
	{
		for(j=0;j<nGels[i];j++)
		{
			GenTypeTrue[k]=0;
			k++;
		}
		for(j=0;j<nVels[i];j++)
		{
			GenTypeTrue[k]=1;
			k++;
		}
	}

	for(i=0;i<p1;i++)
	{
		if(i<nGelsAll)
		{
			k=0;
			for(ipls=0;ipls<npls;ipls++)
			{
				if(i<k+nGels[ipls])
				{
					GenToPls[i]=ipls;
					GenType[i]=0;
					break;
				}
				k+=nGels[ipls];
			}
		}
		else
		{
			k=nGelsAll;
			for(ipls=0;ipls<npls;ipls++)
			{
				if(i<k+nVels[ipls])
				{
					GenToPls[i]=ipls;
					GenType[i]=1;
					break;
				}
				k+=nVels[ipls];
			}
		}
	}


	RecvPlsIgB.resize(p1+1);
	RecvPlsIgE.resize(p1+1);
	RecvPlsIgE0.resize(p1+1);

	inf.open("xyzVectorE0");
	if(inf)
	{
		inf>>npntE0;
		inf.close();
	}
	else
	{
		npntE0=0;
	}
	inf.clear();

	if(npntE0)
	{
		inf.open("TgCompE0");
		if(!inf)
		{
			logfile<<"Error in open file "<<"TgCompE0"<<endl;
			cout<<"Error in open file "<<"TgCompE0"<<endl;
			return 1;
		}
		TgCompFlagX.resize(npntE0);
		TgCompFlagY.resize(npntE0);
		TgCompFlagZ.resize(npntE0);
		TgCompE0.resize(npntE0);
		for(i=0;i<npntE0;i++)
		{
			inf>>TgCompE0[i].x>>TgCompE0[i].y>>TgCompE0[i].z;
			TgCompFlagX[i]=(fabs(TgCompE0[i].x)>d_eps_3);
			TgCompFlagY[i]=(fabs(TgCompE0[i].y)>d_eps_3);
			TgCompFlagZ[i]=(fabs(TgCompE0[i].z)>d_eps_3);
		}
		inf.close();
		inf.clear();
	}

	inf.open("recvsb");
	if(!inf)
	{
		logfile<<"Error in open file "<<"recvsb"<<endl;
		cout<<"Error in open file "<<"recvsb"<<endl;
		return 1;
	}
	RecvPlsIgB[0]=0;
	for(i=0;i<npls;i++)
	{
		inf>>NrecB[i];
		RecvPlsIgB[i+1]=NrecB[i];
	}
	inf.close();
	inf.clear();

	inf.open("recvse");
	if(!inf)
	{
		logfile<<"Error in open file "<<"recvse"<<endl;
		cout<<"Error in open file "<<"recvse"<<endl;
		return 1;
	}
	RecvPlsIgE[0]=0;
	for(i=0;i<npls;i++)
	{
		inf>>NrecE[i];
		RecvPlsIgE[i+1]=NrecE[i];
	}
	inf.close();
	inf.clear();


	RecvPlsIgE0[0]=0;
	for(i=0;i<npls;i++){RecvPlsIgE0[i+1]=npntE0;}


	p2=p1;
	for(i=npls;i>0;i--)
	{
		RecvPlsIgB[p2]=RecvPlsIgB[i];
		for(k=p2-1;k>p2-nGenByPls[i-1];k--){RecvPlsIgB[k]=RecvPlsIgB[p2];}
		p2-=nGenByPls[i-1];
	}

	p2=p1;
	for(i=npls;i>0;i--)
	{
		RecvPlsIgE[p2]=RecvPlsIgE[i];
		for(k=p2-1;k>p2-nGenByPls[i-1];k--){RecvPlsIgE[k]=RecvPlsIgE[p2];}
		p2-=nGenByPls[i-1];
	}

	p2=p1;
	for(i=npls;i>0;i--)
	{
		RecvPlsIgE0[p2]=RecvPlsIgE0[i];
		for(k=p2-1;k>p2-nGenByPls[i-1];k--){RecvPlsIgE0[k]=RecvPlsIgE0[p2];}
		p2-=nGenByPls[i-1];
	}

	npls_main=npls;
	npls=p1;

	for(i=0;i<npls;i++)
	{
		RecvPlsIgB[i+1]+=RecvPlsIgB[i];
		RecvPlsIgE[i+1]+=RecvPlsIgE[i];
		RecvPlsIgE0[i+1]+=RecvPlsIgE0[i];
	}

	npntB=RecvPlsIgB[npls];
	npntE=RecvPlsIgE[npls];
	npntE0=RecvPlsIgE0[npls];

	RecToSourceB.resize(npntB);
	RecToSourceE.resize(npntE);
	RecToSourceE0.resize(npntE0);

	for(i=0;i<npls;i++)
	{
		for(k=RecvPlsIgB[i];k<RecvPlsIgB[i+1];k++)
		{
			RecToSourceB[k]=i;
		}
	}

	for(i=0;i<npls;i++)
	{
		for(k=RecvPlsIgE[i];k<RecvPlsIgE[i+1];k++)
		{
			RecToSourceE[k]=i;
		}
	}

	for(i=0;i<npls;i++)
	{
		for(k=RecvPlsIgE0[i];k<RecvPlsIgE0[i+1];k++)
		{
			RecToSourceE0[k]=i;
		}
	}

	TgCompFlagXAll.resize(npntE0);
	TgCompFlagYAll.resize(npntE0);
	TgCompFlagZAll.resize(npntE0);
	TgCompFlagRAll.resize(npntE0);
	for(i=0;i<npls;i++)
	{
		j=0;
		for(k=RecvPlsIgE0[i];k<RecvPlsIgE0[i+1];k++)
		{
			TgCompFlagXAll[k]=TgCompFlagX[j];
			TgCompFlagYAll[k]=TgCompFlagY[j];
			TgCompFlagZAll[k]=TgCompFlagZ[j];
			TgCompFlagRAll[k]=(TgCompFlagX[j] || TgCompFlagY[j]);
			j++;
		}
	}

	GenSq.resize(npls);
	inf.open("srsclch");
	if(!inf)
	{
		logfile<<"Error in open file "<<"srsclch"<<endl;
		cout<<"Error in open file "<<"srsclch"<<endl;
		return 1;
	}
	for(i=0;i<npls;i++){inf>>GenSq[i].A.x>>GenSq[i].A.y>>GenSq[i].A.z>>GenSq[i].B.x>>GenSq[i].B.y>>GenSq[i].B.z;}
	inf.close();
	inf.clear();

	RepairGeneratorOrder(GenToPls,nGels,nVels,GenSq,npls,npls_main);

	inf.open("timeintervalforprint");
	if(!inf)
	{
		logfile<<"No file "<<"timeintervalforprint"<<'\n';
		return 1;
	}
	inf>>ftime>>ltime;
	inf.close();
	inf.clear();

	time_gaps_for_loop=0;




	SetPath(PathTo2d);
	
	retp=Read(1);
	if(retp)
	{
		logfile<<"Function Read() for Hfi returned "<<retp<<'\n';
		return 1;
	}

	retp=ReadDiscontinues(npls);
	if(retp)
	{
		logfile<<"Function ReadDiscontinues() for Hfi returned "<<retp<<'\n';
		return 1;
	}

	retp=ReadMtr3d2d(".");
	if(retp)
	{
		logfile<<"Function ReadMtr3d2d() for . returned "<<retp<<'\n';
		return 1;
	}

	rc0=0.5*(rm[0]+rm[1]);

	DminSrc=rc0;
	inf.open("dminsrc");
	if(inf)
	{
		inf>>DminSrc;
		inf.close();
	}
	inf.clear();

	inf.open("xyzVectorB");
	if(!inf)
	{
		logfile<<"No file "<<"xyzVectorB"<<'\n';
		return 1;
	}
	inf>>npntBor;
	xyzVectorB=new PointXYZ[npntBor];
	recB=new PointRZ[npntB*2];
	sinfiB=new double[npntB*2];
	cosfiB=new double[npntB*2];
	npc=0;
	l=0;
	for(ipm=0;ipm<npls_main;ipm++)
	{
		m=l;
		for(ipls=0;ipls<nGenByPls[ipm];ipls++)
		{
			for(i=RecvPlsIgB[npc+ipls];i<RecvPlsIgB[npc+ipls+1];i++)
			{
				Av=GenSq[RecToSourceB[i]].A;
				Bv=GenSq[RecToSourceB[i]].B;
				if(!ipls)
				{
					inf>>Tp.x>>Tp.y>>Tp.z;
					xyzVectorB[l]=Tp;
					l++;
				}
				else
				{
					Tp=xyzVectorB[m+i-RecvPlsIgB[npc+ipls]];
				}
				AddReciver(Tp,Av,recB[2*i],sinfiB[2*i],cosfiB[2*i],rc0);
				AddReciver(Tp,Bv,recB[2*i+1],sinfiB[2*i+1],cosfiB[2*i+1],rc0);
			}
		}
		npc+=nGenByPls[ipm];
	}
	inf.close();
	inf.clear();
	

	inf.open("xyzVectorE");
	if(!inf)
	{
		logfile<<"No file "<<"xyzVectorE"<<'\n';
		return 1;
	}
	inf>>npntEor;
	xyzVectorE=new PointXYZ[npntEor];
	recE=new PointRZ[npntE*2];
	sinfiE=new double[npntE*2];
	cosfiE=new double[npntE*2];
	matrec=new int[npntE*2];
	npc=0;
	l=0;
	for(ipm=0;ipm<npls_main;ipm++)
	{
		m=l;
		for(ipls=0;ipls<nGenByPls[ipm];ipls++)
		{
			for(i=RecvPlsIgE[npc+ipls];i<RecvPlsIgE[npc+ipls+1];i++)
			{
				Av=GenSq[RecToSourceE[i]].A;
				Bv=GenSq[RecToSourceE[i]].B;
				if(!ipls)
				{
					inf>>Tp.x>>Tp.y>>Tp.z;
					xyzVectorE[l]=Tp;
					l++;
				}
				else
				{
					Tp=xyzVectorE[m+i-RecvPlsIgE[npc+ipls]];
				}
				AddReciver(Tp,Av,recE[2*i],sinfiE[2*i],cosfiE[2*i],rc0);
				AddReciver(Tp,Bv,recE[2*i+1],sinfiE[2*i+1],cosfiE[2*i+1],rc0);
				matrec[2*i]=0;
				matrec[2*i+1]=0;
			}
		}
		npc+=nGenByPls[ipm];
	}
	inf.close();
	inf.clear();

	npntE0or=0;
	if(npntE0)
	{
		inf.open("xyzVectorE0");
		if(!inf)
		{
			logfile<<"No file "<<"xyzVectorE0"<<'\n';
			return 1;
		}
		inf>>npntE0or;
		xyzVectorE0=new PointXYZ[npntE0or];
		for(i=0;i<npntE0or;i++)
		{
			inf>>Tp.x>>Tp.y>>Tp.z;
			xyzVectorE0[i]=Tp;
		}
		inf.close();
		inf.clear();
	}

	recE0=new PointRZ[npntE0*2];
	sinfiE0=new double[npntE0*2];
	cosfiE0=new double[npntE0*2];
	matrec0=new int[npntE0*2];
	npc=0;
	for(ipm=0;ipm<npls_main;ipm++)
	{
		for(ipls=0;ipls<nGenByPls[ipm];ipls++)
		{
			l=0;
			for(i=RecvPlsIgE0[npc+ipls];i<RecvPlsIgE0[npc+ipls+1];i++)
			{
				Av=GenSq[RecToSourceE0[i]].A;
				Bv=GenSq[RecToSourceE0[i]].B;
				Tp=xyzVectorE0[l];
				l++;
				AddReciver(Tp,Av,recE0[2*i],sinfiE0[2*i],cosfiE0[2*i],rc0);
				AddReciver(Tp,Bv,recE0[2*i+1],sinfiE0[2*i+1],cosfiE0[2*i+1],rc0);
				matrec0[2*i]=0;
				matrec0[2*i+1]=0;
			}
		}
		npc+=nGenByPls[ipm];
	}



	nt=ntimesdec;

	Hphi=new double *[3];
	v2=new double[kpnt*npls];
	DsFlgNode=new int[kpnt];
	B=new double*[nt];
	F=new double *[nt];
	Er=new double*[nt];
	Ez=new double*[nt];
	Er0=new double[2*npntE0];
	Ez0=new double[2*npntE0];

	for(it=0;it<nt;it++)
	{
		B[it]=new double[2*npntB];
		F[it]=new double[2*npntB];
		Er[it]=new double[2*npntE];
		Ez[it]=new double[2*npntE];
	}

	retp=SobjB.Init(npntB*2,kpnt,krect,recB,pnt,rect,0,NULL,sigma,
		sigmaZ,mu,nreg,qr,qz,reg,rm,zm,1);
	if(retp)
	{
		logfile<<"SobjB.Init() returned "<<retp<<endl;
		return retp;
	}

	retp=SobjE.Init(npntE*2,kpnt,krect,recE,pnt,rect,0,NULL,sigma,
		sigmaZ,mu,nreg,qr,qz,reg,rm,zm,1);
	if(retp)
	{
		logfile<<"SobjE.Init() returned "<<retp<<endl;
		return retp;
	}


	double rmin,zmin,rmax,zmax;

	vRecOutData = new RectOutData[npntE0*2];
	if(!vRecOutData)return RETCODE_NOMEM;

	rmin = zmin = 1e+30;
	rmax = zmax = -1e+30;

	for(i=0;i<npntE0*2;i++)
	{
		if(recE0[i].r<rmin){rmin=recE0[i].r; }
		if(recE0[i].r>rmax){rmax=recE0[i].r; }
		if(recE0[i].z<zmin){zmin=recE0[i].z; }
		if(recE0[i].z>zmax){zmax=recE0[i].z; }
	}

	retp=FindElementsForReceivers(npntE0*2,qr,qz,rm,zm,reg,recE0,vRecOutData,pnt,rect,sigma,sigmaZ,rmin,zmin,rmax,zmax);
	if (retp)
	{
		logfile << "FindElementsForReceivers returned " << retp << endl;
		return retp;
	}


	curr=new double[ntimes];

	inf.open("currentfunction");
	if(!inf)
	{
		logfile<<"No file "<<"currentfunction"<<'\n';
		return 1;
	}
	for(it=0;it<ntimes;it++)
	{
		inf>>tmpd>>curr[it];
	}
	inf.close();
	inf.clear();


	THphi=new double[4*krect*npls];
	if(!THphi)return RETCODE_NOMEM;
	for(it=0;it<3;it++)
	{
		Hphi[2-it]=new double[4*krect*npls];
		if(!Hphi[2-it])return RETCODE_NOMEM;
		GetField(it,v2,Hphi[2-it],str,DsFlgNode,THphi,nGels,nVels,GenToPls,npls_main);
	}

	return 0;
}

void task2d_Hfi::output(int it,vector<float> &enor)
{
	if(tbeg+it>2)
	{
		RollFields(Hphi[0],Hphi[1],Hphi[2]);
		GetField(it,v2,Hphi[0],str,DsFlgNode,THphi,nGels,nVels,GenToPls,npls_main);
	}

	double *Hphi0=(tbeg+it>=2)? Hphi[0] : (tbeg+it==1)? Hphi[1] : Hphi[2];

	SobjB.Output_Bphi(Hphi0,B[it],0,RecToSourceB);

	 if(tbeg+it>=1)
	{
		_diff_t2(F[it],B[it],B[it-1],2*npntB,times,it);
	}
	else
	{
		for(i=0;i<2*npntB;i++){F[it][i]=0.0;}
	}

	if(tbeg+it>=0)
	{
		SobjE.Output_Er(Hphi0,Er[it],1,RecToSourceE);
		SobjE.Output_Ez(Hphi0,Ez[it],1,RecToSourceE);
	}

	for(i=0;i<2*npntE0;i++){Er0[i]=0.0;}
	for(i=0;i<2*npntE0;i++){Ez0[i]=0.0;}

	if(tbeg+it>=1)
	{
		Output_ErForHphi_Pr(npntE0*2,Hphi0,Er0,vRecOutData,RecToSourceE0,krect);
		Output_EzForHphi_Pr(npntE0*2,Hphi0,Ez0,vRecOutData,RecToSourceE0,krect,recE0);
	}

	if(npntE0or)
	{
		kk=0;
		for(k=0;k<2;k++)	//    
		{
			for(ipls=0;ipls<npls;ipls++)
			{
				if(GenTypeTrue[ipls]==k)
				{
					j=0;
					for(i=RecvPlsIgE0[ipls];i<RecvPlsIgE0[ipls+1];i++)
					{
						if(recE0[i].r>DminSrc || (fabs(recE0[i].z-GenSq[ipls].A.z)>DminSrc && fabs(recE0[i].z-GenSq[ipls].B.z)>DminSrc))
						{
							x=(Er0[2*i+1]*cosfiE0[2*i+1]*(1-GenTypeTrue[ipls])-Er0[2*i]*cosfiE0[2*i]);
							y=(Er0[2*i+1]*sinfiE0[2*i+1]*(1-GenTypeTrue[ipls])-Er0[2*i]*sinfiE0[2*i]);
							z=Ez0[2*i+1]*(1-GenTypeTrue[ipls])-Ez0[2*i];
							enor[kk]=float(x*TgCompE0[j].x+y*TgCompE0[j].y+z*TgCompE0[j].z);
						}
						else
						{
							enor[kk]=0.0;
						}
						j++;
						kk++;
					}
				}
			}
		}
	}
}

void task2d_Hfi::output_rec()
{
	sprintf_s(buf, "ball.hphi.stat");
	ofpstat.open(buf);
	sprintf_s(buf, "ball.hphi");
	ofp.open(buf);
	ofp<<scientific<<setprecision(14);
	for(k=0;k<2;k++)	//    
	{
		for(ipls=0;ipls<npls;ipls++)
		{
			if(GenTypeTrue[ipls]==k)
			{
				npc=0;
				for(j=GenToPls[ipls]-1;j>=0;j--)
				{
					npc+=NrecB[j];
				}

				for(i=RecvPlsIgB[ipls];i<RecvPlsIgB[ipls+1];i++)
				{
					Tp=xyzVectorB[npc+i-RecvPlsIgB[ipls]];

					ofp<<(int)Tp.x<<" "<<(int)Tp.y<<endl;
					ofp<<Tp.x<<" "<<Tp.y<<" "<<Tp.z<<endl;
					ofp<<"1"<<endl;
					ofp<<"1  1  1  1"<<endl;
					ofp<<"1"<<endl;
					ofp<<"      t (s)         Bx (T)     By (T)     Bz (T)"<<endl;
					for(it=0;it<ntimesdec;it++)
					{
						const double& tj=times[it+tbeg];
						if (tj>=ftime&&tj<=ltime)
						{
							x=(-B[it][2*i+1]*sinfiB[2*i+1]*(1-GenTypeTrue[ipls])+B[it][2*i]*sinfiB[2*i]);
							y=(B[it][2*i+1]*cosfiB[2*i+1]*(1-GenTypeTrue[ipls])-B[it][2*i]*cosfiB[2*i]);
							ofp<<tj*1000<<"   "<<x*1000<<" "<<y*1000<<" "<<0.0<<'\n';
						}
						if(!it /* && !iDec*/)
						{
							x=(-B[it][2*i+1]*sinfiB[2*i+1]*(1-GenTypeTrue[ipls])+B[it][2*i]*sinfiB[2*i]);
							y=(B[it][2*i+1]*cosfiB[2*i+1]*(1-GenTypeTrue[ipls])-B[it][2*i]*cosfiB[2*i]);
							ofpstat<<x*1000<<" "<<y*1000<<" "<<0.0<<'\n';
						}
					}
				}
			}
		}
	}
	ofp.close();
	ofp.clear();
	ofpstat.close();
	ofpstat.clear();

	sprintf_s(buf, "edsall.hphi");
	ofp.open(buf);
	ofp<<scientific<<setprecision(14);
	for(k=0;k<2;k++)	//    
	{
		for(ipls=0;ipls<npls;ipls++)
		{
			if(GenTypeTrue[ipls]==k)
			{
				npc=0;
				for(j=GenToPls[ipls]-1;j>=0;j--)
				{
					npc+=NrecB[j];
				}

				for(i=RecvPlsIgB[ipls];i<RecvPlsIgB[ipls+1];i++)
				{
					Tp=xyzVectorB[npc+i-RecvPlsIgB[ipls]];

					ofp<<(int)Tp.x<<" "<<(int)Tp.y<<endl;
					ofp<<Tp.x<<" "<<Tp.y<<" "<<Tp.z<<endl;
					ofp<<"1"<<endl;
					ofp<<"1  1  1  1"<<endl;
					ofp<<"1"<<endl;
					ofp<<"      t (s)         Emfx (T)     Emfy (T)     Emfz (T)"<<endl;
					for(it=0;it<ntimesdec;it++)
					{
						const double& tj=times[it+tbeg];
						if (tj>=ftime&&tj<=ltime)
						{
							x=(-F[it][2*i+1]*sinfiB[2*i+1]*(1-GenTypeTrue[ipls])+F[it][2*i]*sinfiB[2*i]);
							y=(F[it][2*i+1]*cosfiB[2*i+1]*(1-GenTypeTrue[ipls])-F[it][2*i]*cosfiB[2*i]);
							ofp<<tj*1000<<"   "<<x*1000<<" "<<y*1000<<" "<<0.0<<'\n';
						}
					}
				}
			}
		}
	}
	ofp.close();
	ofp.clear();


	sprintf_s(buf, "eall.hphi.stat");
	ofpstat.open(buf);
	sprintf_s(buf, "eall.hphi");
	ofp.open(buf);
	ofp<<scientific<<setprecision(14);
	for(k=0;k<2;k++)	//    
	{
		for(ipls=0;ipls<npls;ipls++)
		{
			if(GenTypeTrue[ipls]==k)
			{
				npc=0;
				for(j=GenToPls[ipls]-1;j>=0;j--)
				{
					npc+=NrecE[j];
				}

				for(i=RecvPlsIgE[ipls];i<RecvPlsIgE[ipls+1];i++)
				{
					Tp=xyzVectorE[npc+i-RecvPlsIgE[ipls]];

					ofp<<(int)Tp.x<<" "<<(int)Tp.y<<endl;
					ofp<<Tp.x<<" "<<Tp.y<<" "<<Tp.z<<endl;
					ofp<<"1"<<endl;
					ofp<<"1  1  1  1"<<endl;
					ofp<<"1"<<endl;
					ofp<<"      t (s)         Ex (V/m)     Ey (V/m)     Ez (V/m)"<<endl;
					for(it=0;it<ntimesdec;it++)
					{
						const double& tj=times[it+tbeg];
						if (tj>=ftime&&tj<=ltime)
						{
							x=(Er[it][2*i+1]*cosfiE[2*i+1]*(1-GenTypeTrue[ipls])-Er[it][2*i]*cosfiE[2*i]);
							y=(Er[it][2*i+1]*sinfiE[2*i+1]*(1-GenTypeTrue[ipls])-Er[it][2*i]*sinfiE[2*i]);
							z=Ez[it][2*i+1]*(1-GenTypeTrue[ipls])-Ez[it][2*i];
							ofp<<tj*1000<<"   "<<x*1000<<" "<<y*1000<<" "<<z*1000<<'\n';
						}
						if(!it /* && !iDec*/)
						{
							x=(Er[it][2*i+1]*cosfiE[2*i+1]*(1-GenTypeTrue[ipls])-Er[it][2*i]*cosfiE[2*i]);
							y=(Er[it][2*i+1]*sinfiE[2*i+1]*(1-GenTypeTrue[ipls])-Er[it][2*i]*sinfiE[2*i]);
							z=Ez[it][2*i+1]*(1-GenTypeTrue[ipls])-Ez[it][2*i];
							ofpstat<<x*1000<<" "<<y*1000<<" "<<z*1000<<'\n';
						}
					}
				}
			}
		}
	}
	ofp.close();
	ofp.clear();
	ofpstat.close();
	ofpstat.clear();
}

void task2d_Hfi::clear()
{
	if(B)
	{
		for(it=0;it<nt;it++)
		{
			if(B[it]){delete [] B[it]; B[it]=NULL;}
		}
		delete [] B; B=NULL;
	}
	if(Er)
	{
		for(it=0;it<nt;it++)
		{
			if(Er[it]){delete [] Er[it]; Er[it]=NULL;}
		}
		delete [] Er; Er=NULL;
	}
	if(Ez)
	{
		for(it=0;it<nt;it++)
		{
			if(Ez[it]){delete [] Ez[it]; Ez[it]=NULL;}
		}
		delete [] Ez; Ez=NULL;
	}

	if(recE){delete [] recE; recE=NULL;}
	if(recB){delete [] recB; recB=NULL;}
	if(sinfiE){delete [] sinfiE; sinfiE=NULL;}
	if(cosfiE){delete [] cosfiE; cosfiE=NULL;}
	if(sinfiB){delete [] sinfiB; sinfiB=NULL;}
	if(cosfiB){delete [] cosfiB; cosfiB=NULL;}


	if(F)
	{
		for(it=0;it<nt;it++)
		{
			if(F[it]){delete [] F[it]; F[it]=NULL;}
		}
		delete [] F; F=NULL;
	}

	if(THphi){delete [] THphi; THphi=NULL;}

	if(v2){delete [] v2; v2=NULL;}
	if(DsFlgNode){delete [] DsFlgNode; DsFlgNode=NULL;}
}
