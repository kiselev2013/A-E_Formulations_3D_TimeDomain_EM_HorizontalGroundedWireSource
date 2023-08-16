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

task2d::task2d(char *_path)
{
	strcpy(path,_path);

	rect=NULL;
	pnt=NULL;

	nmat3d=0;
	mtr3d2d=NULL;

	nmat2d=0;
	sigma=NULL;

	reg=NULL;
	rm=NULL;
	zm=NULL;

	v2s=NULL;
	v2c=NULL;

	indr = NULL;
	indz = NULL;
}

task2d::~task2d()
{
	if(rect){delete [] rect; rect=NULL;}
	if(pnt){delete [] pnt; pnt=NULL;}
	if(v2s){delete [] v2s; v2s=NULL;}
	if(v2c){delete [] v2c; v2c=NULL;}
	if(mtr3d2d){delete [] mtr3d2d; mtr3d2d=NULL;}
	if(sigma){delete [] sigma; sigma=NULL;}
	if(reg){delete [] reg; reg=NULL;}
	if(rm){delete [] rm; rm=NULL;}
	if(zm){delete [] zm; zm=NULL;}
	if (indr){ delete[] indr; indr = NULL; }
	if (indz){ delete[] indz; indz = NULL; }
}

int task2d::ReadRZind()
{
	int i, j, k, retc;
	char buf[256];
	ifstream inf;
	double sum;

	if (!(indr = new int[qr])) return RETCODE_NOMEM;
	if (!(indz = new int[qz])) return RETCODE_NOMEM;

	sprintf(buf, "%s\\r_ind.dat", path);
	inf.open(buf, ios::binary);
	if (!inf) return RETCODE_NOFILE;
	for (i = 0; i<qr; i++)
	{
		inf > indr[i];
		indr[i]--;
	}
	inf.close();
	inf.clear();

	sprintf(buf, "%s\\z_ind.dat", path);
	inf.open(buf, ios::binary);
	if (!inf) return RETCODE_NOFILE;
	for (i = 0; i<qz; i++)
	{
		inf > indz[i];
		indz[i]--;
	}
	inf.close();
	inf.clear();

	return 0;
}

int task2d::Read()
{
	int i,j,k,l,retc;
	char buf[256];
	ifstream inf, rzf, nvtrf, nvkatf, l1f, nvk1f;
	FILE *fp;
	double sums,sumc;

	sprintf(buf,"%s\\inf2tr.dat",path);
	inf.open(buf);
	if (!inf) return RETCODE_NOFILE;
	inf.ignore(1000, '\n');
	inf.ignore(1000, '='); inf>>kpnt;
	inf.ignore(1000, '='); inf>>krect;
	inf.close();
	inf.clear();

	rmin=1e+30;
	rmax=-1e+30;
	zmin=1e+30;
	zmax=-1e+30;

	sprintf(buf,"%s\\rz.dat",path);
	rzf.open(buf, ios::binary);
	if (!rzf) return RETCODE_NOFILE;
	pnt=new PointRZ[kpnt];
	if (!pnt) return RETCODE_NOMEM;
	for (i=0; i<kpnt; i++)
	{
		rzf > pnt[i].r;
		rzf > pnt[i].z;

		if(pnt[i].r<rmin){rmin=pnt[i].r;}
		if(pnt[i].r>rmax){rmax=pnt[i].r;}
		if(pnt[i].z<zmin){zmin=pnt[i].z;}
		if(pnt[i].z>zmax){zmax=pnt[i].z;}
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

	sprintf(buf,"%s\\rz.txt",path);
	inf.open(buf);
	if(!inf) return RETCODE_NOFILE;
	inf>>qr>>qz;
	inf.close();
	inf.clear();

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
	for(i=0;i<qz;i++){inf>zm[i];}
	inf.close();
	inf.clear();


	return RETCODE_OK;
}

int task2d::ReadSolution(int npls)
{
	char buf[256];
	FILE *fp;

	v2s = new double[kpnt*npls];
	if (!v2s) return RETCODE_NOMEM;
	v2c = new double[kpnt*npls];
	if (!v2c) return RETCODE_NOMEM;

	sprintf(buf, "%s\\v2s.dat", path);
	fp = fopen(buf, "rb");
	if (!fp) return RETCODE_NOFILE;
	fread(v2s, sizeof(double), kpnt*npls, fp);
	fclose(fp);

	sprintf(buf, "%s\\v2c.dat", path);
	fp = fopen(buf, "rb");
	if (!fp) return RETCODE_NOFILE;
	fread(v2c, sizeof(double), kpnt*npls, fp);
	fclose(fp);
	
	return RETCODE_OK;
}

void task2d::GenKU(int nGels,int nVels,vector<double> &zgel,vector<double> &zvela,vector<double> &zvelb)
{
	char buf[256];
	double zbest;

	int i,j,tmpi,ipls,npls;
	ofstream ofp;
	vector<int> l1flg;
	vector<int> l1;
	vector<int> kt1dc;
	vector<vector<int>> l1dc;

	npls=nGels+nVels;



	kt1dc.resize(npls);
	l1dc.resize(npls);
	for(ipls=0;ipls<npls;ipls++){kt1dc[ipls]=0;}

	for(ipls=0;ipls<nGels;ipls++)
	{
		kt1dc[ipls]=qr;
	}

	for(ipls=0;ipls<nVels;ipls++)
	{
		for(j=0;j<qz;j++)
		{
			kt1dc[nGels+ipls]+=((zm[j]+1e-4>zvela[ipls]) && (zm[j]-1e-4<zvelb[ipls]));
		}
	}
	for(ipls=0;ipls<npls;ipls++){l1dc[ipls].resize(kt1dc[ipls]);}
	for(ipls=0;ipls<npls;ipls++){kt1dc[ipls]=0;}

	for(ipls=0;ipls<nGels;ipls++)
	{
		zbest=1e+30;
		for(j=0;j<qz;j++)
		{
			if(fabs(zm[j]-zgel[ipls])<fabs(zbest-zgel[ipls]))
			{
				zbest=zm[j];
			}
		}

		for(j=0;j<qz;j++)
		{
			if(zbest==zm[j])
			{
				break;
			}
		}

		for(i=0;i<qr;i++)
		{
			l1dc[ipls][kt1dc[ipls]]=j*qr+i;
			kt1dc[ipls]++;
		}
	}

	for(ipls=0;ipls<nVels;ipls++)
	{
		i=0;
		for(j=0;j<qz;j++)
		{
			if((zm[j]+1e-4>zvela[ipls]) && (zm[j]-1e-4<zvelb[ipls]))
			{
				l1dc[nGels+ipls][kt1dc[nGels+ipls]]=j*qr+i;
				kt1dc[nGels+ipls]++;
			}
		}
	}

	sprintf(buf, "%s\\kt1dc", path);
	ofp.open(buf);
	for(ipls=0;ipls<npls;ipls++){ofp<<kt1dc[ipls]<<'\n';}
	ofp.close();
	ofp.clear();

	sprintf(buf, "%s\\l1dc.dat", path);
	ofp.open(buf,ios::binary);
	for(ipls=0;ipls<npls;ipls++)
	{
		for(i=0;i<kt1dc[ipls];i++)
		{
			tmpi=l1dc[ipls][i]+1;
			ofp.write((char *)&(tmpi),size_i);
		}
	}
	ofp.close();
	ofp.clear();
}
