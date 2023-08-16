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
 *  This file contains code for smooth outputing desired primary field components functions
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 January 16, 2023
*/

#include "stdafx.h"
#include "SmoothOutput2D.h"
#include <omp.h> 

extern ofstream logfile;

SmoothOutput2D::SmoothOutput2D()
{
	nrec=kpnt=krect=0;
	rec=NULL;
	pnt=NULL;
	rect=NULL;
	OutRecList=NULL;
	RecToElem.clear();
	int nthreads;
	nthreads = omp_get_max_threads();
	nthreads = (nthreads>3)? 3 : nthreads;
	omp_set_num_threads(nthreads);
}

SmoothOutput2D::~SmoothOutput2D()
{
	RecToElem.clear();
	if(OutRecList){delete [] OutRecList; OutRecList=NULL;}
}

int SmoothOutput2D::Init(int p_nrec,int p_kpnt,int p_krect,PointRZ *p_rec,PointRZ *p_pnt,Rect *p_rect,int MatFlag,int *matrec,double *p_sigma,
						 double *p_sigmaZ,double *p_mu,int p_nreg,int p_qr,int p_qz,int *p_reg,double *p_rm,double *p_zm,int stype)
{
	int i,j,k,l,m,ii,jj,ef;
	PointRZ celm;
	
	bool fstop;

	vector<int> TmpIntVec;
	vector<int> ElemForPoints;
	vector<int> ElemForPointsIg;

	nrec=p_nrec;
	kpnt=p_kpnt;
	krect=p_krect;
	rec=p_rec;
	pnt=p_pnt;
	rect=p_rect;
	sigma=p_sigma;
	sigmaZ=p_sigmaZ;
	mu=p_mu;
	nreg=p_nreg;
	qr=p_qr;
	qz=p_qz;
	reg=p_reg;
	rm=p_rm;
	zm=p_zm;

	if(nrec==0)return 0;

	if(!((int)RecToElem.size()))RecToElem.resize(nrec);
	if(!OutRecList)OutRecList=new Output2D[nrec];

	fstop=false;
	for(k=0;k<nrec;k++)
	{
		i=FindIntervalInDoubleMas(rm,qr,rec[k].r);
		j=FindIntervalInDoubleMas(zm,qz,rec[k].z);
		if(i!=-1 && j!=-1)
		{
			m=j*(qr-1)+i;
			l=reg[m]-1;

			if(matrec)
			{
				if(matrec[k] && rect[l].mtr!=matrec[k])
				{
					jj=0;
					if(fabs(rec[k].z-pnt[rect[l].nodes[0]-1].z)<1e-6)
					{
						if(j>0)
						{
							m=(j-1)*(qr-1)+i;
							l=reg[m]-1;
							if(rect[l].mtr!=matrec[k])
							{
								jj=1;
							}
							else
							{
								j--;
							}
						}
						else
						{
							jj=1;
						}					
					}
					else if(fabs(rec[k].z-pnt[rect[l].nodes[3]-1].z)<1e-6)
					{
						if(j<qz-1)
						{
							m=(j+1)*(qr-1)+i;
							l=reg[m]-1;
							if(rect[l].mtr!=matrec[k])
							{
								jj=1;
							}
							else
							{
								j++;
							}
						}
						else
						{
							jj=1;
						}
					}
					else
					{
						jj=1;
					}
				
					if(jj)
					{
						logfile<<"Material problem for reciver "<<k+1<<" r= "<<rec[k].r<<" z= "<<rec[k].z<<'\n';
						fstop=true;
						exit(1);
					}
				}
			}

			RecToElem[k]=l;

			celm.r=0.0;
			celm.z=0.0;
			for(ii=0;ii<4;ii++)
			{
				OutRecList[k].or.lem[ii]=-1;
				OutRecList[k].oz.lem[ii]=-1;
				jj=rect[l].nodes[ii]-1;
				celm.r+=pnt[jj].r;
				celm.z+=pnt[jj].z;
			}
			celm.r*=0.25;
			celm.z*=0.25;

			ef=!(rec[k].r>celm.r);
			GetSmoothForR(i,j,k,ef,l,0,rec[k],OutRecList[k].or.lem,OutRecList[k].or.coef,OutRecList[k].or.v,OutRecList[k].or.ves,OutRecList[k].or.fl,stype);

			ef=!(rec[k].z>celm.z);
			GetSmoothForZ(i,j,k,ef,l,0,rec[k],OutRecList[k].oz.lem,OutRecList[k].oz.coef,OutRecList[k].oz.v,OutRecList[k].oz.ves,OutRecList[k].oz.fl,stype);
		}
		else
		{
			RecToElem[k]=-1;
		}
	}

	return fstop;
}

int SmoothOutput2D::Init(int p_nrec,int p_kpnt,int p_krect,PointRZ *p_rec,PointRZ *p_pnt,Rect *p_rect,int MatFlag,int *matrec,double *p_sigma,
						 double *p_sigmaZ,double *p_mu,int p_nreg,int p_qr,int p_qz,int *p_reg,double *p_rm,double *p_zm,int stype,vector<int> &OutCompFlag)
{
	int i,j,k,l,m,ii,jj,ef;
	PointRZ celm;
	
	bool fstop;

	vector<int> TmpIntVec;			//  
	vector<int> ElemForPoints;		//    
	vector<int> ElemForPointsIg;	//        ig

	nrec=p_nrec;
	kpnt=p_kpnt;
	krect=p_krect;
	rec=p_rec;
	pnt=p_pnt;
	rect=p_rect;
	sigma=p_sigma;
	sigmaZ=p_sigmaZ;
	mu=p_mu;
	nreg=p_nreg;
	qr=p_qr;
	qz=p_qz;
	reg=p_reg;
	rm=p_rm;
	zm=p_zm;

	if(nrec==0)return 0;

	if(!((int)RecToElem.size()))RecToElem.resize(nrec);
	if(!OutRecList)OutRecList=new Output2D[nrec];

	fstop=false;
	for(k=0;k<nrec;k++)
	{
		if(OutCompFlag[k/2])
		{
			i=FindIntervalInDoubleMas(rm,qr,rec[k].r);
			j=FindIntervalInDoubleMas(zm,qz,rec[k].z);
			if(i!=-1 && j!=-1)
			{
				m=j*(qr-1)+i;
				l=reg[m]-1;

				if(matrec)
				{
					if(matrec[k] && rect[l].mtr!=matrec[k])
					{
						jj=0;
						if(fabs(rec[k].z-pnt[rect[l].nodes[0]-1].z)<1e-6)
						{
							if(j>0)
							{
								m=(j-1)*(qr-1)+i;
								l=reg[m]-1;
								if(rect[l].mtr!=matrec[k])
								{
									jj=1;
								}
								else
								{
									j--;
								}
							}
							else
							{
								jj=1;
							}					
						}
						else if(fabs(rec[k].z-pnt[rect[l].nodes[3]-1].z)<1e-6)
						{
							if(j<qz-1)
							{
								m=(j+1)*(qr-1)+i;
								l=reg[m]-1;
								if(rect[l].mtr!=matrec[k])
								{
									jj=1;
								}
								else
								{
									j++;
								}
							}
							else
							{
								jj=1;
							}
						}
						else
						{
							jj=1;
						}
					
						if(jj)
						{
							logfile<<"Material problem for reciver "<<k+1<<" r= "<<rec[k].r<<" z= "<<rec[k].z<<'\n';
							fstop=true;
							exit(1);
						}
					}
				}

				RecToElem[k]=l;

				celm.r=0.0;
				celm.z=0.0;
				for(ii=0;ii<4;ii++)
				{
					OutRecList[k].or.lem[ii]=-1;
					OutRecList[k].oz.lem[ii]=-1;
					jj=rect[l].nodes[ii]-1;
					celm.r+=pnt[jj].r;
					celm.z+=pnt[jj].z;
				}
				celm.r*=0.25;
				celm.z*=0.25;

				ef=!(rec[k].r>celm.r);
				GetSmoothForR(i,j,k,ef,l,0,rec[k],OutRecList[k].or.lem,OutRecList[k].or.coef,OutRecList[k].or.v,OutRecList[k].or.ves,OutRecList[k].or.fl,stype);

				ef=!(rec[k].z>celm.z);
				GetSmoothForZ(i,j,k,ef,l,0,rec[k],OutRecList[k].oz.lem,OutRecList[k].oz.coef,OutRecList[k].oz.v,OutRecList[k].oz.ves,OutRecList[k].oz.fl,stype);
			}
			else
			{
				RecToElem[k]=-1;
			}
		}
	}

	return fstop;
}

void SmoothOutput2D::GetSmoothForR(int i,int j,int k,int ef,int l,int MatFlag,PointRZ &reck,int rlem[4],float rcoef[4][4],float rv[2][3],float &ves,int &fl,int stype)
{
	int n,ii,kk,ll,lt,mm,tt;
	double cr[3],r0,r3,z0,z3,invhrhz;
	bool fpp[2];
	int lef;
	PointRZ clem,lrecz;

	fl=0;
	ves=0.0;

	for(kk=0;kk<4;kk++)rlem[kk]=-1;

	rlem[1+ef]=l;

	if(stype)
	{
		clem.r=0.0;
		clem.z=0.0;
		for(kk=0;kk<4;kk++)
		{
			clem.r+=pnt[rect[l].nodes[kk]-1].r;
			clem.z+=pnt[rect[l].nodes[kk]-1].z;
		}
		clem.r*=0.25;
		clem.z*=0.25;

		lrecz.r=clem.r;
		lrecz.z=reck.z;
		lef=!(lrecz.z>clem.z);

		tt=1+ef;
		GetSmoothForZ(i,j,k,lef,l,0,lrecz,OutRecList[k].ozfr[tt].lem,OutRecList[k].ozfr[tt].coef,OutRecList[k].ozfr[tt].v,OutRecList[k].ozfr[tt].ves,OutRecList[k].ozfr[tt].fl,0);
	}

	lt=l;
	n=0;
	for(ii=i-1;ii>=0 && n<(1+ef);ii--)
	{
		mm=j*(qr-1)+ii;
		ll=reg[mm]-1;
		if(ll!=lt)
		{
			if(MatFlag && rect[l].mtr!=rect[ll].mtr)break;
			n+=1;
			rlem[(1+ef)-n]=ll;
			lt=ll;

			if(stype)
			{
				clem.r=0.0;
				clem.z=0.0;
				for(kk=0;kk<4;kk++)
				{
					clem.r+=pnt[rect[ll].nodes[kk]-1].r;
					clem.z+=pnt[rect[ll].nodes[kk]-1].z;
				}
				clem.r*=0.25;
				clem.z*=0.25;

				lrecz.r=clem.r;
				lrecz.z=reck.z;
				lef=!(lrecz.z>clem.z);

				tt=(1+ef)-n;
				GetSmoothForZ(ii,j,k,lef,ll,0,lrecz,OutRecList[k].ozfr[tt].lem,OutRecList[k].ozfr[tt].coef,OutRecList[k].ozfr[tt].v,OutRecList[k].ozfr[tt].ves,OutRecList[k].ozfr[tt].fl,0);
			}
		}
	}

	if(ef && rlem[1]==-1)
	{
		rlem[1]=rlem[1+ef];
		rlem[1+ef]=-1;

		if(stype)
		{
			OutRecList[k].ozfr[1].fl=OutRecList[k].ozfr[1+ef].fl;
			OutRecList[k].ozfr[1].ves=OutRecList[k].ozfr[1+ef].ves;

			for(ii=0;ii<2;ii++)
			{
				for(kk=0;kk<3;kk++)
				{
					OutRecList[k].ozfr[1].v[ii][kk]=OutRecList[k].ozfr[1+ef].v[ii][kk];
				}
			}

			for(ii=0;ii<4;ii++)
			{
				OutRecList[k].ozfr[1].lem[ii]=OutRecList[k].ozfr[1+ef].lem[ii];
				if(OutRecList[k].ozfr[1].lem[ii]!=-1)
				{
					for(kk=0;kk<4;kk++)
					{
						OutRecList[k].ozfr[1].coef[ii][kk]=OutRecList[k].ozfr[1+ef].coef[ii][kk];
					}
				}
			}
		}
		ef=0;
	}

	lt=l;
	n=0;
	for(ii=i+1;ii<qr-1 && n<(2-ef);ii++)
	{
		mm=j*(qr-1)+ii;
		ll=reg[mm]-1;
		if(ll!=lt)
		{
			if(MatFlag && rect[l].mtr!=rect[ll].mtr)break;
			n+=1;
			rlem[(1+ef)+n]=ll;
			lt=ll;

			if(stype)
			{
				clem.r=0.0;
				clem.z=0.0;
				for(kk=0;kk<4;kk++)
				{
					clem.r+=pnt[rect[ll].nodes[kk]-1].r;
					clem.z+=pnt[rect[ll].nodes[kk]-1].z;
				}
				clem.r*=0.25;
				clem.z*=0.25;

				lrecz.r=clem.r;
				lrecz.z=reck.z;
				lef=!(lrecz.z>clem.z);

				tt=(1+ef)+n;
				GetSmoothForZ(ii,j,k,lef,ll,0,lrecz,OutRecList[k].ozfr[tt].lem,OutRecList[k].ozfr[tt].coef,OutRecList[k].ozfr[tt].v,OutRecList[k].ozfr[tt].ves,OutRecList[k].ozfr[tt].fl,0);
			}
		}
	}

	fpp[0]=(rlem[0]!=-1 && rlem[1]!=-1 && rlem[2]!=-1);
	fpp[1]=(rlem[1]!=-1 && rlem[2]!=-1 && rlem[3]!=-1);

	if(fpp[0])
	{
		for(kk=0;kk<3;kk++)
		{
			cr[kk]=0.0;
			ll=rlem[kk];
			for(ii=0;ii<4;ii++)
			{
				cr[kk]+=pnt[rect[ll].nodes[ii]-1].r;
			}
			cr[kk]*=0.25;
			
			r0=pnt[rect[ll].nodes[0]-1].r;
			z0=pnt[rect[ll].nodes[0]-1].z;
			r3=pnt[rect[ll].nodes[3]-1].r;
			z3=pnt[rect[ll].nodes[3]-1].z;

			invhrhz=1.0/((r3-r0)*(z3-z0));

			rcoef[kk][0]=((r3-cr[kk])*(z3-reck.z))*invhrhz;
			rcoef[kk][1]=((cr[kk]-r0)*(z3-reck.z))*invhrhz;
			rcoef[kk][2]=((r3-cr[kk])*(reck.z-z0))*invhrhz;
			rcoef[kk][3]=((cr[kk]-r0)*(reck.z-z0))*invhrhz;
		}
		if(!fpp[1])
		{
			fl=1;
			ves=1.0;
		}
		else
		{
			fl=2;
			ves=1.0-(reck.r-cr[1])/(cr[2]-cr[1]);
		}
		rv[0][0]=((reck.r-cr[1])*(reck.r-cr[2]))/((cr[1]-cr[0])*(cr[2]-cr[0]));
		rv[0][1]=((reck.r-cr[0])*(reck.r-cr[2]))/((cr[1]-cr[0])*(cr[1]-cr[2]));
		rv[0][2]=((reck.r-cr[0])*(reck.r-cr[1]))/((cr[2]-cr[0])*(cr[2]-cr[1]));
	}

	if(fpp[1])
	{
		if(!fpp[0])
		{
			for(kk=0;kk<3;kk++)
			{
				cr[kk]=0.0;
				ll=rlem[kk+1];
				for(ii=0;ii<4;ii++)
				{
					cr[kk]+=pnt[rect[ll].nodes[ii]-1].r;
				}
				cr[kk]*=0.25;
				
				r0=pnt[rect[ll].nodes[0]-1].r;
				z0=pnt[rect[ll].nodes[0]-1].z;
				r3=pnt[rect[ll].nodes[3]-1].r;
				z3=pnt[rect[ll].nodes[3]-1].z;

				invhrhz=1.0/((r3-r0)*(z3-z0));

				rcoef[kk+1][0]=((r3-cr[kk])*(z3-reck.z))*invhrhz;
				rcoef[kk+1][1]=((cr[kk]-r0)*(z3-reck.z))*invhrhz;
				rcoef[kk+1][2]=((r3-cr[kk])*(reck.z-z0))*invhrhz;
				rcoef[kk+1][3]=((cr[kk]-r0)*(reck.z-z0))*invhrhz;
			}
			fl=1;
			ves=0.0;
		}
		else
		{
			cr[0]=cr[1];
			cr[1]=cr[2];
			cr[2]=0.0;
			ll=rlem[3];
			for(ii=0;ii<4;ii++)
			{
				cr[2]+=pnt[rect[ll].nodes[ii]-1].r;
			}
			cr[2]*=0.25;
			
			r0=pnt[rect[ll].nodes[0]-1].r;
			z0=pnt[rect[ll].nodes[0]-1].z;
			r3=pnt[rect[ll].nodes[3]-1].r;
			z3=pnt[rect[ll].nodes[3]-1].z;

			invhrhz=1.0/((r3-r0)*(z3-z0));

			rcoef[3][0]=((r3-cr[2])*(z3-reck.z))*invhrhz;
			rcoef[3][1]=((cr[2]-r0)*(z3-reck.z))*invhrhz;
			rcoef[3][2]=((r3-cr[2])*(reck.z-z0))*invhrhz;
			rcoef[3][3]=((cr[2]-r0)*(reck.z-z0))*invhrhz;
		}
		rv[1][0]=((reck.r-cr[1])*(reck.r-cr[2]))/((cr[1]-cr[0])*(cr[2]-cr[0]));
		rv[1][1]=((reck.r-cr[0])*(reck.r-cr[2]))/((cr[1]-cr[0])*(cr[1]-cr[2]));
		rv[1][2]=((reck.r-cr[0])*(reck.r-cr[1]))/((cr[2]-cr[0])*(cr[2]-cr[1]));
	}

	if(stype)
	{
		for(kk=0;kk<4;kk++)
		{
			if(rlem[kk]==-1)OutRecList[k].ozfr[kk].fl=0;
		}
	}
}

void SmoothOutput2D::GetSmoothForZ(int i,int j,int k,int ef,int l,int MatFlag,PointRZ &reck,int zlem[4],float zcoef[4][4],float zv[2][3],float &ves,int &fl,int stype)
{
	int n,ii,kk,ll,lt,mm,tt;
	double cz[3],r0,r3,z0,z3,invhrhz;
	bool fpp[2];
	int lef;
	PointRZ clem,lrecr;

	fl=0;
	ves=0.0;

	for(kk=0;kk<4;kk++)zlem[kk]=-1;

	zlem[1+ef]=l;

	if(stype)
	{
		clem.r=0.0;
		clem.z=0.0;
		for(kk=0;kk<4;kk++)
		{
			clem.r+=pnt[rect[l].nodes[kk]-1].r;
			clem.z+=pnt[rect[l].nodes[kk]-1].z;
		}
		clem.r*=0.25;
		clem.z*=0.25;

			lrecr.r=reck.r;
			lrecr.z=clem.z;
			lef=!(lrecr.r>clem.r);

		tt=1+ef;
		GetSmoothForR(i,j,k,lef,l,0,lrecr,OutRecList[k].orfz[tt].lem,OutRecList[k].orfz[tt].coef,OutRecList[k].orfz[tt].v,OutRecList[k].orfz[tt].ves,OutRecList[k].orfz[tt].fl,0);
	}

	lt=l;
	n=0;
	for(ii=j-1;ii>=0 && n<(1+ef);ii--)
	{
		mm=ii*(qr-1)+i;
		ll=reg[mm]-1;
		if(ll!=lt)
		{
			if(MatFlag && rect[l].mtr!=rect[ll].mtr)break;
			n+=1;
			zlem[(1+ef)-n]=ll;
			lt=ll;

			if(stype)
			{
				clem.r=0.0;
				clem.z=0.0;
				for(kk=0;kk<4;kk++)
				{
					clem.r+=pnt[rect[ll].nodes[kk]-1].r;
					clem.z+=pnt[rect[ll].nodes[kk]-1].z;
				}
				clem.r*=0.25;
				clem.z*=0.25;

				lrecr.r=reck.r;
				lrecr.z=clem.z;
				lef=!(lrecr.r>clem.r);

				tt=(1+ef)-n;
				GetSmoothForR(i,ii,k,lef,ll,0,lrecr,OutRecList[k].orfz[tt].lem,OutRecList[k].orfz[tt].coef,OutRecList[k].orfz[tt].v,OutRecList[k].orfz[tt].ves,OutRecList[k].orfz[tt].fl,0);
			}
		}
	}

	if(ef && zlem[1]==-1)
	{
		zlem[1]=zlem[1+ef];
		zlem[1+ef]=-1;
		
		if(stype)
		{
			OutRecList[k].orfz[1].fl=OutRecList[k].orfz[1+ef].fl;
			OutRecList[k].orfz[1].ves=OutRecList[k].orfz[1+ef].ves;

			for(ii=0;ii<2;ii++)
			{
				for(kk=0;kk<3;kk++)
				{
					OutRecList[k].orfz[1].v[ii][kk]=OutRecList[k].orfz[1+ef].v[ii][kk];
				}
			}

			for(ii=0;ii<4;ii++)
			{
				OutRecList[k].orfz[1].lem[ii]=OutRecList[k].orfz[1+ef].lem[ii];
				if(OutRecList[k].orfz[1].lem[ii]!=-1)
				{
					for(kk=0;kk<4;kk++)
					{
						OutRecList[k].orfz[1].coef[ii][kk]=OutRecList[k].orfz[1+ef].coef[ii][kk];
					}
				}
			}
		}
		ef=0;
	}

	lt=l;
	n=0;
	for(ii=j+1;ii<qz-1 && n<(2-ef);ii++)
	{
		mm=ii*(qr-1)+i;
		ll=reg[mm]-1;
		if(ll!=lt)
		{
			if(MatFlag && rect[l].mtr!=rect[ll].mtr)break;
			n+=1;
			zlem[(1+ef)+n]=ll;
			lt=ll;

			if(stype)
			{
				clem.r=0.0;
				clem.z=0.0;
				for(kk=0;kk<4;kk++)
				{
					clem.r+=pnt[rect[ll].nodes[kk]-1].r;
					clem.z+=pnt[rect[ll].nodes[kk]-1].z;
				}
				clem.r*=0.25;
				clem.z*=0.25;

				lrecr.r=reck.r;
				lrecr.z=clem.z;
				lef=!(lrecr.r>clem.r);

				tt=(1+ef)+n;
				GetSmoothForR(i,ii,k,lef,ll,0,lrecr,OutRecList[k].orfz[tt].lem,OutRecList[k].orfz[tt].coef,OutRecList[k].orfz[tt].v,OutRecList[k].orfz[tt].ves,OutRecList[k].orfz[tt].fl,0);
			}
		}
	}

	fpp[0]=(zlem[0]!=-1 && zlem[1]!=-1 && zlem[2]!=-1);
	fpp[1]=(zlem[1]!=-1 && zlem[2]!=-1 && zlem[3]!=-1);

	if(fpp[0])
	{
		for(kk=0;kk<3;kk++)
		{
			cz[kk]=0.0;
			ll=zlem[kk];
			for(ii=0;ii<4;ii++)
			{
				cz[kk]+=pnt[rect[ll].nodes[ii]-1].z;
			}
			cz[kk]*=0.25;
			
			r0=pnt[rect[ll].nodes[0]-1].r;
			z0=pnt[rect[ll].nodes[0]-1].z;
			r3=pnt[rect[ll].nodes[3]-1].r;
			z3=pnt[rect[ll].nodes[3]-1].z;

			invhrhz=1.0/((r3-r0)*(z3-z0));

			zcoef[kk][0]=((r3-reck.r)*(z3-cz[kk]))*invhrhz;
			zcoef[kk][1]=((reck.r-r0)*(z3-cz[kk]))*invhrhz;
			zcoef[kk][2]=((r3-reck.r)*(cz[kk]-z0))*invhrhz;
			zcoef[kk][3]=((reck.r-r0)*(cz[kk]-z0))*invhrhz;
		}
		if(!fpp[1])
		{
			fl=1;
			ves=1.0;
		}
		else
		{
			fl=2;
			ves=1.0-(reck.z-cz[1])/(cz[2]-cz[1]);
		}
		zv[0][0]=((reck.z-cz[1])*(reck.z-cz[2]))/((cz[1]-cz[0])*(cz[2]-cz[0]));
		zv[0][1]=((reck.z-cz[0])*(reck.z-cz[2]))/((cz[1]-cz[0])*(cz[1]-cz[2]));
		zv[0][2]=((reck.z-cz[0])*(reck.z-cz[1]))/((cz[2]-cz[0])*(cz[2]-cz[1]));
	}

	if(fpp[1])
	{
		if(!fpp[0])
		{
			for(kk=0;kk<3;kk++)
			{
				cz[kk]=0.0;
				ll=zlem[kk+1];
				for(ii=0;ii<4;ii++)
				{
					cz[kk]+=pnt[rect[ll].nodes[ii]-1].z;
				}
				cz[kk]*=0.25;
				
				r0=pnt[rect[ll].nodes[0]-1].r;
				z0=pnt[rect[ll].nodes[0]-1].z;
				r3=pnt[rect[ll].nodes[3]-1].r;
				z3=pnt[rect[ll].nodes[3]-1].z;

				invhrhz=1.0/((r3-r0)*(z3-z0));

				zcoef[kk+1][0]=((r3-reck.r)*(z3-cz[kk]))*invhrhz;
				zcoef[kk+1][1]=((reck.r-r0)*(z3-cz[kk]))*invhrhz;
				zcoef[kk+1][2]=((r3-reck.r)*(cz[kk]-z0))*invhrhz;
				zcoef[kk+1][3]=((reck.r-r0)*(cz[kk]-z0))*invhrhz;
			}
			fl=1;
			ves=0.0;
		}
		else
		{
			cz[0]=cz[1];
			cz[1]=cz[2];
			cz[2]=0.0;
			ll=zlem[3];
			for(ii=0;ii<4;ii++)
			{
				cz[2]+=pnt[rect[ll].nodes[ii]-1].z;
			}
			cz[2]*=0.25;
			
			r0=pnt[rect[ll].nodes[0]-1].r;
			z0=pnt[rect[ll].nodes[0]-1].z;
			r3=pnt[rect[ll].nodes[3]-1].r;
			z3=pnt[rect[ll].nodes[3]-1].z;

			invhrhz=1.0/((r3-r0)*(z3-z0));

			zcoef[3][0]=((r3-reck.r)*(z3-cz[2]))*invhrhz;
			zcoef[3][1]=((reck.r-r0)*(z3-cz[2]))*invhrhz;
			zcoef[3][2]=((r3-reck.r)*(cz[2]-z0))*invhrhz;
			zcoef[3][3]=((reck.r-r0)*(cz[2]-z0))*invhrhz;
		}
		zv[1][0]=((reck.z-cz[1])*(reck.z-cz[2]))/((cz[1]-cz[0])*(cz[2]-cz[0]));
		zv[1][1]=((reck.z-cz[0])*(reck.z-cz[2]))/((cz[1]-cz[0])*(cz[1]-cz[2]));
		zv[1][2]=((reck.z-cz[0])*(reck.z-cz[1]))/((cz[2]-cz[0])*(cz[2]-cz[1]));
	}

	if(stype)
	{
		for(kk=0;kk<4;kk++)
		{
			if(zlem[kk]==-1)OutRecList[k].orfz[kk].fl=0;
		}
	}
}


void SmoothOutput2D_Ax::OutputE(double *_v2,double *result,int OutType,vector<int> &RecToSourceE)
{
	int i,j,k,l,m,it,nn[2];
	double s,ls,v,hr,hz,psi,eta;
	int ipls;
	double *v2;

	for(k=0;k<nrec;k++)
	{
		ipls=RecToSourceE[k/NUMBEROFABDIPOLES];
		v2=_v2+kpnt*ipls;

		if(OutType && OutRecList[k].or.fl)
		{
			s=0.0;
			for(it=0;it<2;it++)
			{
				if((!it && OutRecList[k].or.ves) || (it && OutRecList[k].or.ves!=1.0))
				{					
					ls=0.0;
					for(j=0;j<3;j++)
					{
						m=OutRecList[k].or.lem[j+it];
						
						v=0.0;
						for(l=0;l<4;l++)
						{
							v+=OutRecList[k].or.coef[j+it][l]*v2[rect[m].nodes[l]-1];
						}	
						
						ls+=OutRecList[k].or.v[it][j]*v;
					}

					s+=ls*((it)? 1.0-OutRecList[k].or.ves : OutRecList[k].or.ves);
				}
			}
		}
		else
		{
			i=RecToElem[k];

			nn[0]=rect[i].nodes[0]-1;
			nn[1]=rect[i].nodes[3]-1;

			hr=pnt[nn[1]].r-pnt[nn[0]].r;
			hz=pnt[nn[1]].z-pnt[nn[0]].z;
			psi=(rec[k].r-pnt[nn[0]].r)/hr;
			eta=(rec[k].z-pnt[nn[0]].z)/hz;

			s=(1.0-psi)*(1.0-eta)*v2[rect[i].nodes[0]-1]+psi*(1.0-eta)*v2[rect[i].nodes[1]-1]
				+(1.0-psi)*eta*v2[rect[i].nodes[2]-1]+psi*eta*v2[rect[i].nodes[3]-1];
		}

		result[k]=s;
	}
}

void SmoothOutput2D_Ax::OutputB(double *_v2,double *resultx,double *resulty,
	double *resultz,vector<PointXYZ> &Ic,double *sinfi,double *cosfi,int OutType,vector<int> &RecToSourceB)
{
	double dAfidr,dAfidz,drx[2],dry[2],pr,drdx,drdy,dc=0.01;
	int i,j,k,m,it,nn[2];
	double s,ls,v,hr,hz,psi,eta;
	int ipls;
	double *v2;
	double Ix,Iy,Iz;

	for(k=0;k<nrec;k++)
	{
		j=k/NUMBEROFABDIPOLES;
		Ix=Ic[j].x;
		Iy=Ic[j].y;
		Iz=Ic[j].z;
		ipls=RecToSourceB[j];
		v2=_v2+kpnt*ipls;

		if(OutType && OutRecList[k].or.fl)
		{
			s=0.0;
			for(it=0;it<2;it++)
			{
				if((!it && OutRecList[k].or.ves) || (it && OutRecList[k].or.ves!=1.0))
				{
					ls=0.0;				
					for(j=0;j<3;j++)
					{
						m=OutRecList[k].or.lem[j+it];

						nn[0]=rect[m].nodes[0]-1;
						nn[1]=rect[m].nodes[3]-1;

						hr=pnt[nn[1]].r-pnt[nn[0]].r;
						hz=pnt[nn[1]].z-pnt[nn[0]].z;
						eta=(rec[k].z-pnt[nn[0]].z)/hz;

						v=(-(1.0-eta)*v2[rect[m].nodes[0]-1]+(1.0-eta)*v2[rect[m].nodes[1]-1]
							-eta*v2[rect[m].nodes[2]-1]+eta*v2[rect[m].nodes[3]-1])/hr;

						ls+=OutRecList[k].or.v[it][j]*v;
					}
					s+=ls*((it)? 1.0-OutRecList[k].or.ves : OutRecList[k].or.ves);
				}
			}
		}
		else
		{
			i=RecToElem[k];

			nn[0]=rect[i].nodes[0]-1;
			nn[1]=rect[i].nodes[3]-1;

			hr=pnt[nn[1]].r-pnt[nn[0]].r;
			hz=pnt[nn[1]].z-pnt[nn[0]].z;
			eta=(rec[k].z-pnt[nn[0]].z)/hz;

			s=(-(1.0-eta)*v2[rect[i].nodes[0]-1]+(1.0-eta)*v2[rect[i].nodes[1]-1]
				-eta*v2[rect[i].nodes[2]-1]+eta*v2[rect[i].nodes[3]-1])/hr;
		}

		dAfidr=s;


		if(OutType && OutRecList[k].oz.fl)
		{
			s=0.0;
			for(it=0;it<2;it++)
			{
				if((!it && OutRecList[k].oz.ves) || (it && OutRecList[k].oz.ves!=1.0))
				{
					ls=0.0;				
					for(j=0;j<3;j++)
					{
						m=OutRecList[k].oz.lem[j+it];

						nn[0]=rect[m].nodes[0]-1;
						nn[1]=rect[m].nodes[3]-1;

						hr=pnt[nn[1]].r-pnt[nn[0]].r;
						hz=pnt[nn[1]].z-pnt[nn[0]].z;
						psi=(rec[k].r-pnt[nn[0]].r)/hr;

						v=(-(1.0-psi)*v2[rect[m].nodes[0]-1]-psi*v2[rect[m].nodes[1]-1]
								+(1.0-psi)*v2[rect[m].nodes[2]-1]+psi*v2[rect[m].nodes[3]-1])/hz;

						ls+=OutRecList[k].oz.v[it][j]*v;
					}
					s+=ls*((it)? 1.0-OutRecList[k].oz.ves : OutRecList[k].oz.ves);
				}
			}
		}
		else
		{
			i=RecToElem[k];

			nn[0]=rect[i].nodes[0]-1;
			nn[1]=rect[i].nodes[3]-1;

			hr=pnt[nn[1]].r-pnt[nn[0]].r;
			hz=pnt[nn[1]].z-pnt[nn[0]].z;
			psi=(rec[k].r-pnt[nn[0]].r)/hr;

			s=(-(1.0-psi)*v2[rect[i].nodes[0]-1]-psi*v2[rect[i].nodes[1]-1]
				+(1.0-psi)*v2[rect[i].nodes[2]-1]+psi*v2[rect[i].nodes[3]-1])/hz;
		}

		dAfidz=s;

		pr=rec[k].r;
		if(pr<10.0*dc)pr=10.0*dc;

		drx[0]=sqrt(pr*pr+dc*dc+2*pr*dc*cosfi[k])-pr;
		drx[1]=pr-sqrt(pr*pr+dc*dc-2*pr*dc*cosfi[k]);
		dry[0]=sqrt(pr*pr+dc*dc+2*pr*dc*sinfi[k])-pr;
		dry[1]=pr-sqrt(pr*pr+dc*dc-2*pr*dc*sinfi[k]);

		drdx=(drx[0]+drx[1])/(2.0*dc);
		drdy=(dry[0]+dry[1])/(2.0*dc);

		resultx[k]=Iz*drdy*dAfidr-Iy*dAfidz;
		resulty[k]=Ix*dAfidz-Iz*drdx*dAfidr;
		resultz[k]=(Iy*drdx-Ix*drdy)*dAfidr;
	}
}



void SmoothOutput2D_Ar::Output_dArdz(double *_v2,double *result,int OutType,vector<int> &RecToSourceB)
{
	int i,j,k,m,it,nn[2];
	double s,ls,v,hr,hz,psi;
	int ipls;
	double *v2;

	for(k=0;k<nrec;k++)
	{
		ipls=RecToSourceB[k/2];
		v2=_v2+kpnt*ipls;
		if(OutType && OutRecList[k].oz.fl)
		{
			s=0.0;
			for(it=0;it<2;it++)
			{
				if((!it && OutRecList[k].oz.ves) || (it && OutRecList[k].oz.ves!=1.0))
				{
					ls=0.0;				
					for(j=0;j<3;j++)
					{
						m=OutRecList[k].oz.lem[j+it];

						nn[0]=rect[m].nodes[0]-1;
						nn[1]=rect[m].nodes[3]-1;

						hr=pnt[nn[1]].r-pnt[nn[0]].r;
						hz=pnt[nn[1]].z-pnt[nn[0]].z;
						psi=(rec[k].r-pnt[nn[0]].r)/hr;

						v=(-(1.0-psi)*v2[rect[m].nodes[0]-1]-psi*v2[rect[m].nodes[1]-1]
								+(1.0-psi)*v2[rect[m].nodes[2]-1]+psi*v2[rect[m].nodes[3]-1])/hz;

						ls+=OutRecList[k].oz.v[it][j]*v;
					}
					s+=ls*((it)? 1.0-OutRecList[k].oz.ves : OutRecList[k].oz.ves);
				}
			}
		}
		else
		{
			i=RecToElem[k];

			nn[0]=rect[i].nodes[0]-1;
			nn[1]=rect[i].nodes[3]-1;

			hr=pnt[nn[1]].r-pnt[nn[0]].r;
			hz=pnt[nn[1]].z-pnt[nn[0]].z;
			psi=(rec[k].r-pnt[nn[0]].r)/hr;

			s=(-(1.0-psi)*v2[rect[i].nodes[0]-1]-psi*v2[rect[i].nodes[1]-1]
				+(1.0-psi)*v2[rect[i].nodes[2]-1]+psi*v2[rect[i].nodes[3]-1])/hz;
		}

		result[k]=s;
	}
}

void SmoothOutput2D_Ar::Output_Er(double *_v2a,double *result,int OutType,vector<int> &RecToSourceE)
{
	int i,j,k,l,m,it,nn[2];
	double s,ls,v,hr,hz,psi,eta;
	int ipls;
	double *v2a;
	
	for(k=0;k<nrec;k++)
	{
		ipls=RecToSourceE[k/2];
		v2a=_v2a+kpnt*ipls;

		if(OutType && OutRecList[k].or.fl)
		{
			s=0.0;
			for(it=0;it<2;it++)
			{
				if((!it && OutRecList[k].or.ves) || (it && OutRecList[k].or.ves!=1.0))
				{					
					ls=0.0;
					for(j=0;j<3;j++)
					{
						m=OutRecList[k].or.lem[j+it];
						
						v=0.0;
						for(l=0;l<4;l++)
						{
							v+=OutRecList[k].or.coef[j+it][l]*v2a[rect[m].nodes[l]-1];
						}	
						
						ls+=OutRecList[k].or.v[it][j]*v;
					}

					s+=ls*((it)? 1.0-OutRecList[k].or.ves : OutRecList[k].or.ves);
				}
			}
		}
		else
		{
			s=0.0;
			i=RecToElem[k];

			nn[0]=rect[i].nodes[0]-1;
			nn[1]=rect[i].nodes[3]-1;

			hr=pnt[nn[1]].r-pnt[nn[0]].r;
			hz=pnt[nn[1]].z-pnt[nn[0]].z;
			psi=(rec[k].r-pnt[nn[0]].r)/hr;
			eta=(rec[k].z-pnt[nn[0]].z)/hz;

			s+=(1.0-psi)*(1.0-eta)*v2a[rect[i].nodes[0]-1]+psi*(1.0-eta)*v2a[rect[i].nodes[1]-1]
				+(1.0-psi)*eta*v2a[rect[i].nodes[2]-1]+psi*eta*v2a[rect[i].nodes[3]-1];
		}

		result[k]=s;
	}
}

void SmoothOutput2D_Ar::Output_Er_Pr(double *_v2a,double *result,int OutType,vector<int> &RecToSourceE,vector<int> &OutCompFlag)
{
	int i,j,k,l,m,it,nn[2];
	double s,ls,v,hr,hz,psi,eta;
	int ipls;
	double *v2a;

	#pragma omp parallel private(i,j,k,l,m,it,nn,s,ls,v,hr,hz,psi,eta,ipls,v2a)
	{
		#pragma omp for
		for(k=0;k<nrec;k++)
		{
			if(OutCompFlag[k/2])
			{
				ipls=RecToSourceE[k/2];
				v2a=_v2a+kpnt*ipls;

				if(OutType && OutRecList[k].or.fl)
				{
					s=0.0;
					for(it=0;it<2;it++)
					{
						if((!it && OutRecList[k].or.ves) || (it && OutRecList[k].or.ves!=1.0))
						{					
							ls=0.0;
							for(j=0;j<3;j++)
							{
								m=OutRecList[k].or.lem[j+it];
								
								v=0.0;
								for(l=0;l<4;l++)
								{
									v+=OutRecList[k].or.coef[j+it][l]*v2a[rect[m].nodes[l]-1];
								}	
								
								ls+=OutRecList[k].or.v[it][j]*v;
							}

							s+=ls*((it)? 1.0-OutRecList[k].or.ves : OutRecList[k].or.ves);
						}
					}
				}
				else
				{
					s=0.0;
					i=RecToElem[k];

					nn[0]=rect[i].nodes[0]-1;
					nn[1]=rect[i].nodes[3]-1;

					hr=pnt[nn[1]].r-pnt[nn[0]].r;
					hz=pnt[nn[1]].z-pnt[nn[0]].z;
					psi=(rec[k].r-pnt[nn[0]].r)/hr;
					eta=(rec[k].z-pnt[nn[0]].z)/hz;

					s+=(1.0-psi)*(1.0-eta)*v2a[rect[i].nodes[0]-1]+psi*(1.0-eta)*v2a[rect[i].nodes[1]-1]
						+(1.0-psi)*eta*v2a[rect[i].nodes[2]-1]+psi*eta*v2a[rect[i].nodes[3]-1];
				}

				result[k]=s;
			}
		}
	}
}


void SmoothOutput2D_Hfi::Output_Er(double *_v2,double *result,int OutType,vector<int> &RecToSourceE)
{
	int i,j,k,m,it,node0,node1,node2,node3;
	double s,ls,v,hr,hz,psi,eta;
	int ipls;
	double *v2;

	for(k=0;k<nrec;k++)
	{
		if(RecToElem[k]!=-1)
		{
		ipls=RecToSourceE[k/2];
		v2=_v2+4*krect*ipls;

		if(OutType && OutRecList[k].or.fl)
		{
			s=0.0;
			for(it=0;it<2;it++)
			{
				if((!it && OutRecList[k].or.ves) || (it && OutRecList[k].or.ves!=1.0))
				{					
					ls=0.0;
					for(j=0;j<3;j++)
					{
						m=OutRecList[k].or.lem[j+it];

						node0=rect[m].nodes[0]-1;
						node1=rect[m].nodes[1]-1;
						node2=rect[m].nodes[2]-1;
						node3=rect[m].nodes[3]-1;

						hr=pnt[node3].r-pnt[node0].r;
						hz=pnt[node3].z-pnt[node0].z;

						psi=0.5;
						
						v=GetSmoothValForR(k,j+it,v2,OutType);
						ls+=OutRecList[k].or.v[it][j]*v;
					}

					s+=ls*((it)? 1.0-OutRecList[k].or.ves : OutRecList[k].or.ves);
				}
			}
		}
		else
		{
			i=RecToElem[k];

			node0=rect[i].nodes[0]-1;
			node1=rect[i].nodes[1]-1;
			node2=rect[i].nodes[2]-1;
			node3=rect[i].nodes[3]-1;

			hr=pnt[node3].r-pnt[node0].r;
			hz=pnt[node3].z-pnt[node0].z;
			
			psi=(rec[k].r-pnt[node0].r)/hr;
			eta=(rec[k].z-pnt[node0].z)/hz;

			s=-(-(1.0-psi)*v2[4*i+0]-psi*v2[4*i+1]+(1.0-psi)*v2[4*i+2]+psi*v2[4*i+3])/hz;

			s/=sigma[rect[i].mtr-1];
		}

		result[k]=s;
	}
		else
		{
			result[k]=0.0;
		}
	}
}

void SmoothOutput2D_Hfi::Output_Ez(double *_v2,double *result,int OutType,vector<int> &RecToSourceE)
{
	int i,j,k,m,it,node0,node1,node2,node3;
	double s,ls,v,hr,hz,psi,eta,r;
	int ipls;
	double *v2;
	double sig;
	
	for(k=0;k<nrec;k++)
	{
		if(RecToElem[k]!=-1)
		{
		ipls=RecToSourceE[k/2];
		v2=_v2+4*krect*ipls;

		if(OutType && OutRecList[k].oz.fl)
		{
			s=0.0;
			for(it=0;it<2;it++)
			{
				if((!it && OutRecList[k].oz.ves) || (it && OutRecList[k].oz.ves!=1.0))
				{					
					ls=0.0;
					for(j=0;j<3;j++)
					{
						m=OutRecList[k].oz.lem[j+it];
						v=GetSmoothValForZ(k,j+it,v2,OutType);
						ls+=OutRecList[k].oz.v[it][j]*v;
					}
					s+=ls*((it)? 1.0-OutRecList[k].oz.ves : OutRecList[k].oz.ves);
				}
			}
		}
		else
		{
			i=RecToElem[k];

			node0=rect[i].nodes[0]-1;
			node1=rect[i].nodes[1]-1;
			node2=rect[i].nodes[2]-1;
			node3=rect[i].nodes[3]-1;

			hr=pnt[node3].r-pnt[node0].r;
			hz=pnt[node3].z-pnt[node0].z;
			
			psi=(rec[k].r-pnt[node0].r)/hr;
			eta=(rec[k].z-pnt[node0].z)/hz;

			r=rec[k].r;

			s=((1.0-psi)*(1.0-eta)*v2[4*i+0]+psi*(1.0-eta)*v2[4*i+1]+(1.0-psi)*eta*v2[4*i+2]+psi*eta*v2[4*i+3])/r+
				(-(1.0-eta)*v2[4*i+0]+(1.0-eta)*v2[4*i+1]-eta*v2[4*i+2]+eta*v2[4*i+3])/hr;
		}

		sig=sigmaZ[rect[RecToElem[k]].mtr-1];
		result[k]=s/sig;
	}
		else
		{
			result[k]=0.0;
		}
	}
}

void SmoothOutput2D_Hfi::Output_Er_Pr(double *_v2,double *result,int OutType,vector<int> &RecToSourceE,vector<int> &OutCompFlag)
{
	int i,j,k,m,it,node0,node1,node2,node3;
	double s,ls,v,hr,hz,psi,eta;
	int ipls;
	double *v2;

	#pragma omp parallel private(i,j,k,m,it,node0,node1,node2,node3,s,ls,v,hr,hz,psi,eta,ipls,v2)
	{
		#pragma omp for
		for(k=0;k<nrec;k++)
		{
			if(RecToElem[k]!=-1)
			{
			if(OutCompFlag[k/2])
			{
				ipls=RecToSourceE[k/2];
				v2=_v2+4*krect*ipls;

				if(OutType && OutRecList[k].or.fl)
				{
					s=0.0;
					for(it=0;it<2;it++)
					{
						if((!it && OutRecList[k].or.ves) || (it && OutRecList[k].or.ves!=1.0))
						{					
							ls=0.0;
							for(j=0;j<3;j++)
							{
								m=OutRecList[k].or.lem[j+it];

								node0=rect[m].nodes[0]-1;
								node1=rect[m].nodes[1]-1;
								node2=rect[m].nodes[2]-1;
								node3=rect[m].nodes[3]-1;

								hr=pnt[node3].r-pnt[node0].r;
								hz=pnt[node3].z-pnt[node0].z;

								psi=0.5;
								
								v=GetSmoothValForR(k,j+it,v2,OutType);
								ls+=OutRecList[k].or.v[it][j]*v;
							}

							s+=ls*((it)? 1.0-OutRecList[k].or.ves : OutRecList[k].or.ves);
						}
					}
				}
				else
				{
					i=RecToElem[k];

					node0=rect[i].nodes[0]-1;
					node1=rect[i].nodes[1]-1;
					node2=rect[i].nodes[2]-1;
					node3=rect[i].nodes[3]-1;

					hr=pnt[node3].r-pnt[node0].r;
					hz=pnt[node3].z-pnt[node0].z;
					
					psi=(rec[k].r-pnt[node0].r)/hr;
					eta=(rec[k].z-pnt[node0].z)/hz;

					s=-(-(1.0-psi)*v2[4*i+0]-psi*v2[4*i+1]+(1.0-psi)*v2[4*i+2]+psi*v2[4*i+3])/hz;

					s/=sigma[rect[i].mtr-1];
				}

				result[k]=s;
			}
		}
			else
			{
				result[k]=0.0;
			}
		}
	}
}

void SmoothOutput2D_Hfi::Output_Ez_Pr(double *_v2,double *result,int OutType,vector<int> &RecToSourceE,vector<int> &OutCompFlag)
{
	int i,j,k,m,it,node0,node1,node2,node3;
	double s,ls,v,hr,hz,psi,eta,r;
	int ipls;
	double *v2;
	double sig;

	#pragma omp parallel private(i,j,k,m,it,node0,node1,node2,node3,s,ls,v,hr,hz,psi,eta,ipls,v2,sig)
	{
		#pragma omp for
		for(k=0;k<nrec;k++)
		{
			if(RecToElem[k]!=-1)
			{
			if(OutCompFlag[k/2])
			{
				ipls=RecToSourceE[k/2];
				v2=_v2+4*krect*ipls;

				if(OutType && OutRecList[k].oz.fl)
				{
					s=0.0;
					for(it=0;it<2;it++)
					{
						if((!it && OutRecList[k].oz.ves) || (it && OutRecList[k].oz.ves!=1.0))
						{					
							ls=0.0;
							for(j=0;j<3;j++)
							{
								m=OutRecList[k].oz.lem[j+it];
								v=GetSmoothValForZ(k,j+it,v2,OutType);
								ls+=OutRecList[k].oz.v[it][j]*v;
							}
							s+=ls*((it)? 1.0-OutRecList[k].oz.ves : OutRecList[k].oz.ves);
						}
					}
				}
				else
				{
					i=RecToElem[k];

					node0=rect[i].nodes[0]-1;
					node1=rect[i].nodes[1]-1;
					node2=rect[i].nodes[2]-1;
					node3=rect[i].nodes[3]-1;

					hr=pnt[node3].r-pnt[node0].r;
					hz=pnt[node3].z-pnt[node0].z;
					
					psi=(rec[k].r-pnt[node0].r)/hr;
					eta=(rec[k].z-pnt[node0].z)/hz;

					r=rec[k].r;

					s=((1.0-psi)*(1.0-eta)*v2[4*i+0]+psi*(1.0-eta)*v2[4*i+1]+(1.0-psi)*eta*v2[4*i+2]+psi*eta*v2[4*i+3])/r+
						(-(1.0-eta)*v2[4*i+0]+(1.0-eta)*v2[4*i+1]-eta*v2[4*i+2]+eta*v2[4*i+3])/hr;
				}

				sig=sigmaZ[rect[RecToElem[k]].mtr-1];
				result[k]=s/sig;
			}
		}
			else
			{
				result[k]=0.0;
			}
		}
	}
}

void SmoothOutput2D_Hfi::Output_Bphi(double *_v2,double *result,int OutType,vector<int> &RecToSourceB)
{
	int i,j,k,l,m,it,nn[2];
	double s,ls,v,hr,hz,psi,eta,mup;
	int ipls;
	double *v2;
	
	for(k=0;k<nrec;k++)
	{
		if(RecToElem[k]!=-1)
		{
		ipls=RecToSourceB[k/2];
		v2=_v2+4*krect*ipls;

		if(OutType && OutRecList[k].or.fl)
		{
			s=0.0;
			for(it=0;it<2;it++)
			{
				if((!it && OutRecList[k].or.ves) || (it && OutRecList[k].or.ves!=1.0))
				{					
					ls=0.0;
					for(j=0;j<3;j++)
					{
						m=OutRecList[k].or.lem[j+it];
						
						v=0.0;
						for(l=0;l<4;l++)
						{
							v+=OutRecList[k].or.coef[j+it][l]*v2[4*m+l];
						}	
						
						ls+=OutRecList[k].or.v[it][j]*v;
					}

					s+=ls*((it)? 1.0-OutRecList[k].or.ves : OutRecList[k].or.ves);
				}
			}
		}
		else
		{
			s=0.0;
			i=RecToElem[k];

			nn[0]=rect[i].nodes[0]-1;
			nn[1]=rect[i].nodes[3]-1;

			hr=pnt[nn[1]].r-pnt[nn[0]].r;
			hz=pnt[nn[1]].z-pnt[nn[0]].z;
			psi=(rec[k].r-pnt[nn[0]].r)/hr;
			eta=(rec[k].z-pnt[nn[0]].z)/hz;

			s+=(1.0-psi)*(1.0-eta)*v2[4*i+0]+psi*(1.0-eta)*v2[4*i+1]+(1.0-psi)*eta*v2[4*i+2]+psi*eta*v2[4*i+3];
		}

		mup=mu[rect[RecToElem[k]].mtr-1];
		result[k]=mup*MU0*s;
	}
		else
		{
			result[k]=0.0;
		}
	}
}

double SmoothOutput2D_Hfi::GetSmoothValForZ(int k,int jj,double *v2,int OutType)
{
	int i,j,m,it,node0,node1,node2,node3;
	double s,ls,v,hr,hz,psi,eta,r;

	if(OutType || !OutRecList[k].orfz[jj].fl)
	{
		s=0.0;
		for(it=0;it<2;it++)
		{
			if((!it && OutRecList[k].orfz[jj].ves) || (it && OutRecList[k].orfz[jj].ves!=1.0))
			{
				ls=0.0;
				for(j=0;j<3;j++)
				{
					m=OutRecList[k].orfz[jj].lem[j+it];

					node0=rect[m].nodes[0]-1;
					node1=rect[m].nodes[1]-1;
					node2=rect[m].nodes[2]-1;
					node3=rect[m].nodes[3]-1;

					hr=pnt[node3].r-pnt[node0].r;
					hz=pnt[node3].z-pnt[node0].z;
					


					psi=eta=0.5;
					r=0.5*(pnt[node3].r+pnt[node0].r);

					v=((1.0-psi)*(1.0-eta)*v2[4*m+0]+psi*(1.0-eta)*v2[4*m+1]+(1.0-psi)*eta*v2[4*m+2]+psi*eta*v2[4*m+3])/r+
						(-(1.0-eta)*v2[4*m+0]+(1.0-eta)*v2[4*m+1]-eta*v2[4*m+2]+eta*v2[4*m+3])/hr;

					ls+=OutRecList[k].orfz[jj].v[it][j]*v;
				}
				s+=ls*((it)? 1.0-OutRecList[k].orfz[jj].ves : OutRecList[k].orfz[jj].ves);
			}
		}
	}
	else
	{
		i=OutRecList[k].oz.lem[jj];

		node0=rect[i].nodes[0]-1;
		node1=rect[i].nodes[1]-1;
		node2=rect[i].nodes[2]-1;
		node3=rect[i].nodes[3]-1;

		hr=pnt[node3].r-pnt[node0].r;
		hz=pnt[node3].z-pnt[node0].z;
		


		psi=(rec[k].r-pnt[node0].r)/hr;
		eta=0.5;
		r=rec[k].r;

		s=((1.0-psi)*(1.0-eta)*v2[4*i+0]+psi*(1.0-eta)*v2[4*i+1]+(1.0-psi)*eta*v2[4*i+2]+psi*eta*v2[4*i+3])/r+
			(-(1.0-eta)*v2[4*i+0]+(1.0-eta)*v2[4*i+1]-eta*v2[4*i+2]+eta*v2[4*i+3])/hr;
	}

	return s;
}

double SmoothOutput2D_Hfi::GetSmoothValForR(int k,int jj,double *v2,int OutType)
{
	int i,j,m,it,node0,node1,node2,node3;
	double s,ls,v,hr,hz,psi,eta;

	if(OutType && OutRecList[k].ozfr[jj].fl)
	{
		s=0.0;
		for(it=0;it<2;it++)
		{
			if((!it && OutRecList[k].ozfr[jj].ves) || (it && OutRecList[k].ozfr[jj].ves!=1.0))
			{
				ls=0.0;
				for(j=0;j<3;j++)
				{
					m=OutRecList[k].ozfr[jj].lem[j+it];

					node0=rect[m].nodes[0]-1;
					node1=rect[m].nodes[1]-1;
					node2=rect[m].nodes[2]-1;
					node3=rect[m].nodes[3]-1;
					
					hr=pnt[node3].r-pnt[node0].r;
					hz=pnt[node3].z-pnt[node0].z;
					

					psi=eta=0.5;

					v=-(-(1.0-psi)*v2[4*m+0]-psi*v2[4*m+1]+(1.0-psi)*v2[4*m+2]+psi*v2[4*m+3])/hz;

					ls+=(OutRecList[k].ozfr[jj].v[it][j]*v)/sigma[rect[m].mtr-1];
				}
				s+=ls*((it)? 1.0-OutRecList[k].ozfr[jj].ves : OutRecList[k].ozfr[jj].ves);
			}
		}
	}
	else
	{
		i=OutRecList[k].or.lem[jj];

		node0=rect[i].nodes[0]-1;
		node1=rect[i].nodes[1]-1;
		node2=rect[i].nodes[2]-1;
		node3=rect[i].nodes[3]-1;
		
		hr=pnt[node3].r-pnt[node0].r;
		hz=pnt[node3].z-pnt[node0].z;


		psi=0.5;
		eta=(rec[k].z-pnt[node0].z)/hz;

		s=-(-(1.0-psi)*v2[4*i+0]-psi*v2[4*i+1]+(1.0-psi)*v2[4*i+2]+psi*v2[4*i+3])/hz;

		s/=sigma[rect[i].mtr-1];
	}

	return s;
}
