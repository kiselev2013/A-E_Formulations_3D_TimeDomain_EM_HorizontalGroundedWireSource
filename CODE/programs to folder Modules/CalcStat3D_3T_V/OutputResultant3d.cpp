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
 *  This file contains the code of functions for outputting with smoothing in 3D tasks
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 January 16, 2023
*/

#include "stdafx.h"
#include "T_Brick.h"
#include "OutputResultant3d.h"
#include "Portret.h"

extern ofstream logfile;

double phc_scal(pv::Point3D &p1,pv::Point3D &p2) { return p1.x()*p2.x()+p1.y()*p2.y()+p1.z()*p2.z(); }

int phc_sign(double x) { if(fabs(x)<1e-6)return 0;if(x>0)return 1;else return -1;}

bool phc_check(_Plane_ &cPlane,pv::Point3D &cPoint0) {if(phc_sign(phc_scal(cPlane.N,cPoint0)+cPlane.D)==1)return false;else return true;}

int phc_main(const pv::Point3D *P,pv::Point3D pM)
{
	bool flag;
	_Plane_ pl[4];
	pv::Point3D HEXP[8],M;
	pv::Vector V[4];
	int igr,ii,ind[2];

	M.x()=pM.getx();
	M.y()=pM.gety();
	M.z()=pM.getz();

	for(igr=0;igr<6;igr++){

		switch(igr){
			case 0:{HEXP[0]=P[0];HEXP[1]=P[1];HEXP[2]=P[2];HEXP[3]=P[3];HEXP[4]=P[4];HEXP[5]=P[5];HEXP[6]=P[6];HEXP[7]=P[7];break;}
		    case 1:{HEXP[0]=P[4];HEXP[1]=P[5];HEXP[2]=P[6];HEXP[3]=P[7];HEXP[4]=P[0];HEXP[5]=P[1];HEXP[6]=P[2];HEXP[7]=P[3];break;}
		    case 2:{HEXP[0]=P[0];HEXP[1]=P[2];HEXP[2]=P[4];HEXP[3]=P[6];HEXP[4]=P[1];HEXP[5]=P[3];HEXP[6]=P[5];HEXP[7]=P[7];break;}
		    case 3:{HEXP[0]=P[1];HEXP[1]=P[3];HEXP[2]=P[5];HEXP[3]=P[7];HEXP[4]=P[0];HEXP[5]=P[2];HEXP[6]=P[4];HEXP[7]=P[6];break;}
		    case 4:{HEXP[0]=P[0];HEXP[1]=P[4];HEXP[2]=P[1];HEXP[3]=P[5];HEXP[4]=P[2];HEXP[5]=P[6];HEXP[6]=P[3];HEXP[7]=P[7];break;}
	        case 5:{HEXP[0]=P[2];HEXP[1]=P[6];HEXP[2]=P[3];HEXP[3]=P[7];HEXP[4]=P[0];HEXP[5]=P[4];HEXP[6]=P[1];HEXP[7]=P[5];break;}
		}

		V[0]=HEXP[1]-HEXP[0];
		V[1]=HEXP[3]-HEXP[2];
		V[2]=HEXP[2]-HEXP[0];
		V[3]=HEXP[3]-HEXP[1];
		
		pl[0].N= (igr%2) ? V[0].vp(V[2]) : V[2].vp(V[0]);
		pl[0].set_D(HEXP[0]);
		
		if(!phc_sign(phc_scal(pl[0].N,HEXP[3])+pl[0].D)){
			flag=phc_check(pl[0],M);
		}
		else{
			ii=0;
			if(phc_check(pl[0],HEXP[3])){ind[ii]=0;ii++;}
			pl[1].N= (igr%2) ? V[0].vp(V[3]) : V[3].vp(V[0]);
			pl[1].set_D(HEXP[1]);
			if(phc_check(pl[1],HEXP[2])){ind[ii]=1;ii++;}
			if(ii<2){
				pl[2].N= (igr%2) ? V[1].vp(V[2]) : V[2].vp(V[1]);
				pl[2].set_D(HEXP[2]);
				if(phc_check(pl[2],HEXP[1])){ind[ii]=2;ii++;}
				if(ii<2){
					pl[3].N= (igr%2) ? V[1].vp(V[3]) : V[3].vp(V[1]);
					pl[3].set_D(HEXP[3]);
					if(phc_check(pl[3],HEXP[0])){ind[ii]=3;ii++;}
				}
			}
			if(ii!=2)throw logic_error("Point in Hex problem!");
			flag=phc_check(pl[ind[0]],M) && phc_check(pl[ind[1]],M);
		}
		if(!flag)break;
	}
	return flag;
}

void GetGlobalCoordinates(pv::Point3D *HexPnt,double *lc,double *gc)
{
	double ksi,eta,phi;

	ksi=lc[0];
	eta=lc[1];
	phi=lc[2];

	gc[0]=	HexPnt[0].x()*(1-ksi)*(1-eta)*(1-phi)+HexPnt[1].x()*(ksi)*(1-eta)*(1-phi)+
			HexPnt[2].x()*(1-ksi)*(eta)*(1-phi)+HexPnt[3].x()*(ksi)*(eta)*(1-phi)+
			HexPnt[4].x()*(1-ksi)*(1-eta)*(phi)+HexPnt[5].x()*(ksi)*(1-eta)*(phi)+
			HexPnt[6].x()*(1-ksi)*(eta)*(phi)+HexPnt[7].x()*(ksi)*(eta)*(phi);
	gc[1]=	HexPnt[0].y()*(1-ksi)*(1-eta)*(1-phi)+HexPnt[1].y()*(ksi)*(1-eta)*(1-phi)+
			HexPnt[2].y()*(1-ksi)*(eta)*(1-phi)+HexPnt[3].y()*(ksi)*(eta)*(1-phi)+
			HexPnt[4].y()*(1-ksi)*(1-eta)*(phi)+HexPnt[5].y()*(ksi)*(1-eta)*(phi)+
			HexPnt[6].y()*(1-ksi)*(eta)*(phi)+HexPnt[7].y()*(ksi)*(eta)*(phi);
	gc[2]=	HexPnt[0].z()*(1-ksi)*(1-eta)*(1-phi)+HexPnt[1].z()*(ksi)*(1-eta)*(1-phi)+
			HexPnt[2].z()*(1-ksi)*(eta)*(1-phi)+HexPnt[3].z()*(ksi)*(eta)*(1-phi)+
			HexPnt[4].z()*(1-ksi)*(1-eta)*(phi)+HexPnt[5].z()*(ksi)*(1-eta)*(phi)+
			HexPnt[6].z()*(1-ksi)*(eta)*(phi)+HexPnt[7].z()*(ksi)*(eta)*(phi);

}

void FindLocalCoordinates(pv::Point3D &R,pv::Point3D *HexPnt,double *lc)
{
	double EpsForFindLocalCoord = 1e-2;
	double CoeffForDiv = 1.2;
	int MaxDeep = 10;


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
		(R.x()-CentGlob[0])*(R.x()-CentGlob[0])+
		(R.y()-CentGlob[1])*(R.y()-CentGlob[1])+
		(R.z()-CentGlob[2])*(R.z()-CentGlob[2])
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
					(R.x()-CentGlob[0])*(R.x()-CentGlob[0])+
					(R.y()-CentGlob[1])*(R.y()-CentGlob[1])+
					(R.z()-CentGlob[2])*(R.z()-CentGlob[2])
					);

			if(disc<dist){
				CentLoc[m]=LocalCoord[1][m]-ods*h[m];
			
				GetGlobalCoordinates(HexPnt,CentLoc,CentGlob);
				dist=sqrt(
						(R.x()-CentGlob[0])*(R.x()-CentGlob[0])+
						(R.y()-CentGlob[1])*(R.y()-CentGlob[1])+
						(R.z()-CentGlob[2])*(R.z()-CentGlob[2])
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
								(R.x()-CentGlob[0])*(R.x()-CentGlob[0])+
								(R.y()-CentGlob[1])*(R.y()-CentGlob[1])+
								(R.z()-CentGlob[2])*(R.z()-CentGlob[2])
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
								(R.x()-CentGlob[0])*(R.x()-CentGlob[0])+
								(R.y()-CentGlob[1])*(R.y()-CentGlob[1])+
								(R.z()-CentGlob[2])*(R.z()-CentGlob[2])
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
		double norm_R,norm_F;
		norm_R=sqrt(R.x()*R.x()+R.y()*R.y()+R.z()*R.z());
		norm_F=sqrt((R.x()-CentGlob[0])*(R.x()-CentGlob[0])+
				   (R.y()-CentGlob[1])*(R.y()-CentGlob[1])+
				   (R.z()-CentGlob[2])*(R.z()-CentGlob[2]));
		logfile<<"Reciver Eps Fail:"<<'\n';
		logfile<<R.x()<<'\t'<<R.y()<<'\t'<<R.z()<<'\n';
		logfile<<CentGlob[0]<<'\t'<<CentGlob[1]<<'\t'<<CentGlob[2]<<'\n';
		logfile<<norm_R<<'\n';
		logfile<<norm_F<<'\n';
		logfile<<endl;
	}
	

	lc[0]=CentLoc[0];
	lc[1]=CentLoc[1];
	lc[2]=CentLoc[2];
}

bool CheckPointInHex(pv::Point3D &R,pv::Point3D *HexPnt,double &mdist)
{
	bool mnf;
	double mnv;
	double EpsForFindLocalCoord = 1e-2;
	int MaxDeep = 9;

	int p,t,m,deep,crd[3][2]={{1,2},{0,2},{0,1}};
	double LocalCoord[2][3],CentGlob[3],CentLoc[3],DopLoc[3];
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
		(R.x()-CentGlob[0])*(R.x()-CentGlob[0])+
		(R.y()-CentGlob[1])*(R.y()-CentGlob[1])+
		(R.z()-CentGlob[2])*(R.z()-CentGlob[2])
		);
	h[0]=LocalCoord[1][0]-LocalCoord[0][0];
	h[1]=LocalCoord[1][1]-LocalCoord[0][1];
	h[2]=LocalCoord[1][2]-LocalCoord[0][2];

	ddg[0]=sqrt((HexPnt[7].x()-HexPnt[0].x())*(HexPnt[7].x()-HexPnt[0].x())+(HexPnt[7].y()-HexPnt[0].y())*(HexPnt[7].y()-HexPnt[0].y())+(HexPnt[7].z()-HexPnt[0].z())*(HexPnt[7].z()-HexPnt[0].z()));
	ddg[1]=sqrt((HexPnt[6].x()-HexPnt[1].x())*(HexPnt[6].x()-HexPnt[1].x())+(HexPnt[6].y()-HexPnt[1].y())*(HexPnt[6].y()-HexPnt[1].y())+(HexPnt[6].z()-HexPnt[1].z())*(HexPnt[6].z()-HexPnt[1].z()));
	ddg[2]=sqrt((HexPnt[5].x()-HexPnt[2].x())*(HexPnt[5].x()-HexPnt[2].x())+(HexPnt[5].y()-HexPnt[2].y())*(HexPnt[5].y()-HexPnt[2].y())+(HexPnt[5].z()-HexPnt[2].z())*(HexPnt[5].z()-HexPnt[2].z()));
	ddg[3]=sqrt((HexPnt[4].x()-HexPnt[3].x())*(HexPnt[4].x()-HexPnt[3].x())+(HexPnt[4].y()-HexPnt[3].y())*(HexPnt[4].y()-HexPnt[3].y())+(HexPnt[4].z()-HexPnt[3].z())*(HexPnt[4].z()-HexPnt[3].z()));
	
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
					(R.x()-CentGlob[0])*(R.x()-CentGlob[0])+
					(R.y()-CentGlob[1])*(R.y()-CentGlob[1])+
					(R.z()-CentGlob[2])*(R.z()-CentGlob[2])
					);

			if(disc<dist)
			{
				CentLoc[m]=LocalCoord[1][m]-ods*h[m];
			
				GetGlobalCoordinates(HexPnt,CentLoc,CentGlob);
				dist=sqrt(
						(R.x()-CentGlob[0])*(R.x()-CentGlob[0])+
						(R.y()-CentGlob[1])*(R.y()-CentGlob[1])+
						(R.z()-CentGlob[2])*(R.z()-CentGlob[2])
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

void GetSolutionOnHex(pv::Point3D &R,pv::Point3D *CntHex,double *q,double &Val)
{
	double lc[3];

	FindLocalCoordinates(R,CntHex,lc);

	Val=q[0]*(1-lc[0])*(1-lc[1])*(1-lc[2])+q[1]*(lc[0])*(1-lc[1])*(1-lc[2])+
		q[2]*(1-lc[0])*(lc[1])*(1-lc[2])+q[3]*(lc[0])*(lc[1])*(1-lc[2])+
		q[4]*(1-lc[0])*(1-lc[1])*(lc[2])+q[5]*(lc[0])*(1-lc[1])*(lc[2])+
		q[6]*(1-lc[0])*(lc[1])*(lc[2])+q[7]*(lc[0])*(lc[1])*(lc[2]);
}

OutputResultant3d::OutputResultant3d(AbstractFEM3D *pTaskCalcMesh,const Res3DValueType& r_type)
{
	int i,n_neib;
	ifstream inf;
	pv::Point3D TmpPoint;
	TaskCalcMesh=pTaskCalcMesh;
	ValueType=r_type;
	n_pointres=TaskCalcMesh->GetNumberOfResPoints(ValueType);
	pointres = new double[n_pointres][3];
	for(i=0;i<n_pointres;i++){
		TmpPoint=TaskCalcMesh->GetResPoint(ValueType,i);
		pointres[i][0]=TmpPoint.x();
		pointres[i][1]=TmpPoint.y();
		pointres[i][2]=TmpPoint.z();
	}
	n_neib=5;
	inf.open("n_neib");
	if(inf)
	{
		inf>>n_neib;
		inf.close();
	}
	inf.clear();
	SetLevelNeighbors(n_neib);
}
OutputResultant3d::~OutputResultant3d()
{
	if(pointres) {delete [] pointres; pointres=NULL;}
}
void OutputResultant3d::SetLevelNeighbors(int val)
{
	levelNeighbors = val;
}
int OutputResultant3d::Prepare()
{
	PointresSort();
	FindPointsForElems();
	InitSubdomains();

	return 0;
}
int OutputResultant3d::InitSubdomains()
{
	int i;
	vector<int> mat;      // mat[i] ==    nvkat,   i-  
	int n_sub;

	int emat;
	int kpar=TaskCalcMesh->GetNumberOfElements();
	
	mat.clear();

	for (i=0; i<kpar; i++)
	{
		if (PointsForElem[i].size() > 0) //   ,    ...
		{
			emat=TaskCalcMesh->GetElementMaterial(i);

			if (!binary_search(mat.begin(), mat.end(), emat)) 
			{
				mat.push_back(emat);
				sort(mat.begin(), mat.end());
			}
		}
	}

	if(ValueType==vtWithoutDiscontinuity)
	{
		n_sub = 1;
		sub.resize(n_sub);
		sub[0].Init(-1, levelNeighbors, PointsForElem, TaskCalcMesh);
	}
	else
	{
		n_sub = (int)mat.size();
		sub.resize(n_sub);
		for (i=0; i<n_sub; i++)
		{
			sub[i].Init(mat[i], levelNeighbors, PointsForElem, TaskCalcMesh);
		}
	}

	return 0;
}
int OutputResultant3d::Output(int itime)
{
	int i, j, k, n_sub;
	double var;

	n_sub=(int)sub.size();

	for (k=0; k<n_sub; k++)
	{
		Subdomain &mainSub=sub[k];

		for(j=0;j<mainSub.n_elem;j++)
		{
			mainSub.ValueInCenter[j]=
			TaskCalcMesh->GetValueInElementCenter(mainSub.renumElemFromNewToOld[j],
			ValueType);
		}
	
		mainSub.AsmGlobalVector();
	
		for (i=0; i<mainSub.n_nodes_c; i++)
			mainSub.x[i] = 0;

		mainSub.prds.solve_nrhs(1,mainSub.pr,mainSub.x);

		for (i=0; i<mainSub.n_elem; i++)
		{
			T_Brick L(i, mainSub.nver, mainSub.xyz);
			double ves[8];

			for (j=0; j<8; j++)
				ves[j] = mainSub.x[mainSub.nver[i][j]];

			int elem = mainSub.renumElemFromNewToOld[i];
			int sz = (int)PointsForElem[elem].size();

			pv::Point3D R;

			for (j=0; j<sz; j++)
			{
				int pnt = PointsForElem[elem][j];
				R.x()=pointres[pnt][0];
				R.y()=pointres[pnt][1];
				R.z()=pointres[pnt][2];
				var = L.ScalarFieldOnPar(R.x(), R.y(), R.z(), ves);
				TaskCalcMesh->SaveResult(ValueType,var,pnt,itime);
			}
		}
	}

	return 0;
}
int OutputResultant3d::OutputDirect(int itime)
{
	int i, j, k, n_sub;
	double var;

	n_sub=(int)sub.size();

	for(k=0;k<n_sub;k++)
	{
		Subdomain &mainSub=sub[k];

		for (i=0; i<mainSub.n_elem; i++)
		{
			T_Brick L(i, mainSub.nver, mainSub.xyz);
			double ves[8];

			TaskCalcMesh->GetValueInElementNodes(mainSub.renumElemFromNewToOld[i],ValueType,ves);

			int elem = mainSub.renumElemFromNewToOld[i];
			int sz = (int)PointsForElem[elem].size();

			pv::Point3D R;

			for (j=0; j<sz; j++)
			{
				int pnt = PointsForElem[elem][j];
				R.x()=pointres[pnt][0];
				R.y()=pointres[pnt][1];
				R.z()=pointres[pnt][2];
				var = L.ScalarFieldOnPar(R.x(), R.y(), R.z(), ves);
				TaskCalcMesh->SaveResult(ValueType,var,pnt,itime);
			}

		}
	}

	return 0;
}
void OutputResultant3d::PointresSort()
{
	vector<PointRes2> p;
	int i;


	logfile << "xyz-Sorting of " << n_pointres << " receiver(s) start\n";

	p.resize(n_pointres);

	for (i=0; i<n_pointres; i++)
	{
		p[i].num = i;
		p[i].point[0] = pointres[i][0];
		p[i].point[1] = pointres[i][1];
		p[i].point[2] = pointres[i][2];
	}
	sort(p.begin(), p.end(), PointRes_less_x());

	PointresXsorted.resize(n_pointres);
	for (i=0; i<n_pointres; i++)
	{
		PointresXsorted[i].i = p[i].num;
		PointresXsorted[i].d = p[i].point[0];
	}

	for (i=0; i<n_pointres; i++)
	{
		p[i].num = i;
		p[i].point[0] = pointres[i][0];
		p[i].point[1] = pointres[i][1];
		p[i].point[2] = pointres[i][2];
	}
	sort(p.begin(), p.end(), PointRes_less_y());
	PointresYsorted.resize(n_pointres);
	for (i=0; i<n_pointres; i++)
	{
		PointresYsorted[i].i = p[i].num;
		PointresYsorted[i].d = p[i].point[1];
	}

	for (i=0; i<n_pointres; i++)
	{
		p[i].num = i;
		p[i].point[0] = pointres[i][0];
		p[i].point[1] = pointres[i][1];
		p[i].point[2] = pointres[i][2];
	}
	sort(p.begin(), p.end(), PointRes_less_z());
	PointresZsorted.resize(n_pointres);
	for (i=0; i<n_pointres; i++)
	{
		PointresZsorted[i].i = p[i].num;
		PointresZsorted[i].d = p[i].point[2];
	}
}
int OutputResultant3d::FindPointsForElems()
{
	int i, j, sz, t;
	int typeOfHex;

	vector<long_double>::iterator beg_x, beg_y, beg_z, end_x, end_y, end_z; 

	long_double coordMin[3], coordMax[3]; 

	vector<long_double> pntInBrick, pntInBrickTmp;
	vector<long_double> tmp_x, tmp_y, tmp_z;


	int kpar=TaskCalcMesh->GetNumberOfElements();

	logfile << "Start searching elements for receivers. ";
	logfile << "Elements=" << kpar << " receivers=" << n_pointres << endl;
	fflush(stdout);


	PointsForElem.resize(kpar);

	ElemForPoint.resize(n_pointres);
	for(i=0; i<n_pointres; i++)
		ElemForPoint[i] = -1;

	double mdist;
	vector<int> edist;
	vector<double> vdist;

	edist.resize(n_pointres);
	vdist.resize(n_pointres);

	for(i=0;i<n_pointres;i++)
	{
		edist[i]=-1;
		vdist[i]=1e+30;
	}

	double recbrd[3][2];
	for(j=0;j<3;j++)
	{
		recbrd[j][0]=recbrd[j][1]=0.0;
	}
	if(n_pointres)
	{
		for(j=0;j<3;j++)
		{
			recbrd[j][0]=recbrd[j][1]=pointres[0][j];
		}
	}
	for(i=1;i<n_pointres;i++)
	{
		for(j=0;j<3;j++)
		{
			if(pointres[i][j]<recbrd[j][0])recbrd[j][0]=pointres[i][j];
			if(pointres[i][j]>recbrd[j][1])recbrd[j][1]=pointres[i][j];
		}
	}
	for(j=0;j<3;j++)
	{
		recbrd[j][0]-=1e-3;
		recbrd[j][1]+=1e-3;
	}

	for(i=0; i<kpar; i++)
	{
		typeOfHex = TaskCalcMesh->GetTypeOfElement(i);

		if(typeOfHex <= 30) // 
		{
			pv::Point3D TmpPoint;

			TmpPoint=TaskCalcMesh->GetNodeTrue(TaskCalcMesh->GetNodeNumberOnElement(i,0));
			coordMin[0].d=TmpPoint.x()-1e-6;
			coordMin[1].d=TmpPoint.y()-1e-6;
			coordMin[2].d=TmpPoint.z()-1e-6;
			TmpPoint=TaskCalcMesh->GetNodeTrue(TaskCalcMesh->GetNodeNumberOnElement(i,7));
			coordMax[0].d=TmpPoint.x()+1e-6;
			coordMax[1].d=TmpPoint.y()+1e-6;
			coordMax[2].d=TmpPoint.z()+1e-6;


			beg_x = lower_bound(PointresXsorted.begin(), PointresXsorted.end(), coordMin[0], Long_double_less());
			if(beg_x == PointresXsorted.end())
				continue;

			beg_y = lower_bound(PointresYsorted.begin(), PointresYsorted.end(), coordMin[1], Long_double_less());
			if(beg_y == PointresYsorted.end())
				continue;

			beg_z = lower_bound(PointresZsorted.begin(), PointresZsorted.end(), coordMin[2], Long_double_less());
			if(beg_z == PointresZsorted.end())
				continue;

			end_x = upper_bound(PointresXsorted.begin(), PointresXsorted.end(), coordMax[0], Long_double_less());
			if(end_x == PointresXsorted.begin())
				continue;

			end_y = upper_bound(PointresYsorted.begin(), PointresYsorted.end(), coordMax[1], Long_double_less());
			if(end_y == PointresYsorted.begin())
				continue;

			end_z = upper_bound(PointresZsorted.begin(), PointresZsorted.end(), coordMax[2], Long_double_less());
			if(end_z == PointresZsorted.begin())
				continue;

			if(distance(beg_x, end_x) <= 0)
				continue;

			if(distance(beg_y, end_y) <= 0)
				continue;

			if(distance(beg_z, end_z) <= 0)
				continue;

			tmp_x.clear();
			copy(beg_x, end_x, back_inserter(tmp_x));
			sort(tmp_x.begin(), tmp_x.end(), Long_double_num_less());

			tmp_y.clear();
			copy(beg_y, end_y, back_inserter(tmp_y));
			sort(tmp_y.begin(), tmp_y.end(), Long_double_num_less());

			tmp_z.clear();
			copy(beg_z, end_z, back_inserter(tmp_z));
			sort(tmp_z.begin(), tmp_z.end(), Long_double_num_less());

			set_intersection(tmp_x.begin(), tmp_x.end(), tmp_y.begin(), tmp_y.end(), back_inserter(pntInBrickTmp), Long_double_num_less()); 
			set_intersection(tmp_z.begin(), tmp_z.end(), pntInBrickTmp.begin(), pntInBrickTmp.end(), back_inserter(pntInBrick), Long_double_num_less()); 

			sz = (int)pntInBrick.size();
			for (j=0; j<sz; j++)
			{
				t = pntInBrick[j].i;

				if (ElemForPoint[t] == -1)
				{
					ElemForPoint[t] = i;
					PointsForElem[i].push_back(t);
				}
			}

			pntInBrick.clear();
			pntInBrickTmp.clear();
		}
		else // 
		{
			double hexbrd[3][2],pc[3];

			pv::Point3D P[8],pM;

			for(j=0;j<8;j++)P[j]=TaskCalcMesh->GetNodeTrue(TaskCalcMesh->GetNodeNumberOnElement(i,j));


			hexbrd[0][0] = hexbrd[0][1] = P[0].getx();
			hexbrd[1][0] = hexbrd[1][1] = P[0].gety();
			hexbrd[2][0] = hexbrd[2][1] = P[0].getz();

			for(t=1;t<8;t++)
			{
				pc[0]=P[t].getx();
				pc[1]=P[t].gety();
				pc[2]=P[t].getz();

				if(pc[0]<hexbrd[0][0])hexbrd[0][0]=pc[0];
				if(pc[0]>hexbrd[0][1])hexbrd[0][1]=pc[0];

				if(pc[1]<hexbrd[1][0])hexbrd[1][0]=pc[1];
				if(pc[1]>hexbrd[1][1])hexbrd[1][1]=pc[1];

				if(pc[2]<hexbrd[2][0])hexbrd[2][0]=pc[2];
				if(pc[2]>hexbrd[2][1])hexbrd[2][1]=pc[2];
			}

			for(j=0;j<3;j++)
			{
				hexbrd[j][0]-=0.5;
				hexbrd[j][1]+=0.5;
			}

			if( !(hexbrd[0][1]<recbrd[0][0] || hexbrd[0][0]>recbrd[0][1] ||
				hexbrd[1][1]<recbrd[1][0] || hexbrd[1][0]>recbrd[1][1] ||
				hexbrd[2][1]<recbrd[2][0] || hexbrd[2][0]>recbrd[2][1]) &&
				!(recbrd[0][1]<hexbrd[0][0] || recbrd[0][0]>hexbrd[0][1] ||
				recbrd[1][1]<hexbrd[1][0] || recbrd[1][0]>hexbrd[1][1] ||
				recbrd[2][1]<hexbrd[2][0] || recbrd[2][0]>hexbrd[2][1])	)
			{
				for(t=0; t<n_pointres; t++)
				{
					pM.x()=pointres[t][0];
					pM.y()=pointres[t][1];
					pM.z()=pointres[t][2];


					if (ElemForPoint[t] == -1)
					{
						if(CheckPointInHex(pM,P,mdist))
						{
							ElemForPoint[t] = i;
							PointsForElem[i].push_back(t);
						}
						else
						{
							if(edist[t]==-1 || vdist[t]>mdist)
							{
								edist[t]=i;
								vdist[t]=mdist;
							}
						}
					}
				}
			}
		}
	}

	fflush(stdout);

	{

		char str[256],buf[256],resList[256];

		bool flag = false;

		for (t=0; t<n_pointres; t++)
		{
			if (ElemForPoint[t] == -1)
			{

				i=edist[t];
				if(i==-1)
				{
					sprintf(buf,"%d",t);
					if(strlen(resList)+strlen(buf)<256)
					{
						strcat(resList,buf);
					}
					flag = true;
				}
				else
				{
					ElemForPoint[t]=i;
					PointsForElem[i].push_back(t);
				}
			}
		}

		if (flag)
		{
			sprintf(str,"There are no elements for receiver(s) N");
			strcat(str,resList);
			logfile<<str<<endl;
			exit(1);
		}
	}

	return 0;
}