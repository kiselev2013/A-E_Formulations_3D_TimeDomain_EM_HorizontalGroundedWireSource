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
 *  This file contains code for outputting stationary 2D task results
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 January 16, 2023
*/

#include"stdafx.h"
#include"mesh2d.h"
#include"mesh3d.h"
#include"retcode.h"

ofstream logfile;

/*!       2d  */
int ProcessResCalcPoints(MeshForVP* ptrMesh)
{
	int i, j;
	
	for (i=0; i<ptrMesh->ResCalcRectNum; i++)
	{
		RectForResultCalculation& crect=ptrMesh->ResCalcRectArray[i];
		const Rect& elrect=ptrMesh->rect[crect.rnum];
		for (j=0; j<4; j++)
			crect.q[j]=ptrMesh->q_all[elrect.nodes[j]-1];
	}


	for (i=ptrMesh->boundsVA.beginP; i<ptrMesh->boundsVA.endP; i++)
	{
		PointForResultCalculation& cpnt=ptrMesh->ResCalcPointArray[i];
		if (cpnt.rnum==-1)
		{
			cpnt.res->x=0;
			continue;
		}
		const RectForResultCalculation& crect=ptrMesh->ResCalcRectArray[cpnt.rnum];
		cpnt.res->x=(ptrMesh->meshoptions&moForCED||ptrMesh->meshoptions&moForVEL||ptrMesh->meshoptions&moForPED?1:-1)*(crect.q[0]*cpnt.rb*cpnt.zb+ 
					 crect.q[1]*cpnt.ra*cpnt.zb+
					 crect.q[2]*cpnt.rb*cpnt.za+
					 crect.q[3]*cpnt.ra*cpnt.za)/(crect.hr*crect.hz);
	}


	for (i=ptrMesh->boundsVB.beginP; i<ptrMesh->boundsVB.endP; i++)
	{
		PointForResultCalculation& cpnt=ptrMesh->ResCalcPointArray[i];
		if (ptrMesh->meshoptions&moForVEL||ptrMesh->meshoptions&moForCED)
		{
			cpnt.res->x=0;
			continue;
		}
		if (cpnt.rnum==-1)
		{
			cpnt.res->x=0;
			continue;
		}
		const RectForResultCalculation& crect=ptrMesh->ResCalcRectArray[cpnt.rnum];
		cpnt.res->x=(crect.q[0]*cpnt.rb*cpnt.zb+ 
					  crect.q[1]*cpnt.ra*cpnt.zb+
					  crect.q[2]*cpnt.rb*cpnt.za+
					  crect.q[3]*cpnt.ra*cpnt.za)/(crect.hr*crect.hz);
	}
	for (i=ptrMesh->boundsVforMNfromA.beginP; i<ptrMesh->boundsVforMNfromA.endP; i++)
	{
		PointForResultCalculation& cpnt=ptrMesh->ResCalcPointArray[i];
		if (cpnt.rnum==-1)
		{
			cpnt.res->x=0;
			continue;
		}
		const RectForResultCalculation& crect=ptrMesh->ResCalcRectArray[cpnt.rnum];
		cpnt.res->x=(ptrMesh->meshoptions&moForCED||ptrMesh->meshoptions&moForVEL?1:-1)*(crect.q[0]*cpnt.rb*cpnt.zb+ 
					crect.q[1]*cpnt.ra*cpnt.zb+
					crect.q[2]*cpnt.rb*cpnt.za+
					crect.q[3]*cpnt.ra*cpnt.za)/(crect.hr*crect.hz);
	}

	for (i=ptrMesh->boundsVforMNfromB.beginP; i<ptrMesh->boundsVforMNfromB.endP; i++)
	{
		PointForResultCalculation& cpnt=ptrMesh->ResCalcPointArray[i];
		if (ptrMesh->meshoptions&moForVEL||ptrMesh->meshoptions&moForCED)
		{
			cpnt.res->x=0;
			continue;
		}
		if (cpnt.rnum==-1)
		{
			cpnt.res->x=0;
			continue;
		}
		const RectForResultCalculation& crect=ptrMesh->ResCalcRectArray[cpnt.rnum];
		cpnt.res->x=(crect.q[0]*cpnt.rb*cpnt.zb+ 
					crect.q[1]*cpnt.ra*cpnt.zb+
					crect.q[2]*cpnt.rb*cpnt.za+
					crect.q[3]*cpnt.ra*cpnt.za)/(crect.hr*crect.hz);
	}
	return RETCODE_OK;
}

struct gen
{
	int type;
	virtual void read(ifstream &inf){};
};

struct gel : public gen
{
	PointXYZ A,B;
	gel(){type=1;}
	void read(ifstream &inf){inf>>A.x>>A.y>>A.z>>B.x>>B.y>>B.z;}
};

struct vel : public gen
{
	PointXYZ A,B;
	vel(){type=2;}
	void read(ifstream &inf){inf>>A.x>>A.y>>A.z>>B.x>>B.y>>B.z;}
};

struct ced : public gen
{
	PointXYZ A;
	double r;
	ced(){type=3;}
	void read(ifstream &inf){inf>>A.x>>A.y>>A.z>>r;}
};

struct ped : public gen
{
	PointXYZ A;
	int n;
	vector<PointXYZ> P;
	ped(){type=4;}
	void read(ifstream &inf)
	{
		int i;
		inf>>A.x>>A.y>>A.z>>n;
		P.resize(n);
		for(i=0;i<n;i++)
		{
			PointXYZ &B=P[i];
			inf>>B.x>>B.y>>B.z;
		}
	}
};

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

void SummingFiles(vector<int> &vNrec,vector<int> &GenToPls,int ipls,int n,char *dv,char *TmpStr,vector<gen *> &gens)
{
	int i,j;
	ifstream inf;
	ofstream ofp;
	vector<double> res;
	double tmpd;
	int nrec;

	nrec=vNrec[ipls];

	res.resize(nrec);
	for(j=0;j<nrec;j++){res[j]=0.0;}

	for(i=0;i<n;i++)
	{
		if(GenToPls[i]==ipls)
		{
			sprintf(TmpStr,"%s_2d.%d",dv,i+1);
			inf.open(TmpStr);
			for(j=0;j<nrec;j++)
			{
				inf>>tmpd;
				if(gens[i]->type==2)
				{
					res[j]-=tmpd;
				}
				else
				{
					res[j]+=tmpd;
				}
			}
			inf.close();
			inf.clear();
		}
	}

	sprintf(TmpStr,"%s_2d.000.%d",dv,ipls+1);
	ofp.open(TmpStr);
	for(j=0;j<nrec;j++){ofp<<res[j]<<'\n';}
	ofp.close();
	ofp.clear();
}

int main()
{
	int i,j,k,l,n,m,retc,ipls,npls,igloc;
	char PathTo2D[256],TmpStr[256];
	vector<double> v2;
	ifstream inf;
	FILE *fp;
	vector<gen *> gens;
	ofstream outfMN,outfMOON;
	vector<int> nGens,nGels,nVels,vNrec,vIGrec,vJGrec,GenToPls;
	int nGelsAll;

	logfile.open("LogOutputStat2D");

	inf.open("PathTo2D");
	if(inf)
	{
		inf>>PathTo2D;
		inf.close();
	}
	else
	{
		sprintf(PathTo2D,".");
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

	nGens.resize(npls);
	inf.open("srsclcsz");
	if(!inf)
	{
		logfile<<"Error in reading file srsclcsz"<<endl;
		return 1;
	}
	for(ipls=0;ipls<npls;ipls++){inf>>nGens[ipls];}
	inf.close();
	inf.clear();

	nGels.resize(npls);
	inf.open("srsclcgsz");
	if(!inf)
	{
		logfile<<"Error in open file "<<"srsclcgsz"<<endl;
		cout<<"Error in open file "<<"srsclcgsz"<<endl;
		return 1;
	}
	for(i=0;i<npls;i++){inf>>nGels[i];}
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
	for(i=0;i<npls;i++){inf>>nVels[i];}
	inf.close();
	inf.clear();

	inf.open("gens");
	if(!inf)
	{
		logfile<<"Error in reading file gens"<<endl;
		return 1;
	}
	inf>>n;
	gens.resize(n);
	for(i=0;i<n;i++)
	{
		inf>>retc;
		if(retc==1){gens[i]=(gen *)(new gel());}
		else if(retc==2){gens[i]=(gen *)(new vel());}
		else if(retc==3){gens[i]=(gen *)(new ced());}
		else if(retc==4){gens[i]=(gen *)(new ped());}
		else 
		{
			cout<<"Unknown generator type= "<<retc<<endl;
			logfile<<"Unknown generator type= "<<retc<<endl;
			return 1;
		}
		gens[i]->type=retc;
		gens[i]->read(inf);
	}
	inf.close();
	inf.clear();

	GenToPls.resize(n);
	vIGrec.resize(n+1);
	vJGrec.clear();
	vIGrec[0]=0;
	for(i=0;i<n;i++)
	{
		if(i<nGelsAll)
		{
			k=0;
			m=0;
			for(ipls=0;ipls<npls;ipls++)
			{
				if(i<k+nGels[ipls])
				{
					GenToPls[i]=ipls;
					vIGrec[i+1]=vIGrec[i]+vNrec[ipls];
					for(l=0;l<vNrec[ipls];l++)
					{
						vJGrec.push_back(m+l);
					}
					break;
				}
				k+=nGels[ipls];
				m+=vNrec[ipls];
			}
		}
		else
		{
			k=nGelsAll;
			m=0;
			for(ipls=0;ipls<npls;ipls++)
			{
				if(i<k+nVels[ipls])
				{
					GenToPls[i]=ipls;
					vIGrec[i+1]=vIGrec[i]+vNrec[ipls];
					for(l=0;l<vNrec[ipls];l++)
					{
						vJGrec.push_back(m+l);
					}
					break;
				}
				k+=nVels[ipls];
				m+=vNrec[ipls];
			}
		}
	}


	MeshForVP meshvp(moUnloadSolutions);

	if ((retc=meshvp.Read(PathTo2D))!=0)
		return retc;

	inf.open("xyzmn");
	if(!inf)
	{
		logfile<<"Error in reading file xyzmn"<<endl;
		return 1;
	}
	inf>>meshvp.qMN;
	inf.close();
	inf.clear();

	v2.resize(meshvp.nc*n);

	sprintf(TmpStr,"%s\\v2.dat",PathTo2D);
	fp=fopen(TmpStr,"rb");
	fread(&(v2.front()),sizeof(double),meshvp.nc*n,fp);
	fclose(fp);

	igloc=0;
	ipls=0;
	for(k=0;k<n;k++)
	{
		if(gens[k]->type==1)
		{
			meshvp.A=((gel *)(gens[k]))->A;
			meshvp.B=((gel *)(gens[k]))->B;
			meshvp.meshoptions=moUnloadSolutions;
		}
		else if(gens[k]->type==2)
		{
			meshvp.A=((vel *)(gens[k]))->A;
			meshvp.B=((vel *)(gens[k]))->B;
			meshvp.meshoptions=moUnloadSolutions|moForVEL;
		}
		else if(gens[k]->type==3)
		{
			meshvp.A=((ced *)(gens[k]))->A;
			meshvp.cedr=((ced *)(gens[k]))->r;
			meshvp.meshoptions=moUnloadSolutions|moForCED;
		}
		else if(gens[k]->type==4)
		{
			meshvp.A=((ped *)(gens[k]))->A;
			meshvp.NvB=((ped *)(gens[k]))->n;
			meshvp.vB.resize(meshvp.NvB);
			for(i=0;i<meshvp.NvB;i++){meshvp.vB[i]=((ped *)(gens[k]))->P[i];}
			meshvp.meshoptions=moUnloadSolutions|moForPED;
		}
		else 
		{
			cout<<"Unknown generator type= "<<gens[k]->type<<endl;
			logfile<<"Unknown generator type= "<<gens[k]->type<<endl;
			return 1;
		}

		if((retc=meshvp.ReadRecievers(".", 0, NULL))!=0)
			return retc;

		if ((retc=meshvp.GetWeights(&(v2[meshvp.nc*k])))!=0)
			return retc;
			
		if ((retc=ProcessResCalcPoints(&meshvp))!=0)
			return retc;

		sprintf(TmpStr,"dV_2d.%d",k+1);
		outfMN.open(TmpStr);
		sprintf(TmpStr,"ddV_2d.%d",k+1);
		outfMOON.open(TmpStr);
		for(l=vIGrec[k];l<vIGrec[k+1];l++)
		{
			double vM_a,vN_a,vO_a,vM_b,vN_b,vO_b;

			i=vJGrec[l];

			vM_a=meshvp.ResCalcPointArray[meshvp.boundsVforMNfromA.beginP+3*i  ].res[0].x;
			vN_a=meshvp.ResCalcPointArray[meshvp.boundsVforMNfromA.beginP+3*i+1].res[0].x;
			vO_a=meshvp.ResCalcPointArray[meshvp.boundsVforMNfromA.beginP+3*i+2].res[0].x;
			
			if(gens[k]->type==4)
			{
				double coeff;
				vM_b=vN_b=vO_b=0.0;
				if(meshvp.NvB)
				{
					coeff=1.0/meshvp.NvB;
					for(j=0;j<meshvp.NvB;j++)
					{
						vM_b+=meshvp.ResCalcPointArray[meshvp.boundsVforMNfromB.beginP+3*(i+meshvp.qMN*j) ].res[0].x;
						vN_b+=meshvp.ResCalcPointArray[meshvp.boundsVforMNfromB.beginP+3*(i+meshvp.qMN*j)+1].res[0].x;
						vO_b+=meshvp.ResCalcPointArray[meshvp.boundsVforMNfromB.beginP+3*(i+meshvp.qMN*j)+2].res[0].x;
					}
					vM_b*=coeff;
					vN_b*=coeff;
					vO_b*=coeff;
				}
			}
			else
			{
				vM_b=meshvp.ResCalcPointArray[meshvp.boundsVforMNfromB.beginP+3*i  ].res[0].x;
				vN_b=meshvp.ResCalcPointArray[meshvp.boundsVforMNfromB.beginP+3*i+1].res[0].x;
				vO_b=meshvp.ResCalcPointArray[meshvp.boundsVforMNfromB.beginP+3*i+2].res[0].x;
			}
			outfMN<<1000*(vM_a+vM_b-vN_a-vN_b)<<endl;
			outfMOON<<1000*((vM_a+vM_b-vO_a-vO_b)-(vO_a+vO_b-vN_a-vN_b))<<endl;
		}
		outfMN.close();
		outfMOON.close();

		igloc++;
		if(igloc==nGens[ipls])
		{
			igloc=0;
			ipls++;
		}
	}

	for(ipls=0;ipls<npls;ipls++)
	{
		SummingFiles(vNrec,GenToPls,ipls,n,"dV",TmpStr,gens);
		SummingFiles(vNrec,GenToPls,ipls,n,"ddV",TmpStr,gens);
	}

	logfile.close();
	logfile.clear();

	return 0;
}
