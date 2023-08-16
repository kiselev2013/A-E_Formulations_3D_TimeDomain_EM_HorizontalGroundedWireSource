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
 *  This file contains the code for extracting anomal edges for 3D VFEM task
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin                         
 *  Novosibirsk State Technical University,                 
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia       
 *  vdv_wk@mail.ru                                          
 *  Version 2.0 January 16, 2023                            
*/

#include "stdafx.h"   
#include "time_approx.h"
#include "in_out.h"
#include "vfem_const.h"

extern ofstream logfile;
extern void CopyToFile(const char* _from, const char* _to);
extern void Memory_allocation_error(const char *var, const char *func);
extern void Cannot_open_file(const char *fname, const char *func);
extern double Norm_Euclid(double *a, int n);
extern double Projection_On_Axis(double *v,double *o);
extern double Scal(double *a, double *b, int n);
extern void Mult_Plot(double *a, double *x, double *y, int n);

extern bool CheckStop(void);

bool IsFileExist(char *fname)
{
	bool flag;
	ifstream inf;

	flag=false;
	inf.open(fname);
	if(inf)
	{
		flag=true;
		inf.close();
	}
	inf.clear();

	return flag;
}

int inputvec(int p_kpar,int *p_n,int *p_nc,int (**ed)[25],int (**edges)[2])
{
	int i,j;
	FILE *fp;
	ifstream inf;
	const int size_i=sizeof(int);
	const int size_d=sizeof(double);

	int p_vc,mv;

	inf.open("tsize3d_.dat");
	if(!inf){
		printf("Error open file tsize3d.dat");
		return 1;
	}
	inf>>j;
	inf>>*p_n;
	inf.close();
	inf.clear();

	*p_nc=*p_n;

	if(!(fp=fopen("nodesforedges.dat","rb"))){
		printf("Error open file nodesforedges.dat");
		return 1;
	}

	if(!((*edges)=new int[*p_n][2]))return 1;
	for(i=0;i<*p_n;i++){
		fread((*edges)[i],size_i,2,fp);
		(*edges)[i][0]--;(*edges)[i][1]--;
	}
	fclose(fp);

	if(!(fp=fopen("edges.dat","rb"))){
		printf("Error open file edges.dat");
		return 1;
	}

	if(!((*ed)=new int[p_kpar][25]))return 1;
	for(i=0;i<p_kpar;i++){
		fread((*ed)[i],size_i,25,fp);
		for(j=0;j<25;j++)(*ed)[i][j]--;
	}
	fclose(fp);
	return 0;
}

Time_approx_for_vfem::Time_approx_for_vfem(bool _flag_vfem_for_polygon, const void* _NormBFromPolygon, bool flag_vfem_for_line, const void* _NormResFromLine)
{
	anomal_edges = NULL;
	d = NULL;
	nodes_f = NULL;
	is_field_in_node = NULL;
	Tmap = NULL;
}
Time_approx_for_vfem::~Time_approx_for_vfem()
{
	if(anomal_edges) {delete [] anomal_edges; anomal_edges=NULL;}
	if(d) {delete d; d=NULL;}
	if(nodes_f) {delete [] nodes_f; nodes_f=NULL;}
	if(is_field_in_node) {delete [] is_field_in_node; is_field_in_node=NULL;}
	if(Tmap) {delete Tmap; Tmap=NULL;}
}
int Time_approx_for_vfem::Read_data()
{
	int i,n_anom_mu;

	d = new Vec_Prep_Data();
	if(d == NULL)
		Memory_allocation_error("d", "Time_approx_for_vfem::Read_data");
	d->Read_mesh_for_nonstat_problem("");
	Unload_anomal_nodes("xyzVectorE0");

	n_anom_mu=0;
	for(i=0;i<d->n_materials;i++)
	{
		n_anom_mu+=(fabs(d->mu3d[i]-d->mu0[i])/MU_0>1e-6);
	}
	if(n_anom_mu)
	{
		Unload_anomal_nodesA0("xyzVectorA0");
	}

	return 0;
}
int Time_approx_for_vfem::Unload_anomal_nodes(char *fname)
{
	FILE *fp=NULL;
	int i, j, k, m, nn[2], li, lj;
	const double eps = 1e-6;
	double px,py,pz,len;
	bool fcont;

	Tmap = new T_Mapping_Vec(d->nver, d->xyz, d->kuzlov, d->kpar);
	if(Tmap == 0)
		Memory_allocation_error("Tmap", "Time_approx_for_vfem::Two_Layers_Schema");

	if(inputvec(Tmap->kpar,&(Tmap->n),&(Tmap->n_c),&(Tmap->ed),&(Tmap->edges)))return 1;
	n = Tmap->n_c;

	isEdgesAnomal.resize(Tmap->n);
	EdgesCrd[0].resize(Tmap->n);
	EdgesCrd[1].resize(Tmap->n);
	EdgesCrd[2].resize(Tmap->n);

	for(i=0;i<Tmap->n;i++){isEdgesAnomal[i]=false;}

	elemRenum.resize(d->kpar, -1);
	isElemAnomal.resize(d->kpar, false);

	k=0;
	for(i=0;i<d->kpar;i++)
	{
		if (fabs(d->sigma0[d->nvkat[i]]-d->sigma3d[d->nvkat[i]])>eps)
		{
			elemRenum[i] = k;
			isElemAnomal[i] = true;

			for(j=0;j<12;j++)
			{
				isEdgesAnomal[Tmap->ed[i][j]]=true;
			}

			k++;
		}
	}

	double z0;
	z0=zlay[nzl-1]+1e-3;
	for(i=0;i<Tmap->n;i++)
	{
		if(isEdgesAnomal[i])
		{
			nn[0]=Tmap->edges[i][0];
			nn[1]=Tmap->edges[i][1];
			if(d->xyz[nn[0]][2]>z0 || d->xyz[nn[1]][2]>z0)
			{
				isEdgesAnomal[i]=false;
			}
		}
	}

	n_nodes_f=0;
	for(i=0;i<Tmap->n;i++)
	{
		if(isEdgesAnomal[i])
		{
			nn[0]=Tmap->edges[i][0];
			nn[1]=Tmap->edges[i][1];
			EdgesCrd[0][i]=0.5*(d->xyz[nn[0]][0]+d->xyz[nn[1]][0]);
			EdgesCrd[1][i]=0.5*(d->xyz[nn[0]][1]+d->xyz[nn[1]][1]);
			EdgesCrd[2][i]=0.5*(d->xyz[nn[0]][2]+d->xyz[nn[1]][2]);
			n_nodes_f++;
		}
	}

	for(i=0;i<d->kpar;i++)
	{
		for(j=0;j<12;j++)
		{
			if(isEdgesAnomal[Tmap->ed[i][j]])
			{
				isElemAnomal[i]=true;
				break;
			}
		}
	}

	ofstream fout;

	fout.open("isElemAnomal");
	for(i=0;i<d->kpar;i++)
	{
		fout<<isElemAnomal[i]<<'\n';
	}
	fout.close();
	fout.clear();

	fout.open("mtrsizeforE0");
	for(i=0;i<Tmap->n;i++)
	{
		fout<<isEdgesAnomal[i]<<'\n';
	}
	fout.close();
	fout.clear();

	fout.open("xyzVectorE0");
	fout<<scientific<<setprecision(14);
	fout<<n_nodes_f<<'\n';
	for(i=0;i<Tmap->n;i++)
	{
		if(isEdgesAnomal[i])
		{
			fout<<EdgesCrd[0][i]<<' '<<EdgesCrd[1][i]<<' '<<EdgesCrd[2][i]<<'\n';
		}
	}
	fout.close();
	fout.clear();

	fout.open("TgCompE0");
	for(i=0;i<Tmap->n;i++)
	{
		if(isEdgesAnomal[i])
		{
			nn[0]=Tmap->edges[i][0];
			nn[1]=Tmap->edges[i][1];
			px=d->xyz[nn[1]][0]-d->xyz[nn[0]][0];
			py=d->xyz[nn[1]][1]-d->xyz[nn[0]][1];
			pz=d->xyz[nn[1]][2]-d->xyz[nn[0]][2];
			len=sqrt(px*px+py*py+pz*pz);
			fout<<px/len<<' '<<py/len<<' '<<pz/len<<'\n';
		}
	}
	fout.close();
	fout.clear();

	return 0;
}

int Time_approx_for_vfem::Unload_anomal_nodesA0(char *fname)
{
	FILE *fp=NULL;
	int i, j, k, m, nn[2];
	const double eps = 1e-6;
	double px,py,pz,len;
	bool BoundAnomal;

	for(i=0;i<Tmap->n;i++){isEdgesAnomal[i]=false;}
	for(i=0;i<d->kpar;i++){isElemAnomal[i]=false;}

	k=0;
	for(i=0;i<d->kpar;i++)
	{
		if (fabs(d->mu0[d->nvkat[i]]-d->mu3d[d->nvkat[i]])>eps)
		{
			isElemAnomal[i] = true;

			for(j=0;j<12;j++)
			{
				isEdgesAnomal[Tmap->ed[i][j]]=true;
			}

			k++;
		}
	}

	n_nodes_f=0;
	for(i=0;i<Tmap->n;i++)
	{
		if(isEdgesAnomal[i])
		{
			nn[0]=Tmap->edges[i][0];
			nn[1]=Tmap->edges[i][1];
			EdgesCrd[0][i]=0.5*(d->xyz[nn[0]][0]+d->xyz[nn[1]][0]);
			EdgesCrd[1][i]=0.5*(d->xyz[nn[0]][1]+d->xyz[nn[1]][1]);
			EdgesCrd[2][i]=0.5*(d->xyz[nn[0]][2]+d->xyz[nn[1]][2]);
			n_nodes_f++;
		}
	}

	ofstream fout;

	fout.open("isElemAnomal_A");
	for(i=0;i<d->kpar;i++)
	{
		fout<<isElemAnomal[i]<<'\n';
	}
	fout.close();
	fout.clear();

	fout.open("mtrsizeforA0");
	for(i=0;i<Tmap->n;i++)
	{
		fout<<isEdgesAnomal[i]<<'\n';
	}
	fout.close();
	fout.clear();

	fout.open("xyzVectorA0");
	fout<<n_nodes_f<<'\n';
	for(i=0;i<Tmap->n;i++)
	{
		if(isEdgesAnomal[i])
		{
			fout<<EdgesCrd[0][i]<<' '<<EdgesCrd[1][i]<<' '<<EdgesCrd[2][i]<<'\n';
		}
	}
	fout.close();
	fout.clear();

	fout.open("TgCompA0");
	for(i=0;i<Tmap->n;i++)
	{
		if(isEdgesAnomal[i])
		{
			nn[0]=Tmap->edges[i][0];
			nn[1]=Tmap->edges[i][1];
			px=d->xyz[nn[1]][0]-d->xyz[nn[0]][0];
			py=d->xyz[nn[1]][1]-d->xyz[nn[0]][1];
			pz=d->xyz[nn[1]][2]-d->xyz[nn[0]][2];
			len=sqrt(px*px+py*py+pz*pz);
			fout<<px/len<<' '<<py/len<<' '<<pz/len<<'\n';
		}
	}
	fout.close();
	fout.clear();

	return 0;
}
