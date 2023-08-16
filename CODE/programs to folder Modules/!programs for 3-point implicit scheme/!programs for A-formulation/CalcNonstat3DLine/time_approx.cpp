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
 *  This file contains the code for calculating nonstationary 3D VFEM task
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
#include "OutputArbitrary.h"
#include "in_out.h"
#include "vfem_const.h"

#include"TimeWatcher.h"


extern void CloseProgramm(int rcode);

extern void  mult_symmetr(int *ig, int *jg, double *gg, double *di, double* x, double *y, int n);

extern ofstream logfile;

extern int output_v3_vec(double *R, int L, int dimrec);
extern void CopyToFile(const char* _from, const char* _to);

extern void Memory_allocation_error(const char *var, const char *func);
extern void Cannot_open_file(const char *fname, const char *func);

extern double Norm_Euclid(double *a, int n);
extern double Projection_On_Axis(double *v,double *o);
extern double Scal(double *a, double *b, int n);
extern void Mult_Plot(double *a, double *x, double *y, int n);

extern int FindInIntForElem(double *vec,double elem,int size);
extern int FindInIntMass(int *vec,int elem,int size);

extern void FindLocalCoordinates(pv::Point3D &R,pv::Point3D *HexPnt,double *lc);
extern void GetNodeSolutionByLocalcCoords(double *lc,double *q,double &Val);

extern bool CheckStop(void);

extern void CloseProgramm(int code);

extern double GlobalStatus;

Time_approx_for_vfem::Time_approx_for_vfem(int _npls)
{
	eps = 1e-6;
	n_max_iter_per_layer = 10000;

	npls=_npls;



	forLine = true;


	give_out = NULL;
	anomal_edges = NULL;

	u_j = NULL;
	u_j1 = NULL;
	u_j2 = NULL;
	cuj = NULL;
	cujsum = NULL;
	Slae = NULL;
	b_gg_all = NULL;
	c_gg_all = NULL;
	c_gg_all_dsig = NULL;
	b_di_all = NULL;
	c_di_all = NULL;
	c_di_all_dsig = NULL;
	d = NULL;
	nodes_f = NULL;
	En_nodes = NULL;
	is_field_in_node = NULL;
	En_edges = NULL;
	Pt = NULL;
	Tmap = NULL;
	Bc = NULL;
	help = NULL;
	help1 = NULL;
	help2 = NULL;
	help3 = NULL;
	help4 = NULL;
	help5 = NULL;
	help6 = NULL;
	help7 = NULL;
	help8 = NULL;	

	An_edges = NULL;
	d2An_edges = NULL;

	Buff=NULL;

	AnomalType=0;
	
	nthreads=1;

	eps_reg=1e-12;

	npntE0=0;

	ia=-1;

	ForInv=0;

	forCED = false;

	iVP=false;

	y_omp = NULL;

}
Time_approx_for_vfem::~Time_approx_for_vfem()
{
	if(give_out) {delete give_out; give_out=NULL;}
	if(anomal_edges) {delete [] anomal_edges; anomal_edges=NULL;}
	if(u_j) {delete [] u_j; u_j=NULL;}
	if(u_j1) {delete [] u_j1; u_j1=NULL;}
	if(u_j2) {delete [] u_j2; u_j2=NULL;}
	if(cuj) {delete [] cuj; cuj=NULL;}
	if(cujsum) {delete [] cujsum; cujsum=NULL;}
	if(Slae) {delete Slae; Slae=NULL;}
	if(b_gg_all) {delete [] b_gg_all; b_gg_all=NULL;}
	if(c_gg_all) {delete [] c_gg_all; c_gg_all=NULL;}
	if(c_gg_all_dsig) {delete [] c_gg_all_dsig; c_gg_all_dsig=NULL;}
	if(b_di_all) {delete [] b_di_all; b_di_all=NULL;}
	if(c_di_all) {delete [] c_di_all; c_di_all=NULL;}
	if(c_di_all_dsig) {delete [] c_di_all_dsig; c_di_all_dsig=NULL;}
	if(d) {delete d; d=NULL;}
	if(nodes_f) {delete [] nodes_f; nodes_f=NULL;}
	if(En_nodes) {delete [] En_nodes; En_nodes=NULL;}
	if(is_field_in_node) {delete [] is_field_in_node; is_field_in_node=NULL;}
	if(En_edges) {delete [] En_edges; En_edges=NULL;}
	if(Pt) {delete Pt; Pt=NULL;}
	if(Tmap) {delete Tmap; Tmap=NULL;}
	if(Bc) {delete Bc; Bc=NULL;}
	if(help) {delete [] help; help=NULL;}
	if(help1) {delete [] help1; help1=NULL;}
	if(help2) {delete [] help2; help2=NULL;}
	if(help3) {delete [] help3; help3=NULL;}
	if(help4) {delete [] help4; help4=NULL;}
	if(help5) {delete [] help5; help5=NULL;}
	if(help6) {delete [] help6; help6=NULL;}
	if(help7) {delete [] help7; help7=NULL;}
	if(help8) {delete [] help8; help8=NULL;}
	if(d2An_edges) {delete [] d2An_edges; d2An_edges=NULL;}
	if(An_edges) {delete [] An_edges; An_edges=NULL;}

	if(y_omp) {delete [] y_omp; y_omp=NULL;}

	if(Buff) {delete [] Buff;Buff=NULL;}

}
int Time_approx_for_vfem::Read_data(vector<int> &_RecvPlsIgB,vector<int> &_RecvPlsIgE)
{
	d = new Vec_Prep_Data();
	if(d == NULL)
		Memory_allocation_error("d", "Time_approx_for_vfem::Read_data");

		d->Read_mesh_for_nonstat_problem("0");


	d->Read_infite0();
	ntime = d->ntime;
	time = d->time;
	finSP = time[ntime-1];

	RecvPlsIgB=_RecvPlsIgB;
	RecvPlsIgE=_RecvPlsIgE;

	return 0;
}

int Time_approx_for_vfem::LoadVectorE0()
{
	int i,j;
	ifstream inf;

	inf.open("xyzVectorE0");
	if(!inf)
	{
		cout<<"Can't open file "<<"xyzVectorE0"<<'\n';
		logfile<<"Can't open file "<<"xyzVectorE0"<<'\n';
		return 1;
	}
	inf>>n_nodes_f;
	inf.close();
	inf.clear();

	nodes_f = new int[n_nodes_f];
	if(nodes_f == 0)
		return 1;

	is_field_in_node=new bool[d->kuzlov];
	if(is_field_in_node == 0)
		return 1;

	inf.open("xyzVectorE0n");
	if(!inf)
	{
		cout<<"Can't open file "<<"xyzVectorE0n"<<'\n';
		logfile<<"Can't open file "<<"xyzVectorE0n"<<'\n';
		return 1;
	}
	for(i=0;i<n_nodes_f;i++)
	{
		inf>>j;
		is_field_in_node[j-1]=true;
	}
	inf.close();
	inf.clear();
	
	j = 0;
	for(i=0; i<d->kuzlov; i++)
		if (is_field_in_node[i])
		{
			nodes_f[j] = i;
			j++;
		}

	En_nodes = new double[d->kuzlov][3];
	if(En_nodes == 0)
		return 1;

	return 0;
}

int inputvec(int p_kpar,int *p_n,int *p_nc,int (**ed)[25],int (**edges)[2])
{
	int i,j;
	FILE *fp;
	ifstream inf;

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

int Time_approx_for_vfem::Prepare_To_Schema()
{
	int i,j,m;
	const char v3[] = {"v3."};
	In_Out r;
	ifstream inf;
	int tmpi;
	double tmpd;

	if(Tmap == NULL)
	{
		Tmap = new T_Mapping_Vec(d->nver, d->xyz, d->kuzlov, d->kpar);
		if(Tmap == 0)
			Memory_allocation_error("Tmap", "Time_approx_for_vfem::Two_Layers_Schema");


		if(inputvec(Tmap->kpar,&(Tmap->n),&(Tmap->n_c),&(Tmap->ed),&(Tmap->edges)))return 1;
		n = Tmap->n_c;
	}

	if(u_j == NULL)
	{
		u_j = new double[Tmap->n*npls];
		if(u_j == NULL)
			Memory_allocation_error("u_j", "Time_approx_for_vfem::Two_Layers_Schema");
	}

	if(u_j1 == NULL)
	{
		u_j1 = new double[Tmap->n*npls];
		if(u_j1 == NULL)
			Memory_allocation_error("u_j1", "Time_approx_for_vfem::Two_Layers_Schema");
	}

	if(u_j2 == NULL)
	{
		u_j2 = new double[Tmap->n*npls];
		if(u_j2 == NULL)
			Memory_allocation_error("u_j2", "Time_approx_for_vfem::Two_Layers_Schema");
	}

	if(cuj == NULL)
	{
		cuj = new double[Tmap->n];
		if(cuj == NULL)
			Memory_allocation_error("cuj", "Time_approx_for_vfem::Two_Layers_Schema");
	}

	if(cujsum == NULL)
	{
		cujsum = new double[Tmap->n];
		if(cujsum == NULL)
			Memory_allocation_error("cuj", "Time_approx_for_vfem::Two_Layers_Schema");
	}

	if(y_omp == NULL)
	{
		if((y_omp = new double[n]) == NULL)
			Memory_allocation_error("y_omp", "Time_approx_for_vfem::Two_Layers_Schema");
	}

	
	isEdgesAnomal.resize(Tmap->n);
	for(i=0;i<Tmap->n;i++){isEdgesAnomal[i]=false;}

	int k;
	ifstream inf0;
	npntE0=0;
	inf0.open("mtrsizeforE0");
	if(inf0)
	{
		for(i=0;i<Tmap->n;i++)
		{
			inf0>>k;
			isEdgesAnomal[i]=k;
				npntE0+=k;
		}
		inf0.close();
	}
	inf0.clear();

	if(npntE0 && En_edges == NULL)
	{
		En_edges = new double[Tmap->n*npls]; //     ( -,  -)
		if(En_edges == NULL)
			Memory_allocation_error("En_edges", "Time_approx_for_vfem::Two_Layers_Schema");

		for(i=0; i<Tmap->n*npls; i++)
			En_edges[i] = 0.0;
	}

	if(Pt == NULL)
	{
		Pt = new T_Portrait((int*)Tmap->ed, Tmap->n, Tmap->n_c, d->kpar);
		if(Pt == 0)
			Memory_allocation_error("Pt", "Time_approx_for_vfem::Two_Layers_Schema");

		Pt->Gen_Portrait();
	}

	if(Slae == NULL)
	{
		Slae = new T_Global_SLAE_Vec(d, Tmap, Pt, npls);
		if(Slae == NULL)
			Memory_allocation_error("Slae", "Time_approx_for_vfem::Two_Layers_Schema");
	}

	if (b_gg_all == NULL)
	{
		b_gg_all = new double[Slae->ig_n_1];
		if(b_gg_all == 0)
			Memory_allocation_error("b_gg_all", "Time_approx_for_vfem::Two_Layers_Schema");
	}

	if (c_gg_all == NULL)
	{
		c_gg_all = new double[Slae->ig_n_1];
		if(c_gg_all == 0)
			Memory_allocation_error("c_gg_all", "Time_approx_for_vfem::Two_Layers_Schema");
	}

	if (c_di_all == NULL)
	{
		c_di_all = new double[Slae->n];
		if(c_di_all == 0)
			Memory_allocation_error("c_di_all", "Time_approx_for_vfem::Two_Layers_Schema");
	}

	if (b_di_all == NULL)
	{
		b_di_all = new double[Slae->n];
		if(b_di_all == 0)
			Memory_allocation_error("b_di_all", "Time_approx_for_vfem::Two_Layers_Schema");
	}

	if(npntE0)
	{
		if (c_gg_all_dsig == NULL)
		{
			c_gg_all_dsig = new double[Slae->ig_n_1];
			if(c_gg_all_dsig == 0)
				Memory_allocation_error("c_gg_all_dsig", "Time_approx_for_vfem::Two_Layers_Schema");
		}

	if (c_di_all_dsig == NULL)
	{
		c_di_all_dsig = new double[Slae->n];
		if(c_di_all_dsig == 0)
				Memory_allocation_error("c_di_all_dsig", "Time_approx_for_vfem::Two_Layers_Schema");
		}
	}

	Slae->Assembling(1, En_edges);
	for (i=0; i<Slae->ig_n_1; i++)
	{
		c_gg_all[i] = Slae->ggl[i];
		Slae->ggl[i] = 0.0;
	}
	for (i=0; i<Slae->n; i++)
	{
		c_di_all[i] = Slae->di[i];
		Slae->di[i] = 0.0;
	}

	if(npntE0)
	{
	Slae->Assembling(4, En_edges);
	for (i=0; i<Slae->ig_n_1; i++)
	{
		c_gg_all_dsig[i] = Slae->ggl[i];
		Slae->ggl[i] = 0.0;
	}
	for (i=0; i<Slae->n; i++)
	{
		c_di_all_dsig[i] = Slae->di[i];
		Slae->di[i] = 0.0;
	}
	}

	Slae->Assembling(2, En_edges);
	for (i=0; i<Slae->ig_n_1; i++){b_gg_all[i] = Slae->ggl[i];}
	for (i=0; i<Slae->n; i++){b_di_all[i] = Slae->di[i];}


	if(Bc==NULL) if((Bc = new T_Bound_Cond(d->kuzlov, Tmap->n_c, d->kt1, d->l13d, d->xyz, Tmap->edges, d->kpar, Tmap->ed, d->nver))==NULL)
		Memory_allocation_error("Bc", "Time_approx_for_vfem::Two_Layers_Schema");

	for(i=0; i<npls*Tmap->n; i++){u_j2[i]=u_j1[i]=u_j[i]=0.0;}

	if(give_out==NULL)
		give_out = new Give_out_vec_loop(d, Tmap, RecvPlsIgB, RecvPlsIgE, npls);
	if(give_out==NULL)
		Memory_allocation_error("give_out", "Time_approx_for_vfem::Two_Layers_Schema");

	Subdomain::npls=npls;
	Subdomain::ntimes=max_dec_size;

	double (*xyz0)[3];
	xyz0=NULL;
	d->xyzt=NULL;
	
	inf.open("xyz0.dat",ios::binary);
	if(inf)
	{
		xyz0=new double[d->kuzlov][3];
		for(i=0;i<d->kuzlov;i++)
		{
			inf>xyz0[i][0]>xyz0[i][1]>xyz0[i][2];
		}
		inf.close();
		d->xyzt=d->xyz;
		d->xyz=xyz0;
	}
	else
	{
		d->xyzt=d->xyz;
	}
	inf.clear();

	if(d->n_pointresB)
	{
		give_out->resultantB = new OutputResultant3d(this,vtWithoutDiscontinuity,nthreads);
		give_out->resultantB->Prepare(3);
		give_out->vvtb.resize(3);
		give_out->vvtb[0]=vtRotxA;
		give_out->vvtb[1]=vtRotyA;
		give_out->vvtb[2]=vtRotzA;
	}

	if(d->n_pointresE)
	{
		give_out->resultantA = new OutputResultant3d(this,vtWithDiscontinuity,nthreads);
		give_out->resultantA->Prepare(3);
		give_out->vvta.resize(3);
		give_out->vvta[0]=vtAx;
		give_out->vvta[1]=vtAy;
		give_out->vvta[2]=vtAz;
	}

	if(xyz0)
	{
		d->xyz=d->xyzt;
		delete [] xyz0;
		xyz0=NULL;
		d->xyzt=NULL;
	}

	give_out->vPre1.resize(ntime);
	give_out->vPre2.resize(ntime);

	give_out->vPre1[0]=-1;
	give_out->vPre2[0]=-1;

	inf.open("eps_reg");
	if(inf)
	{
		inf>>eps_reg;
		inf.close();
	}
	inf.clear();

	if(d->fdirect)
	{
		tmpi=d->LoadAndInitCurrent(npls);
		if(tmpi){return tmpi;}
	}

	return 0;
}

int ReadField3d(double *u,int n,int ind)
{
	FILE *fp;
	char buff[256];
	if(ind>=0)
	{
		sprintf(buff,"v3.%d",ind);
		cout<<"Reading file "<<buff<<endl;
		fp=fopen(buff,"rb");
		if(!fp){return 1;}
		fread(u,sizeof(double),n,fp);
		fclose(fp);
		fflush(fp);
	}
	return 0;
}

int Time_approx_for_vfem::Two_Layers_Schema(int nlst)
{
	int i,j,l,m,ret;
	int itime;
	char fname[30];
	char buf[20];
	const char v3[] = {"v3."};
	In_Out r;
	int ipls;
	double value,binn[3],bout[3];
	T_Brick brk;
	
	ifstream inf;
	ofstream ofp;

	if(Bc==NULL) if((Bc = new T_Bound_Cond(d->kuzlov, Tmap->n_c, d->kt1, d->l13d, d->xyz, Tmap->edges, d->kpar, Tmap->ed, d->nver))==NULL)
		Memory_allocation_error("Bc", "Time_approx_for_vfem::Two_Layers_Schema");

	bool a0Exist = false;

	inf.open("a0.edge",ios::binary);
	if(inf)
	{
		cout<<"Reading file "<<"a0.edge"<<endl;
		logfile<<"Reading open file "<<"a0.edge"<<endl;
		a0Exist=true;
		inf.read((char *)u_j,sizeof(double)*Slae->n*npls);
		inf.close();
		for(i=0;i<Tmap->n*npls;i++){u_j2[i]=u_j1[i]=u_j[i];}
	}
	inf.clear();

	if(!a0Exist)
	{
		double *a0Node=NULL;

		if((a0Node=new double[d->kuzlov*3*npls])==NULL) Memory_allocation_error("a0Node", "Time_approx_for_vfem::Two_Layers_Schema");
		

		inf.open("a0.dat",ios::binary);
		if(!inf)
		{
			cout<<"Can't open file "<<"a0.dat"<<endl;
			logfile<<"Can't open file "<<"a0.dat"<<endl;
			exit(1);
		}
		inf.read((char *)a0Node,sizeof(double)*d->kuzlov*3*npls);
		inf.close();
		inf.clear();


		a0Exist = true;

		if(a0Exist)
		{
			int v1,v2;
			double len;
			double axis[3];
			double mid_value[3];
			for(ipls=0;ipls<npls;ipls++)
			{
				for(i=0;i<n;i++)
				{
					v1=Tmap->edges[i][0];
					v2=Tmap->edges[i][1];
					for(j=0;j<3;j++)
					{
						axis[j]=d->xyz[v2][j]-d->xyz[v1][j];
						mid_value[j]=0.5*(a0Node[v1*3+j+ipls*d->kuzlov*3]+a0Node[v2*3+j+ipls*d->kuzlov*3]);
					}
					len=Norm_Euclid(axis,3);
					u_j2[i+ipls*Tmap->n]=u_j1[i+ipls*Tmap->n]=u_j[i+ipls*Tmap->n]=Projection_On_Axis(mid_value,axis)*len*0.5;
				}
			}
		}
		else
		{
			for(i=0;i<Tmap->n*npls;i++){u_j2[i]=u_j1[i]=u_j[i]=0.0;}
		}
		
		if(a0Node){delete [] a0Node; a0Node=NULL;}
	}



	double *h=new double[npls*Tmap->n];

	{
		for(itime=1;itime<3;itime++)
		{
			if(itime==1)
			{
				dt=time[itime]-time[itime-1]; 
			}

			if(npntE0 && !d->fdirect)
			{
				LoadE0FromRZ(itime,En_edges);
			}

			for(ipls=0;ipls<npls;ipls++)
			{
				ipls_cur=ipls;

				if(npntE0 && !d->fdirect)
				{
					mult_symmetr(Pt->ig, Pt->jg, c_gg_all_dsig, c_di_all_dsig, En_edges+(Tmap->n*ipls), Slae->pr+(Slae->n*ipls), n);
				}
				else
				{
					for(i=0;i<n;i++){Slae->pr[Slae->n*ipls+i] = 0.0;}
				}

				for(i=0;i<n;i++){cujsum[i]=u_j[Tmap->n*ipls+i]/dt;}
				mult_symmetr(Pt->ig, Pt->jg, c_gg_all, c_di_all, cujsum, cuj, n);
				for(i=0;i<n;i++){Slae->pr[Slae->n*ipls+i] += cuj[i];}
			}


			if(itime==1)
			{
				for(i=0; i<n; i++)
				{
					Slae->di[i] = c_di_all[i]/dt + b_di_all[i];
					for(j=Pt->ig[i]; j<=Pt->ig[i+1]-1; j++)
						Slae->ggl[j] = c_gg_all[j]/dt + b_gg_all[j];
				}

				for(i=0;i<Slae->n;i++){Slae->di[i]+=Slae->di[i]*1e-12;}

				Bc->Set_Dirichlet_Cond(Pt->ig, Pt->jg, Slae->di, Slae->ggl, Slae->pr, time[itime],0);
			}

			if(d->fdirect && d->nSrsRazb)
			{
				if(d->currentfunction[itime])
				{
					for(ipls=0;ipls<npls;ipls++)
					{
						m=(int)d->srs[ipls].size();
						for(l=0;l<m;l++)
						{
							loc_source &ls=d->srs[ipls][l];
							
							binn[0]=2.0*ls.ksi-1.0;
							binn[1]=2.0*ls.eta-1.0;
							binn[2]=2.0*ls.dzeta-1.0;

							for(j=0; j<12; j++)
							{
								int ii,i_mu,it,nn0,nn1;
								double edlen;
								
								ii=Tmap->ed[ls.elem][j];

								nn0=Tmap->edges[ii][0];
								nn1=Tmap->edges[ii][1];

								edlen=sqrt((d->xyz[nn1][0]-d->xyz[nn0][0])*(d->xyz[nn1][0]-d->xyz[nn0][0])+
									(d->xyz[nn1][1]-d->xyz[nn0][1])*(d->xyz[nn1][1]-d->xyz[nn0][1])+
									(d->xyz[nn1][2]-d->xyz[nn0][2])*(d->xyz[nn1][2]-d->xyz[nn0][2]));

								brk.Basis_func_on_reference_vec_par(j,binn,bout);
								value=(bout[0]*ls.dx+bout[1]*ls.dy+bout[2]*ls.dz)*d->currentfunction[itime]*2.0/edlen;
													
								Slae->pr[Slae->n*ipls+ii] += value;
							}
						}
					}
				}
			}

			for(ipls=0;ipls<npls;ipls++)
			{
				Bc->Set_Dirichlet_Cond_Pr(Slae->pr, time[itime],0,ipls);
			}

			if(itime==1)
			{
				prds.factorize(Slae->n,Pt->ig,Pt->jg,Slae->ggl,Slae->di,nthreads);
			}

			prds.solve_nrhs(npls,Slae->pr,u_j);


			if(itime>=nlst)
			{
				_itoa_s(itime, buf, 10);
				strcpy(fname, v3);
				strcat(fname, buf);

				ofp.open(fname,ios::binary);
				ofp.write((char *)u_j,npls*size_d*Slae->n);
				ofp.close();
				ofp.clear();
			}

			give_out->vPre1[itime]=itime-1;
			give_out->vPre2[itime]=itime-2;

			for(i=0; i<npls*Tmap->n; i++)
			{
				u_j2[i] = u_j1[i];
				u_j1[i] = u_j[i];
			}

			cout<<"time layer - "<<itime<<" done"<<endl;
		}
	}

	prds.stop_solver();

	return 0;
}
int Time_approx_for_vfem::Three_Layers_Schema(int tbeg,int tend,int i_u_1,int i_u_2,int nlst,int tbeg_pre,bool f_last_dec)
{
	int i,j,l,m,ret;
	int itime;
	char fname[30];
	char buf[256];
	const char v3[] = {"v3."};
	In_Out r;
	double mult, mult1, mult2;
	double t_j, t_j1, t_j2;  //        
	int tnum, tnum1, tnum2;

	__time64_t time_total, time_beg, time_end; //   
	__time64_t timestruct;
	FILE *f_time=NULL;
	int msg_iter;

	ifstream inf;
	ofstream ofp;
	bool fcalc;

	int ipls;
	double value,binn[3],bout[3];
	T_Brick brk;

	fcalc=false;

	cur_dec_strt=tbeg;
	give_out->cur_dec_strt_res=cur_dec_strt;
	Subdomain::ntimes=tend-tbeg;

	if(1>=tbeg && 1<tend)
	{
		sprintf(buf,"Start two-layer schema");
		TimeWatcher::getInstance().AddTimeSpot(clock(),buf);

		Two_Layers_Schema(nlst);

		//for(i=0; i<npls*Tmap->n; i++)
		//	u_j1[i] = u_j[i];

		//t_j1 = time[0];
		//t_j = time[1];

		//tnum1 = 0;
		//tnum = 1;

		t_j2 = time[0];
		t_j1 = time[1];
		t_j = time[2];

		tnum2 = 0;
		tnum1 = 1;
		tnum = 2;

		sprintf(buf,"Finish two-layer schema");
		TimeWatcher::getInstance().AddTimeSpot(clock(),buf);
	}
	else
	{
		for(i=0; i<npls*Tmap->n; i++){u_j2[i]=u_j1[i]=u_j[i]=0.0;}
	}

	if ((f_time=fopen("time_solver","w"))==0)
		Cannot_open_file("time_solver", "Time_approx_for_vfem::Three_Layers_Schema");

	int _ntime = ntime;

	bool FirstDecadeStep;

	FirstDecadeStep=true;

	double *h=new double[npls*Tmap->n];

	for(itime=3;itime<_ntime;itime++)//////////////////////////
	{
		give_out->fsdiff=false;

		if(itime>=tbeg && itime<tend)
		{
			if(CheckStop())break;

			if(!fcalc && itime!=3)
			{
				sprintf(buf,"Start preparing data for new decade");
				TimeWatcher::getInstance().AddTimeSpot(clock(),buf);

				if((ReadField3d(u_j1,npls*Slae->n,i_u_1))!=0)return 1;
				if((ReadField3d(u_j2,npls*Slae->n,i_u_2))!=0)return 1;
				give_out->fsdiff=true;
				sprintf(buf,"Finish preparing data for new decade");
				TimeWatcher::getInstance().AddTimeSpot(clock(),buf);
			}

			if (time[itime]>finSP)
			{
				break;
			}


			_time64(&timestruct); //  
			time_beg = timestruct;

			msg_iter = 0;

			if(!(give_out->fsdiff))
			{
				t_j2 = t_j1;
				t_j1 = t_j;
				t_j = time[itime];

				tnum2 = tnum1;
				tnum1 = tnum;
				tnum = itime;

				Subdomain::ntimes=tend-tbeg;
				give_out->cur_dec_strt_res=cur_dec_strt;
				if (tend >= give_out->time_layer_last)
				{
					Subdomain::ftout = 1;
				}
				else
				{
					Subdomain::ftout = 2;
				}
			}
			else
			{
				t_j2 = time[i_u_2];
				t_j1 = time[i_u_1];
				t_j = time[itime];

				tnum2 = i_u_2;
				tnum1 = i_u_1;
				tnum = itime;

				Subdomain::ntimes=tbeg-tbeg_pre;
				give_out->cur_dec_strt_res=tbeg_pre;
				Subdomain::ftout=1;
			}

			give_out->vPre1[itime]=tnum1;
			give_out->vPre2[itime]=tnum2;

			if(FirstDecadeStep)
			{
				dt =  t_j  - t_j2;
				dt0 = t_j  - t_j1;
				dt1 = t_j1 - t_j2;

				FirstDecadeStep=false;
			}

			mult  = (dt + dt0)/(dt*dt0);
			mult1 = dt/(dt1*dt0);
			mult2 = -dt0/(dt*dt1);


			if(!fcalc)
			{
				sprintf(buf,"Start prepare matrix for new decade");
				TimeWatcher::getInstance().AddTimeSpot(clock(),buf);

				for(i=0; i<n; i++)
				{
					Slae->di[i] = c_di_all[i]*mult + b_di_all[i];
					for(j=Pt->ig[i]; j<=Pt->ig[i+1]-1; j++)
						Slae->ggl[j] = c_gg_all[j]*mult + b_gg_all[j];
				}

				for(i=0;i<Slae->n;i++){Slae->di[i]+=Slae->di[i]*1e-12;}

				Bc->Set_Dirichlet_Cond(Pt->ig,Pt->jg,Slae->di,Slae->ggl,Slae->pr,time[itime],0);


				sprintf(buf,"Start factorization");
				TimeWatcher::getInstance().AddTimeSpot(clock(),buf);

				prds.factorize(Slae->n,Pt->ig,Pt->jg,Slae->ggl,Slae->di,nthreads);

				sprintf(buf,"Finish factorization");
				TimeWatcher::getInstance().AddTimeSpot(clock(),buf);
			}


				sprintf(buf,"Start loading and adding Enorm");
				TimeWatcher::getInstance().AddTimeSpot(clock(),buf);


				if(npntE0 && !d->fdirect)
				{
					LoadE0FromRZ(itime,En_edges);
				}

				for(ipls=0;ipls<npls;ipls++)
				{
					ipls_cur=ipls;
				
					if(npntE0 && !d->fdirect)
					{
						mult_symmetr(Pt->ig, Pt->jg, c_gg_all_dsig, c_di_all_dsig, En_edges+(Tmap->n*ipls), Slae->pr+(Slae->n*ipls), n);
					}
					else
					{
						for(i=0;i<n;i++){Slae->pr[Slae->n*ipls+i] = 0.0;}
					}

					for(i=0;i<n;i++){cujsum[i]=u_j1[Tmap->n*ipls+i]*mult1+u_j2[Tmap->n*ipls+i]*mult2;}
					mult_symmetr(Pt->ig, Pt->jg, c_gg_all, c_di_all, cujsum, cuj, n);
					for(i=0; i<n; i++){Slae->pr[Slae->n*ipls+i] += cuj[i];}
				}

				sprintf(buf,"Finish loading and adding Enorm");
				TimeWatcher::getInstance().AddTimeSpot(clock(),buf);

				if(d->fdirect && d->nSrsRazb)
				{
					if(d->currentfunction[itime])
					{
						for(ipls=0;ipls<npls;ipls++)
						{
							m=(int)d->srs[ipls].size();
							for(l=0;l<m;l++)
							{
								loc_source &ls=d->srs[ipls][l];
								
								binn[0]=2.0*ls.ksi-1.0;
								binn[1]=2.0*ls.eta-1.0;
								binn[2]=2.0*ls.dzeta-1.0;

								for(j=0; j<12; j++)
								{
									int ii,i_mu,it,nn0,nn1;
									double edlen;
									
									ii=Tmap->ed[ls.elem][j];

									nn0=Tmap->edges[ii][0];
									nn1=Tmap->edges[ii][1];

									edlen=sqrt((d->xyz[nn1][0]-d->xyz[nn0][0])*(d->xyz[nn1][0]-d->xyz[nn0][0])+
										(d->xyz[nn1][1]-d->xyz[nn0][1])*(d->xyz[nn1][1]-d->xyz[nn0][1])+
										(d->xyz[nn1][2]-d->xyz[nn0][2])*(d->xyz[nn1][2]-d->xyz[nn0][2]));

									brk.Basis_func_on_reference_vec_par(j,binn,bout);
									value=(bout[0]*ls.dx+bout[1]*ls.dy+bout[2]*ls.dz)*d->currentfunction[itime]*2.0/edlen;
					
									Slae->pr[Slae->n*ipls+ii] += value;
								}
							}
						}
					}
				}

				for(ipls=0;ipls<npls;ipls++)
				{
					Bc->Set_Dirichlet_Cond_Pr(Slae->pr,time[itime],0,ipls);
				}


				sprintf(buf,"Start Direct-Inverse");
				TimeWatcher::getInstance().AddTimeSpot(clock(),buf);

				prds.solve_nrhs(npls,Slae->pr,u_j);


				sprintf(buf,"Finish Direct-Inverse");
				TimeWatcher::getInstance().AddTimeSpot(clock(),buf);

				if(itime>=nlst)
				{
					_itoa_s(itime, buf, 10);
					strcpy(fname, v3);
					strcat(fname, buf);

					ofp.open(fname,ios::binary);
					ofp.write((char *)u_j,npls*size_d*Slae->n);
					ofp.close();
					ofp.clear();
				}
			
			if(itime==_ntime-1)
			{
				ofp.open("a.last",ios::binary);
				ofp.write((char *)u_j,npls*size_d*Slae->n);
				ofp.close();
				ofp.clear();
			}

			sprintf(buf,"Start Resultant");
			TimeWatcher::getInstance().AddTimeSpot(clock(),buf);

			for(ipls=0;ipls<npls;ipls++)
			{
				give_out->ipls_cur=ipls_cur=ipls;

				if(tnum1 < give_out->time_layer_last)
				{
					give_out->Work2((u_j2+ipls*Tmap->n), u_r=(u_j1+ipls*Tmap->n), (u_j+ipls*Tmap->n), tnum2, tnum1, tnum, tnum1);
				}
				if(tnum == give_out->time_layer_last)
				{
					give_out->Work2((u_j2+ipls*Tmap->n), (u_j1+ipls*Tmap->n), u_r=(u_j+ipls*Tmap->n), tnum2, tnum1, tnum, tnum);
				}
			}

			sprintf(buf,"Finish Resultant");
			TimeWatcher::getInstance().AddTimeSpot(clock(),buf);

			for(i=0; i<npls*Tmap->n; i++)
			{
				u_j2[i] = u_j1[i];
				u_j1[i] = u_j[i];
			}


			cout<<"time layer - "<<itime<<" done"<<endl;
			fcalc=true;
		}
	}

	fclose(f_time);



	prds.stop_solver();

	return 0;
}

void Time_approx_for_vfem::FinishCalculation()
{
	if(!CheckStop())
	{
		int retv;

		if(d->n_pointresB)
		{
			retv=give_out->resultantB->StopSolvers();
			if(retv)
			{
				cout<<"Error in function "<<"resultantB->StopSolvers()"<<endl;
				CloseProgramm(retv);
			}
		}

		if(d->n_pointresE)
		{
			retv=give_out->resultantA->StopSolvers();
			if(retv)
			{
				cout<<"Error in function "<<"resultantA->StopSolvers()"<<endl;
				CloseProgramm(retv);
			}
		}

		logfile << "Begin Gather_edsall...\n";
		cout << "Begin Gather_edsall...\n";
		fflush(stdout);
		give_out->Gather_edsall(true);
		logfile << "End Gather_edsall.\n";
	}
}
int Time_approx_for_vfem::Projection_from_nodes_to_edges(int n_edges_c, 
	double (*En_nodes)[3], double *En_edges, double (*xyz)[3],
	int (*edges)[2], bool *is_field_in_node, int n_elem, int (*nver)[14],
	int (*ed)[25])
{
	int i, k;
	int v1, v2;
	double mid_value[3];
	int edge;

	double axis[3];
	double len;
	for (k=0; k<n_anomal_edges; k++)
	{
		edge = anomal_edges[k];
		v1 = edges[edge][0];
		v2 = edges[edge][1];

		for (i=0; i<3; i++)
		{
			axis[i] = xyz[v2][i] - xyz[v1][i];
			mid_value[i] = 0.5*(En_nodes[ipls_cur*d->kuzlov+v1][i] + En_nodes[ipls_cur*d->kuzlov+v2][i]);
		}

		len = Norm_Euclid(axis, 3);
		En_edges[edge] = Projection_On_Axis(mid_value, axis)*len*0.5;
	}

	return 0;
}

int Time_approx_for_vfem::OpenEnorm(int nlayer,ifstream &fp)
{
	char fname[30];
	char buf[20];
	char buffer2[20];
	const char en[] = {"enor."};

	_itoa_s(nlayer, buf, 10);

	switch(strlen(buf))
	{
	case 1:
		strcpy(buffer2, "000");
		strcat(buffer2, buf);
		strcpy(buf, buffer2);
		break;
	case 2:
		strcpy(buffer2, "00");
		strcat(buffer2, buf);
		strcpy(buf, buffer2);
		break;
	case 3:
		strcpy(buffer2, "0");
		strcat(buffer2, buf);
		strcpy(buf, buffer2);
	}

	strcpy(fname, en);
	strcat(fname, buf);

	fp.open(fname, ios::binary);
	if(!fp)
	{
		return 1;
	}

	return 0;
}
wchar_t *convertCharArrayToLPCWSTR(const char* charArray)
{
	wchar_t* wString = new wchar_t[4096];
	MultiByteToWideChar(CP_ACP, 0, charArray, -1, wString, 4096);
	return wString;
}

int Time_approx_for_vfem::LoadE0FromRZ(int tnum,double *En_edges)
{
	ifstream inf;
	int i,j,v1,v2,ipls;
	double len,axis[3];
	char buf[256],PName[256];
	float tmpf;
	inf.open("PName");
	if(!inf)
	{
		cout << "Error in open file " << "PName" << '\n';
		logfile << "Error in open file " << "PName" << '\n';
		return 1;
	}
	inf>>PName;
	inf.close();
	inf.clear();
	if(!AnomalType)
	{
		sprintf(buf,"m%s.enor.%04d",PName,tnum);
	}
	else
	{
		sprintf(buf,"m%s.enor.%04d_%d",PName,tnum,ia);
	}

	int MemSize;
	HANDLE hMapFile;
	LPCTSTR pBuf;
	int nnn,fff;
	const int size_f=sizeof(float);

	MemSize=0;
	for(i=0;i<Tmap->n;i++)
	{
		MemSize+=isEdgesAnomal[i];
	}
	MemSize*=npls*size_f;

	if(MemSize)
	{
		if(!Buff)
		{
			Buff=new char[MemSize];
		}

		hMapFile = OpenFileMapping(FILE_MAP_READ,TRUE,convertCharArrayToLPCWSTR(buf));
		if (hMapFile == NULL || hMapFile == INVALID_HANDLE_VALUE)
		{
			printf("        (%d).\n", GetLastError());
			return 1;
		}
		pBuf = (LPTSTR)MapViewOfFile(hMapFile, FILE_MAP_READ, 0, 0, MemSize);
		if (pBuf == NULL)
		{
			printf("     (%d).\n", GetLastError());
			return 1;
		}

		CopyMemory(Buff, (PVOID)pBuf, MemSize);

		UnmapViewOfFile(pBuf);
		CloseHandle(hMapFile);
	}

	nnn=0;
	for(ipls=0;ipls<npls;ipls++)
	{
		for(i=0;i<Tmap->n;i++)
		{
			En_edges[i+ipls*Tmap->n]=0.0;
			if(isEdgesAnomal[i])
			{
				for(fff=0;fff<size_f;fff++){*(((char *)&tmpf)+fff)=Buff[nnn+fff];}
				nnn+=size_f;
				v1=Tmap->edges[i][0];
				v2=Tmap->edges[i][1];
				for(j=0;j<3;j++)
				{
					axis[j]=d->xyz[v2][j]-d->xyz[v1][j];
				}
				len = Norm_Euclid(axis, 3);
				En_edges[i+ipls*Tmap->n]=tmpf*len*0.5;
			}
		}
	}


	return 0;
}
double Time_approx_for_vfem::Calc_dof(double *J, double *func, int n_local_edge)
{
	double Jt[3];

	Mult_Plot(J, (double*)TANGENT_VECTORS_ON_REFERENCE_CUBE[n_local_edge], Jt, 3);

	return Scal(func, Jt, 3);
}
double Time_approx_for_vfem::dA_dt(double t, double u_j, double u_j1, double u_j2,
								double dt, double dt0, double dt1, double t_j, double t_j1, double t_j2)
{
	double du_dt;

	if(t < t_j2) t = t_j2;
	if(t > t_j)  t = t_j;

	du_dt = u_j2*(2.0*t - t_j  - t_j1)/(dt1*dt) - 
		u_j1*(2.0*t - t_j  - t_j2)/(dt1*dt0) + 
		u_j*(2.0*t - t_j1 - t_j2)/(dt*dt0);

	return du_dt;
}
int Time_approx_for_vfem::ReadV3(int nlayer, double *ax, double *ay, double *az)
{
	int i, j, k;
	double *v3=NULL;
	In_Out R;
	char fname[20];
	double x[8],y[8],z[8];
	double ves[12];
	int ipls;

	if(nlayer==0)
	{
		v3=new double[d->kuzlov*3*npls];
		if (forLine)
			R.Read_Bin_File_Of_Double("a0.dat", v3, d->kuzlov*npls, 3);
		else
			for(i=0;i<d->kuzlov*3*npls;i++){v3[i]=0;}
		k=0;
		for(i=0;i<d->kpar;i++)
		{
			if(!isElemAnomal[i])
				continue;
			for(j=0;j<8;j++)
			{
				for(ipls=0;ipls<npls;ipls++)
				{
					ax[k*8+j+ipls*nAnomalElem*8]=v3[d->nver[i][j]*3+ipls*d->kuzlov*3];
					ay[k*8+j+ipls*nAnomalElem*8]=v3[d->nver[i][j]*3+1+ipls*d->kuzlov*3];
					az[k*8+j+ipls*nAnomalElem*8]=v3[d->nver[i][j]*3+2+ipls*d->kuzlov*3];
				}
			}
			k++;
		}
	} 
	else
	{
		v3=new double[Tmap->n*npls];

		sprintf_s(fname, "v3.%d", nlayer);
		R.Read_Bin_File_Of_Double(fname, v3, Tmap->n_c*npls, 1);

		k=0;
		for(i=0;i<d->kpar;i++)
		{
			if(!isElemAnomal[i])
				continue;
			for (j=0; j<8; j++)
			{
				x[j]=d->xyz[d->nver[i][j]][0];
				y[j]=d->xyz[d->nver[i][j]][1];
				z[j]=d->xyz[d->nver[i][j]][2];
			}
			T_Brick hex(x,y,z);
			for(ipls=0;ipls<npls;ipls++)
			{
				for(j=0;j<12;j++){ves[j]=v3[Tmap->ed[i][j]+ipls*Tmap->n_c];}
				hex.GetVectorFieldNodes(ves,&ax[k*8+ipls*nAnomalElem*8],&ay[k*8+ipls*nAnomalElem*8],&az[k*8+ipls*nAnomalElem*8]);
			}
			k++;
		}
	}

	if(v3) {delete [] v3; v3=NULL;}

	return 0;
}
int Time_approx_for_vfem::GetNumberOfNodes()
{
	return d->kuzlov;
}
int Time_approx_for_vfem::GetNumberOfElements()
{
	return d->kpar;
}
int Time_approx_for_vfem::GetElementNodesNumber()
{
	return 8;
}
const pv::Point3D Time_approx_for_vfem::GetNode(const int& i_node)
{
	pv::Point3D Point;
	Point.x()=d->xyz[i_node][0];
	Point.y()=d->xyz[i_node][1];
	Point.z()=d->xyz[i_node][2];
	return Point;
}
const pv::Point3D Time_approx_for_vfem::GetNodeTrue(const int& i_node)
{
	pv::Point3D Point;
	Point.x()=d->xyzt[i_node][0];
	Point.y()=d->xyzt[i_node][1];
	Point.z()=d->xyzt[i_node][2];
	return Point;
}
int Time_approx_for_vfem::GetNodeNumberOnElement(const int& i_element, const int& i_node)
{
	return d->nver[i_element][i_node];
}
int Time_approx_for_vfem::GetElementMaterial(const int& i_element)
{
	return d->nvkat[i_element];
}
int Time_approx_for_vfem::GetTypeOfElement(const int& i_element)
{
	return d->nver[i_element][13];
}
double Time_approx_for_vfem::GetValueInElementCenter(const int& i_element, const Res3DValueType& r_type)
{
	int i;
	double f[12];
	double x[8],y[8],z[8];
	double in[3],out[3],FieldOnElem;

	FieldOnElem=0;

	for (i=0; i<12; i++)f[i] = u_r[Slae->ed[i_element][i]];

	for (i=0; i<8; i++){
		x[i] = d->xyz[d->nver[i_element][i]][0];
		y[i] = d->xyz[d->nver[i_element][i]][1];
		z[i] = d->xyz[d->nver[i_element][i]][2];
	}
	
	T_Brick L(x, y, z);

	in[0]=in[1]=in[2]=0.0;

	if(r_type==vtAx ||r_type==vtAy || r_type==vtAz){
		L.Calc_value_inside_hex(f,in,out);

		if(r_type==vtAx)FieldOnElem=out[0];
		if(r_type==vtAy)FieldOnElem=out[1];
		if(r_type==vtAz)FieldOnElem=out[2];
	}
	else{
		L.Calc_rotor_inside_hex(f,in,out);

		if(r_type==vtRotxA)FieldOnElem=out[0];
		if(r_type==vtRotyA)FieldOnElem=out[1];
		if(r_type==vtRotzA)FieldOnElem=out[2];
	}

	return FieldOnElem;
}
int Time_approx_for_vfem::GetNumberOfResPoints(const Res3DValueType& r_type)
{
	return (r_type==vtWithDiscontinuity)? d->n_pointresE : d->n_pointresB;
}
pv::Point3D Time_approx_for_vfem::GetResPoint(const Res3DValueType& r_type, const int& i_point)
{
	pv::Point3D Point;
	Point.x()=(r_type==vtWithDiscontinuity)? d->pointresE[i_point][0] : d->pointresB[i_point][0];
	Point.y()=(r_type==vtWithDiscontinuity)? d->pointresE[i_point][1] : d->pointresB[i_point][1];
	Point.z()=(r_type==vtWithDiscontinuity)? d->pointresE[i_point][2] : d->pointresB[i_point][2];
	return Point;
}
int * Time_approx_for_vfem::GetPointerToRegular()
{
	return &(d->regular[0]);
}
int Time_approx_for_vfem::GetXSize()
{
	return d->N_X;
}
int Time_approx_for_vfem::GetYSize()
{
	return d->N_Y;
}
int Time_approx_for_vfem::GetZSize()
{
	return d->N_Z;
}

double * Time_approx_for_vfem::GetPointerToX()
{
	return &(d->Xcrd[0]);
}

double * Time_approx_for_vfem::GetPointerToY()
{
	return &(d->Ycrd[0]);
}

double * Time_approx_for_vfem::GetPointerToZ()
{
	return &(d->Zcrd[0]);
}
void Time_approx_for_vfem::SaveResult(const Res3DValueType& r_type, const double& r_value, const int& j, const int& tnum, int ipls)
{
	if(tnum >= give_out->time_layer_first && tnum<=give_out->time_layer_last)
	{
		if(r_type==vtRotxA)
				give_out->Bx_all[((tnum - give_out->time_layer_first)*npls + ipls)*d->n_pointresB + j]=r_value;
		else if(r_type==vtRotyA)
				give_out->By_all[((tnum - give_out->time_layer_first)*npls + ipls)*d->n_pointresB + j]=r_value;
		else if(r_type==vtRotzA)
				give_out->Bz_all[((tnum - give_out->time_layer_first)*npls + ipls)*d->n_pointresB + j]=r_value;
		else if(r_type==vtAx)
				give_out->Ex_all[((tnum - give_out->time_layer_first)*npls + ipls)*d->n_pointresE + j]=r_value;
		else if(r_type==vtAy)
				give_out->Ey_all[((tnum - give_out->time_layer_first)*npls + ipls)*d->n_pointresE + j]=r_value;
		else if(r_type==vtAz)
				give_out->Ez_all[((tnum - give_out->time_layer_first)*npls + ipls)*d->n_pointresE + j]=r_value;
		else
		{
			cout<<"Unknown Res3DValueType in Time_approx_for_vfem::SaveResult"<<'\n';
			logfile<<"Unknown Res3DValueType in Time_approx_for_vfem::SaveResult"<<'\n';
		}
	}
}

void Time_approx_for_vfem::WriteTemporaryAllFile(int nrec,vector<int> &RecToSource,int tbeg,int tend,double *FieldNorm,double *FieldAnom,char *fname)
{
	ofstream ofp;
	int i,j,ind,ipls;
	double tmpd;
	ofp.open(fname,ios::binary);
		for(i=0;i<nrec;i++)
		{
		ipls=RecToSource[i];

			for(j=tbeg;j<tend;j++)
			{
				if(j>=give_out->time_layer_first && j<=give_out->time_layer_last)
				{
					ind = ((j - give_out->time_layer_first)*npls+ipls)*nrec + i;	
					tmpd=FieldNorm[ind];
					ofp.write((char *)&tmpd,size_d);
					tmpd=FieldAnom[ind];
					ofp.write((char *)&tmpd,size_d);
				}
		}
	}
	ofp.close();
	ofp.clear();
}

void Time_approx_for_vfem::SaveCurrentResults(int iDec,int tbeg,int tend,vector<int> &RecToSourceB,vector<int> &RecToSourceE)
{
	char str[32];
	
	sprintf(str,"bxall_dec_%d",iDec+1);
	WriteTemporaryAllFile(d->n_pointresB,RecToSourceB,tbeg,tend,give_out->Bxn_all,give_out->Bx_all,str);
	sprintf(str,"byall_dec_%d",iDec+1);
	WriteTemporaryAllFile(d->n_pointresB,RecToSourceB,tbeg,tend,give_out->Byn_all,give_out->By_all,str);
	sprintf(str,"bzall_dec_%d",iDec+1);
	WriteTemporaryAllFile(d->n_pointresB,RecToSourceB,tbeg,tend,give_out->Bzn_all,give_out->Bz_all,str);

	sprintf(str,"axall_dec_%d",iDec+1);
	WriteTemporaryAllFile(d->n_pointresE,RecToSourceE,tbeg,tend,give_out->Exn_all,give_out->Ex_all,str);
	sprintf(str,"ayall_dec_%d",iDec+1);
	WriteTemporaryAllFile(d->n_pointresE,RecToSourceE,tbeg,tend,give_out->Eyn_all,give_out->Ey_all,str);
	sprintf(str,"azall_dec_%d",iDec+1);
	WriteTemporaryAllFile(d->n_pointresE,RecToSourceE,tbeg,tend,give_out->Ezn_all,give_out->Ez_all,str);
}

int Time_approx_for_vfem::ReadTemporaryAllFile(int nrec,vector<int> &RecToSource,int tbeg,int tend,double *FieldNorm,double *FieldAnom,char *fname)
{
	ifstream inf;
	int i,j,ind,ipls;
	double tmpd;
	inf.open(fname,ios::binary);
	if(!inf)
	{
		logfile<<"Error: can't open file "<<fname<<endl;
		cout<<"Error: can't open file "<<fname<<endl;
		return 1;
	}
		for(i=0;i<nrec;i++)
		{
		ipls=RecToSource[i];

			for(j=tbeg;j<tend;j++)
			{
				if(j>=give_out->time_layer_first && j<=give_out->time_layer_last)
				{
					ind = ((j - give_out->time_layer_first)*npls+ipls)*nrec + i;
					inf.read((char *)&tmpd,size_d);
					FieldNorm[ind]=tmpd;
					inf.read((char *)&tmpd,size_d);
					FieldAnom[ind]=tmpd;
				}
			}
		}
	inf.close();
	inf.clear();
	return 0;
}

int Time_approx_for_vfem::LoadPreviousResults(int iDec,int tbeg,int tend,vector<int> &RecToSourceB,vector<int> &RecToSourceE)
{
	char str[32];
	int retv;
	
	sprintf(str,"bxall_dec_%d",iDec+1);
	retv=ReadTemporaryAllFile(d->n_pointresB,RecToSourceB,tbeg,tend,give_out->Bxn_all,give_out->Bx_all,str);
	if(retv){return retv;}
	sprintf(str,"byall_dec_%d",iDec+1);
	retv=ReadTemporaryAllFile(d->n_pointresB,RecToSourceB,tbeg,tend,give_out->Byn_all,give_out->By_all,str);
	if(retv){return retv;}
	sprintf(str,"bzall_dec_%d",iDec+1);
	retv=ReadTemporaryAllFile(d->n_pointresB,RecToSourceB,tbeg,tend,give_out->Bzn_all,give_out->Bz_all,str);
	if(retv){return retv;}

	sprintf(str,"axall_dec_%d",iDec+1);
	retv=ReadTemporaryAllFile(d->n_pointresE,RecToSourceE,tbeg,tend,give_out->Exn_all,give_out->Ex_all,str);
	if(retv){return retv;}
	sprintf(str,"ayall_dec_%d",iDec+1);
	retv=ReadTemporaryAllFile(d->n_pointresE,RecToSourceE,tbeg,tend,give_out->Eyn_all,give_out->Ey_all,str);
	if(retv){return retv;}
	sprintf(str,"azall_dec_%d",iDec+1);
	retv=ReadTemporaryAllFile(d->n_pointresE,RecToSourceE,tbeg,tend,give_out->Ezn_all,give_out->Ez_all,str);
	if(retv){return retv;}

	return 0;
}

int Time_approx_for_vfem::InitvPre(int tbeg,int tend,int i_u_1,int i_u_2)
{
	int itime,tnum2,tnum1,tnum;
	bool fcalc;
	
	int _ntime = ntime;

	if(1>=tbeg && 1<tend)
	{
	for(itime=0; itime<1; itime++)
	{
		give_out->vPre1[itime]=itime-1;
		give_out->vPre2[itime]=itime-2;
	}

	tnum1 = 0;
	tnum = 1;
	}

	fcalc=false;

	for(itime=2; itime<_ntime; itime++)
	{
		give_out->fsdiff=false;

		if(itime>=tbeg && itime<tend)
		{
		if(!fcalc && itime!=2)
		{
			give_out->fsdiff=true;
		}

		if(!(give_out->fsdiff))
		{
			tnum2 = tnum1;
			tnum1 = tnum;
			tnum = itime;
		}
		else
		{
			tnum2 = i_u_2;
			tnum1 = i_u_1;
			tnum = itime;
		}

		give_out->vPre1[itime]=tnum1;
		give_out->vPre2[itime]=tnum2;

		fcalc=true;
		}
	}

	return 0;
}
