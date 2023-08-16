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
 *  This file contains the code for outputting E and EMF fields from solution in 3D VFEM
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin                         
 *  Novosibirsk State Technical University,                 
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia       
 *  vdv_wk@mail.ru                                          
 *  Version 2.0 January 16, 2023                            
*/

#include "stdafx.h" 
#include "Give_out_vec_loop.h"
extern ofstream logfile;
extern void Memory_allocation_error(const char *var, const char *func);
extern void Cannot_open_file(const char *fname, const char *func);
extern bool CheckStop();
Give_out_vec_loop::Give_out_vec_loop(Vec_Prep_Data *d, T_Mapping_Vec *tmap,
									 vector<int> &_RecvPlsIgB,vector<int> &_RecvPlsIgE, int _npls)
{
	FILE *fp=NULL;

	resultantA=NULL;
	resultantB=NULL;

	int i;

	this->d = d;
	this->tmap = tmap;

	npls=_npls;



	forLine=true;

	ifstream inf;
	inf.open("timeintervalforprint");
	inf>>time_interval_for_print_left>>time_interval_for_print_right;
	inf.close();


	time_layer_first = time_layer_last = 2; //  0  1     
	for (i=0; i<d->ntime; i++)
	{
		if(d->time[i] >= time_interval_for_print_left)
		{
			time_layer_first = i;
			break;
		}
	}
	if(time_layer_first > d->ntime-1)
		time_layer_first = d->ntime-1;

	for (i=0; i<d->ntime; i++)
	{
		if(d->time[i] >= time_interval_for_print_right-1e-6)
		{
			time_layer_last = i;
			break;
		}
	}
	if(time_interval_for_print_right>d->time[d->ntime-1]-1e-6 )
		time_layer_last = d->ntime-1;

	Edsx_all = NULL;
	Edsy_all = NULL;
	Edsz_all = NULL;

	Bx_all = NULL;
	By_all = NULL;
	Bz_all	= NULL;

	Edsxn_all = NULL;
	Edsyn_all = NULL;
	Edszn_all = NULL;

	Bxn_all = NULL;
	Byn_all = NULL;
	Bzn_all = NULL;

	if(d->n_pointresB)
	{
	if((Edsx_all = new double[npls*(time_layer_last - time_layer_first + 1)*d->n_pointresB])==NULL)
		Memory_allocation_error("Edsx_all", "Give_out_vec_loop::Give_out_vec_loop");
	Clear(Edsx_all, npls*(time_layer_last - time_layer_first + 1)*d->n_pointresB);

	if((Edsy_all = new double[npls*(time_layer_last - time_layer_first + 1)*d->n_pointresB])==NULL)
		Memory_allocation_error("Edsy_all", "Give_out_vec_loop::Give_out_vec_loop");
	Clear(Edsy_all, npls*(time_layer_last - time_layer_first + 1)*d->n_pointresB);

	if((Edsz_all = new double[npls*(time_layer_last - time_layer_first + 1)*d->n_pointresB])==NULL)
		Memory_allocation_error("Edsz_all", "Give_out_vec_loop::Give_out_vec_loop");
	Clear(Edsz_all, npls*(time_layer_last - time_layer_first + 1)*d->n_pointresB);


	if((Bx_all = new double[npls*(time_layer_last - time_layer_first + 1)*d->n_pointresB])==NULL)
		Memory_allocation_error("Bx_all", "Give_out_vec_loop::Give_out_vec_loop");
	Clear(Bx_all, npls*(time_layer_last - time_layer_first + 1)*d->n_pointresB);

	if((By_all = new double[npls*(time_layer_last - time_layer_first + 1)*d->n_pointresB])==NULL)
		Memory_allocation_error("By_all", "Give_out_vec_loop::Give_out_vec_loop");
	Clear(By_all, npls*(time_layer_last - time_layer_first + 1)*d->n_pointresB);

	if((Bz_all = new double[npls*(time_layer_last - time_layer_first + 1)*d->n_pointresB])==NULL)
		Memory_allocation_error("Bz_all", "Give_out_vec_loop::Give_out_vec_loop");
	Clear(Bz_all, npls*(time_layer_last - time_layer_first + 1)*d->n_pointresB);


	if((Edsxn_all = new double[npls*(time_layer_last-time_layer_first+1)*d->n_pointresB])==0)
		Memory_allocation_error("Edsxn_all", "Give_out_vec_loop::Give_out_vec_loop");
	Clear(Edsxn_all, npls*(time_layer_last - time_layer_first + 1)*d->n_pointresB);

	if((Edsyn_all = new double[npls*(time_layer_last-time_layer_first+1)*d->n_pointresB])==0)
		Memory_allocation_error("Edsyn_all", "Give_out_vec_loop::Give_out_vec_loop");
	Clear(Edsyn_all, npls*(time_layer_last - time_layer_first + 1)*d->n_pointresB);

	if((Edszn_all = new double[npls*(time_layer_last-time_layer_first+1)*d->n_pointresB])==0)
		Memory_allocation_error("Edszn_all", "Give_out_vec_loop::Give_out_vec_loop");
	Clear(Edszn_all, npls*(time_layer_last - time_layer_first + 1)*d->n_pointresB);

	if((Bxn_all = new double[npls*(time_layer_last-time_layer_first+1)*d->n_pointresB])==0)
		Memory_allocation_error("Bxn_all", "Give_out_vec_loop::Give_out_vec_loop");
	Clear(Bxn_all, npls*(time_layer_last - time_layer_first + 1)*d->n_pointresB);

	if((Byn_all = new double[npls*(time_layer_last-time_layer_first+1)*d->n_pointresB])==0)
		Memory_allocation_error("Byn_all", "Give_out_vec_loop::Give_out_vec_loop");
	Clear(Byn_all, npls*(time_layer_last - time_layer_first + 1)*d->n_pointresB);

	if((Bzn_all = new double[npls*(time_layer_last-time_layer_first+1)*d->n_pointresB])==0)
		Memory_allocation_error("Bzn_all", "Give_out_vec_loop::Give_out_vec_loop");
	Clear(Bzn_all, npls*(time_layer_last - time_layer_first + 1)*d->n_pointresB);
	}


	Ex_all = NULL;
	Ey_all = NULL;
	Ez_all	= NULL;

	Exn_all = NULL;
	Eyn_all = NULL;
	Ezn_all = NULL;


	if(d->n_pointresE)
	{
		if((Ex_all = new double[npls*(time_layer_last-time_layer_first+1)*d->n_pointresE])==0)
			Memory_allocation_error("Ex_all", "Give_out_vec_loop::Give_out_vec_loop");
		Clear(Ex_all, npls*(time_layer_last - time_layer_first + 1)*d->n_pointresE);

		if((Ey_all = new double[npls*(time_layer_last-time_layer_first+1)*d->n_pointresE])==0)
			Memory_allocation_error("Ey_all", "Give_out_vec_loop::Give_out_vec_loop");
		Clear(Ey_all, npls*(time_layer_last - time_layer_first + 1)*d->n_pointresE);

		if((Ez_all = new double[npls*(time_layer_last-time_layer_first+1)*d->n_pointresE])==0)
			Memory_allocation_error("Ez_all", "Give_out_vec_loop::Give_out_vec_loop");
		Clear(Ez_all, npls*(time_layer_last - time_layer_first + 1)*d->n_pointresE);


		if((Exn_all = new double[npls*(time_layer_last-time_layer_first+1)*d->n_pointresE])==0)
			Memory_allocation_error("Exn_all", "Give_out_vec_loop::Give_out_vec_loop");
		Clear(Exn_all, npls*(time_layer_last - time_layer_first + 1)*d->n_pointresE);

		if((Eyn_all = new double[npls*(time_layer_last-time_layer_first+1)*d->n_pointresE])==0)
			Memory_allocation_error("Eyn_all", "Give_out_vec_loop::Give_out_vec_loop");
		Clear(Eyn_all, npls*(time_layer_last - time_layer_first + 1)*d->n_pointresE);

		if((Ezn_all = new double[npls*(time_layer_last-time_layer_first+1)*d->n_pointresE])==0)
			Memory_allocation_error("Ezn_all", "Give_out_vec_loop::Give_out_vec_loop");
		Clear(Ezn_all, npls*(time_layer_last - time_layer_first + 1)*d->n_pointresE);
	}

	RecvPlsIgB=_RecvPlsIgB;
	RecvPlsIgE=_RecvPlsIgE;
}
Give_out_vec_loop::~Give_out_vec_loop()
{
	if(resultantA){delete resultantA;resultantA = NULL;}
	if(resultantB){delete resultantB;resultantB = NULL;}

	if(Edsxn_all) {delete [] Edsxn_all; Edsxn_all=NULL;}
	if(Edsyn_all) {delete [] Edsyn_all; Edsyn_all=NULL;}
	if(Edszn_all) {delete [] Edszn_all; Edszn_all=NULL;}
	
	if(Bxn_all) {delete [] Bxn_all; Bxn_all=NULL;}
	if(Byn_all) {delete [] Byn_all; Byn_all=NULL;}
	if(Bzn_all) {delete [] Bzn_all; Bzn_all=NULL;}

	if(Edsx_all) {delete [] Edsx_all; Edsx_all=NULL;}
	if(Edsy_all) {delete [] Edsy_all; Edsy_all=NULL;}
	if(Edsz_all) {delete [] Edsz_all; Edsz_all=NULL;}

	if(Bx_all) {delete [] Bx_all; Bx_all=NULL;}
	if(By_all) {delete [] By_all; By_all=NULL;}
	if(Bz_all) {delete [] Bz_all; Bz_all=NULL;}


	if(Exn_all) {delete [] Exn_all; Exn_all=NULL;}
	if(Eyn_all) {delete [] Eyn_all; Eyn_all=NULL;}
	if(Ezn_all) {delete [] Ezn_all; Ezn_all=NULL;}

	if(Ex_all) {delete [] Ex_all; Ex_all=NULL;}
	if(Ey_all) {delete [] Ey_all; Ey_all=NULL;}
	if(Ez_all) {delete [] Ez_all; Ez_all=NULL;}
}
void Give_out_vec_loop::Clear(double *u, const int& usize)
{
	for (int i=0; i<usize; i++)
		u[i]=0;
}
void Give_out_vec_loop::Work2(double *u_j2, double *u_j1, double *u_j,
							  int t_2, int t_1, int t_0, int t)
{
	this->v3_j2 = u_j2;
	this->v3_j1 = u_j1;
 	this->v3_j = u_j;

	tnum0 = t_2;
	tnum1 = t_1;
	tnum2 = t_0;

	tnum = t;

	t_j2 = d->time[tnum0];
	t_j1 = d->time[tnum1];
	t_j = d->time[tnum2];

	if(t >= time_layer_first && t<= time_layer_last)
		Work_on_current_time_layer();
}
void Give_out_vec_loop::Gather_edsall(bool forCED)
{
	int i, j, j_plus_1, j_plus_2;
	ofstream outf, outf2;
	double sqrloop=1/*NormResFromLoop->nloops*_PI_*NormResFromLoop->rloop*NormResFromLoop->rloop*/; //   
	int ind;

	int id[3];
	double dt, dt0, dt1;
	double mt0,mt1,mt2;
	double F[3][3];

	int ipls;
	char fname[256];

	double *Ax_all,*Ay_all,*Az_all;

	Ax_all=Ay_all=Az_all=NULL;

	if(d->n_pointresE)
	{
		if((Ax_all = new double[(time_layer_last-time_layer_first+1)*npls*d->n_pointresE])==0)
			Memory_allocation_error("Ex_all", "Give_out_vec_loop::Give_out_vec_loop");
		Clear(Ax_all, (time_layer_last - time_layer_first + 1)*d->n_pointresE);

		if((Ay_all = new double[(time_layer_last-time_layer_first+1)*npls*d->n_pointresE])==0)
			Memory_allocation_error("Ey_all", "Give_out_vec_loop::Give_out_vec_loop");
		Clear(Ay_all, (time_layer_last - time_layer_first + 1)*d->n_pointresE);

		if((Az_all = new double[(time_layer_last-time_layer_first+1)*npls*d->n_pointresE])==0)
			Memory_allocation_error("Ez_all", "Give_out_vec_loop::Give_out_vec_loop");
		Clear(Az_all, (time_layer_last - time_layer_first + 1)*d->n_pointresE);

		for (i=0; i<((time_layer_last-time_layer_first+1))*npls*d->n_pointresE; i++)
		{
			Ax_all[i]=Ex_all[i];
			Ay_all[i]=Ey_all[i];
			Az_all[i]=Ez_all[i];
		}
	}

	
	for(ipls=0;ipls<npls;ipls++)
	{
		for(i=0;i<d->n_pointresB;i++)
		{
			if(i<RecvPlsIgB[ipls] || i>=RecvPlsIgB[ipls+1]){continue;}

			j=time_layer_first;
			id[0]=((j - time_layer_first)*npls+ipls)*d->n_pointresB + i;

			j_plus_1=j+1;
			id[1]=((j_plus_1 - time_layer_first)*npls+ipls)*d->n_pointresB + i;

			j_plus_2=j_plus_1+1;
			id[2]=((j_plus_2 - time_layer_first)*npls+ipls)*d->n_pointresB + i;

			F[0][0]=Bx_all[id[0]];
			F[0][1]=Bx_all[id[1]];
			F[0][2]=Bx_all[id[2]];

			F[1][0]=By_all[id[0]];
			F[1][1]=By_all[id[1]];
			F[1][2]=By_all[id[2]];
			
			F[2][0]=Bz_all[id[0]];
			F[2][1]=Bz_all[id[1]];
			F[2][2]=Bz_all[id[2]];

			dt=d->time[j_plus_2]-d->time[j];
			dt1=d->time[j_plus_1]-d->time[j];
			dt0=d->time[j_plus_2]-d->time[j_plus_1];

			mt0=-(dt+dt1)/(dt1*dt);
			mt1=(dt)/(dt1*dt0);
			mt2=(-dt1)/(dt*dt0);

			Edsx_all[id[0]]=(F[0][0]*mt0+F[0][1]*mt1+F[0][2]*mt2);
			Edsy_all[id[0]]=(F[1][0]*mt0+F[1][1]*mt1+F[1][2]*mt2);
			Edsz_all[id[0]]=(F[2][0]*mt0+F[2][1]*mt1+F[2][2]*mt2);

			do{
				mt0=(-dt0)/(dt1*dt);
				mt1=-(dt1-dt0)/(dt1*dt0);
				mt2=(dt1)/(dt*dt0);
			
				Edsx_all[id[1]]=(F[0][0]*mt0+F[0][1]*mt1+F[0][2]*mt2);
				Edsy_all[id[1]]=(F[1][0]*mt0+F[1][1]*mt1+F[1][2]*mt2);
				Edsz_all[id[1]]=(F[2][0]*mt0+F[2][1]*mt1+F[2][2]*mt2);

				if(j_plus_2==time_layer_last)break;

				j=j_plus_1;
				j_plus_2++;

				j_plus_1=vPre1[j_plus_2];
				


				id[2]=((j_plus_2 - time_layer_first)*npls+ipls)*d->n_pointresB + i;

				id[0]=((j - time_layer_first)*npls+ipls)*d->n_pointresB + i;
				id[1]=((j_plus_1 - time_layer_first)*npls+ipls)*d->n_pointresB + i;

				F[0][2]=Bx_all[id[2]];

				F[0][0]=Bx_all[id[0]];
				F[0][1]=Bx_all[id[1]];

				F[1][2]=By_all[id[2]];

				F[1][0]=By_all[id[0]];
				F[1][1]=By_all[id[1]];

				F[2][2]=Bz_all[id[2]];

				F[2][0]=Bz_all[id[0]];
				F[2][1]=Bz_all[id[1]];

				dt=d->time[j_plus_2]-d->time[j];
				dt1=d->time[j_plus_1]-d->time[j];
				dt0=d->time[j_plus_2]-d->time[j_plus_1];

			}while(true);

			mt0=(dt0)/(dt1*dt);
			mt1=-(dt)/(dt1*dt0);
			mt2=(dt+dt0)/(dt*dt0);

			Edsx_all[id[2]]=(F[0][0]*mt0+F[0][1]*mt1+F[0][2]*mt2);
			Edsy_all[id[2]]=(F[1][0]*mt0+F[1][1]*mt1+F[1][2]*mt2);
			Edsz_all[id[2]]=(F[2][0]*mt0+F[2][1]*mt1+F[2][2]*mt2);
		}

		for(i=0;i<d->n_pointresE;i++)
		{
			if(i<RecvPlsIgE[ipls] || i>=RecvPlsIgE[ipls+1]){continue;}

			j=time_layer_first;
			id[0]=((j - time_layer_first)*npls+ipls)*d->n_pointresE + i;

			j_plus_1=j+1;
			id[1]=((j_plus_1 - time_layer_first)*npls+ipls)*d->n_pointresE + i;

			j_plus_2=j_plus_1+1;
			id[2]=((j_plus_2 - time_layer_first)*npls+ipls)*d->n_pointresE + i;


			

			F[0][0]=Ax_all[id[0]];
			F[0][1]=Ax_all[id[1]];
			F[0][2]=Ax_all[id[2]];

			F[1][0]=Ay_all[id[0]];
			F[1][1]=Ay_all[id[1]];
			F[1][2]=Ay_all[id[2]];
			
			F[2][0]=Az_all[id[0]];
			F[2][1]=Az_all[id[1]];
			F[2][2]=Az_all[id[2]];

			dt=d->time[j_plus_2]-d->time[j];
			dt1=d->time[j_plus_1]-d->time[j];
			dt0=d->time[j_plus_2]-d->time[j_plus_1];

			mt0=-(dt+dt1)/(dt1*dt);
			mt1=(dt)/(dt1*dt0);
			mt2=(-dt1)/(dt*dt0);

			Ex_all[id[0]]=-(F[0][0]*mt0+F[0][1]*mt1+F[0][2]*mt2);
			Ey_all[id[0]]=-(F[1][0]*mt0+F[1][1]*mt1+F[1][2]*mt2);
			Ez_all[id[0]]=-(F[2][0]*mt0+F[2][1]*mt1+F[2][2]*mt2);
			
			do{
				mt0=(-dt0)/(dt1*dt);
				mt1=-(dt1-dt0)/(dt1*dt0);
				mt2=(dt1)/(dt*dt0);
			
				Ex_all[id[1]]=-(F[0][0]*mt0+F[0][1]*mt1+F[0][2]*mt2);
				Ey_all[id[1]]=-(F[1][0]*mt0+F[1][1]*mt1+F[1][2]*mt2);
				Ez_all[id[1]]=-(F[2][0]*mt0+F[2][1]*mt1+F[2][2]*mt2);



				if(j_plus_2==time_layer_last)break;

				j=j_plus_1;
				j_plus_2++;

				j_plus_1=vPre1[j_plus_2];

				id[2]=((j_plus_2 - time_layer_first)*npls+ipls)*d->n_pointresE + i;

				id[0]=((j - time_layer_first)*npls+ipls)*d->n_pointresE + i;
				id[1]=((j_plus_1 - time_layer_first)*npls+ipls)*d->n_pointresE + i;

				F[0][2]=Ax_all[id[2]];

				F[0][0]=Ax_all[id[0]];
				F[0][1]=Ax_all[id[1]];

				F[1][2]=Ay_all[id[2]];

				F[1][0]=Ay_all[id[0]];
				F[1][1]=Ay_all[id[1]];


				F[2][2]=Az_all[id[2]];

				F[2][0]=Az_all[id[0]];
				F[2][1]=Az_all[id[1]];

				dt=d->time[j_plus_2]-d->time[j];
				dt1=d->time[j_plus_1]-d->time[j];
				dt0=d->time[j_plus_2]-d->time[j_plus_1];
			}while(true);

			mt0=(dt0)/(dt1*dt);
			mt1=-(dt)/(dt1*dt0);
			mt2=(dt+dt0)/(dt*dt0);

			Ex_all[id[2]]=-(F[0][0]*mt0+F[0][1]*mt1+F[0][2]*mt2);
			Ey_all[id[2]]=-(F[1][0]*mt0+F[1][1]*mt1+F[1][2]*mt2);
			Ez_all[id[2]]=-(F[2][0]*mt0+F[2][1]*mt1+F[2][2]*mt2);
		}




























































		sprintf(fname,"edsall_anom.%d",ipls+1);
		outf.open(fname);
		outf<<scientific<<setprecision(16);
		for (i=0; i<d->n_pointresB; i++)
		{
			if(i<RecvPlsIgB[ipls] || i>=RecvPlsIgB[ipls+1]){continue;}

			outf<<(int)d->pointresB[i][0]<<' '<<(int)d->pointresB[i][1]<<endl;		
			outf<<d->pointresB[i][0]<<' '<<d->pointresB[i][1]<<' '<<d->pointresB[i][2]<<endl;		
			outf<<"1\n1 1 1 1\n1\n";
			outf<<"      t (ms)         Emfx_anom (mV)     Emfy_anom (mV)     Emfz_anom (mV)\n";

			for(j=time_layer_first;j<=time_layer_last;j++)
			{

				ind = ((j - time_layer_first)*npls+ipls)*d->n_pointresB + i;	
				outf<<d->time[j]*1000<<"\t";
				outf<<sqrloop*Edsx_all[ind]*(-1000)<<"\t";
				outf<<sqrloop*Edsy_all[ind]*(-1000)<<"\t";
				outf<<sqrloop*Edsz_all[ind]*(-1000)<<"\n";
			}
		}
		if(d->n_pointresB==0){outf<<' '<<endl;}
		outf.close();
		outf.clear();

		sprintf(fname,"ball_anom.%d",ipls+1);
		outf.open(fname);
		outf<<scientific<<setprecision(16);
		for(i=0;i<d->n_pointresB;i++)
		{
			if(i<RecvPlsIgB[ipls] || i>=RecvPlsIgB[ipls+1]){continue;}

			outf<<(int)d->pointresB[i][0]<<' '<<(int)d->pointresB[i][1]<<endl;		
			outf<<d->pointresB[i][0]<<' '<<d->pointresB[i][1]<<' '<<d->pointresB[i][2]<<endl;		
			outf<<"1\n1 1 1 1\n1\n";
			outf<<"      t (ms)        Bx_anom (mTl)     By_anom (mTl)     Bz_anom (mTl)\n";

			for(j=time_layer_first;j<=time_layer_last;j++)
			{

				ind = ((j - time_layer_first)*npls+ipls)*d->n_pointresB + i;	
				outf<<d->time[j]*1000<<"\t";
				outf<<Bx_all[ind]*1000<<"\t";
				outf<<By_all[ind]*1000<<"\t";
				outf<<Bz_all[ind]*1000<<"\n";
			}
		}
		if(d->n_pointresB==0){outf<<' '<<endl;}
		outf.close();
		outf.clear();


		sprintf(fname,"eall_anom.%d",ipls+1);
		outf.open(fname);
		outf<<scientific<<setprecision(16);
		for(i=0;i<d->n_pointresE;i++)
		{
			if(i<RecvPlsIgE[ipls] || i>=RecvPlsIgE[ipls+1]){continue;}

			outf<<(int)d->pointresE[i][0]<<' '<<(int)d->pointresE[i][1]<<endl;		
			outf<<d->pointresE[i][0]<<' '<<d->pointresE[i][1]<<' '<<d->pointresE[i][2]<<endl;		
			outf<<"1\n1 1 1 1\n1\n";
			outf<<"      t (ms)         Ex_anom n (mV)     Ey_anom a (mV)     Ez_anom s (mV)\n";

			for(j=time_layer_first;j<=time_layer_last;j++)
			{

				ind = ((j - time_layer_first)*npls+ipls)*d->n_pointresE + i;
				outf<<d->time[j]*1000<<"\t";
				outf<<Ex_all[ind]*1000<<"\t";
				outf<<Ey_all[ind]*1000<<"\t";
				outf<<Ez_all[ind]*1000<<"\n";
			}
		}
		if(d->n_pointresE==0){outf<<' '<<endl;}
		outf.close();
		outf.clear();
	}

	if(Ax_all){delete [] Ax_all; Ax_all=NULL;}
	if(Ay_all){delete [] Ay_all; Ay_all=NULL;}
	if(Az_all){delete [] Az_all; Az_all=NULL;}
}
void Give_out_vec_loop::Work_on_current_time_layer()
{
	int i, j, e, num_hex;
	double Bx, By, Bz;
	
	double x[8], y[8], z[8];
	
	double ves_j[12], ves_j1[12], ves_j2[12];
	double val_j[3], val_j1[3], val_j2[3];

	double local_coords[3];


	dt = t_j - t_j2; 
	dt0 = t_j - t_j1;
	dt1 = t_j1 - t_j2;


	if(d->n_pointresB)
	{
		if(CheckStop())return;
		resultantB->ValueType=vtRotxA;
		resultantB->Output(tnum,ipls_cur,cur_dec_strt_res,0,vvtb);
		if(CheckStop())return;
		resultantB->ValueType=vtRotyA;
		resultantB->Output(tnum,ipls_cur,cur_dec_strt_res,1,vvtb);
		if(CheckStop())return;
		resultantB->ValueType=vtRotzA;
		resultantB->Output(tnum,ipls_cur,cur_dec_strt_res,2,vvtb);
	}

	if(d->n_pointresE)
	{
		if(CheckStop())return;
		resultantA->ValueType=vtAx;
		resultantA->Output(tnum,ipls_cur,cur_dec_strt_res,0,vvta);
		if(CheckStop())return;
		resultantA->ValueType=vtAy;
		resultantA->Output(tnum,ipls_cur,cur_dec_strt_res,1,vvta);
		if(CheckStop())return;
		resultantA->ValueType=vtAz;
		resultantA->Output(tnum,ipls_cur,cur_dec_strt_res,2,vvta);
	}

}
double Give_out_vec_loop::dA_dt(double t, double u_j, double u_j1, double u_j2,
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
