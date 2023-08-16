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
 *  This file contains the code for reading and storing edge mesh in 3D VFEM
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin                         
 *  Novosibirsk State Technical University,                 
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia       
 *  vdv_wk@mail.ru                                          
 *  Version 2.0 January 16, 2023                            
*/

#include "stdafx.h"
#include "in_out.h"
#include "vec_prep_data.h"
#include "vfem_const.h"


extern void Memory_allocation_error(const char *var, const char *func);
extern void Cannot_open_file(const char *fname, const char *func);

extern ofstream logfile;

extern bool IsFileExist(char *fname);

Vec_Prep_Data::Vec_Prep_Data()
{
	nver = NULL;
	nvkat = NULL;
	xyz = NULL;
	mu3d = NULL;
	mu0 = NULL;
    sigma3d = NULL;
	sigma0 = NULL;
	n_pointresB=0;
	pointresB = NULL;
	n_pointresE=0;
	pointresE = NULL;
	mesh_regular_x = NULL;
	mesh_regular_y = NULL;
	l13d = NULL;
	time = NULL;
	tasktype = 0;
	dpr3d = NULL;
	dpr0 = NULL;
	nobj=0;
}
Vec_Prep_Data::~Vec_Prep_Data()
{
	if(mu0) {delete [] mu0; mu0=NULL;}
	if(xyz) {delete [] xyz; xyz=NULL;}
	if(nver) {delete [] nver; nver=NULL;}
	if(mu3d) {delete [] mu3d; mu3d=NULL;}
	if(l13d) {delete [] l13d; l13d=NULL;}
	if(time) {delete [] time; time=NULL;}
	if(nvkat) {delete [] nvkat; nvkat=NULL;}
	if(sigma0) {delete [] sigma0; sigma0=NULL;}
	if(sigma3d) {delete [] sigma3d; sigma3d=NULL;}
	if(pointresB) {delete [] pointresB; pointresB=NULL;}
	if(pointresE) {delete [] pointresE; pointresE=NULL;}
	if(mesh_regular_x) {delete [] mesh_regular_x; mesh_regular_x=NULL;}
	if(mesh_regular_y) {delete [] mesh_regular_y; mesh_regular_y=NULL;}
	if(dpr3d) {delete [] dpr3d; dpr3d=NULL;}
	if(dpr0) {delete [] dpr0; dpr0=NULL;}
}
int Vec_Prep_Data::Read_mesh_for_nonstat_problem(char *pointres_fname)
{
	In_Out R;
	FILE *fp=NULL;
	double temp1, temp2;
	int i, j;
	int max_material;
	int tmp;

	int n_reg;
	ifstream inf;

	inf.open("sig3d");
	if(!inf)
		Cannot_open_file("sig3d", "Vec_Prep_Data::Read_mesh_for_nonstat_problem");

	max_material = 0;
	while(!inf.eof())
	{
		inf>>tmp>>temp1>>temp2;
		if (!inf.good())
			break;
		if(tmp > max_material)
			max_material = tmp;
	}
	n_materials = max_material;
	inf.close();
	inf.clear();

	mu3d = new double[n_materials];
	if(mu3d == 0)
		Memory_allocation_error("mu3d", "Vec_Prep_Data::Read_mesh_for_nonstat_problem");

	mu0 = new double[n_materials];
	if(mu0 == 0)
		Memory_allocation_error("mu0", "Vec_Prep_Data::Read_mesh_for_nonstat_problem");

	inf.open("mu3d");
	if(!inf)
		Cannot_open_file("mu3d", "Vec_Prep_Data::Read_prep_data");

	for (i=0; i<n_materials; i++)
	{
		int n_of_current_material;
		inf>>n_of_current_material;
		n_of_current_material--;
		inf>>mu3d[n_of_current_material]>>mu0[n_of_current_material];
		mu3d[n_of_current_material]*=MU_0;
		mu0[n_of_current_material]*=MU_0;
	}
	inf.close();
	inf.clear();


	sigma3d = new double[n_materials];
	if(sigma3d == 0)
		Memory_allocation_error("sigma3d", "Vec_Prep_Data::Read_mesh_for_nonstat_problem");

	sigma0 = new double[n_materials];
	if(sigma0 == 0)
		Memory_allocation_error("sigma0", "Vec_Prep_Data::Read_mesh_for_nonstat_problem");

	inf.open("sig3d");
	if(!inf)
		Cannot_open_file("sig3d", "Vec_Prep_Data::Read_mesh_for_nonstat_problem");

	for (i=0; i<n_materials; i++)
	{
		int n_of_current_material;
		inf>>n_of_current_material;
		n_of_current_material--;
		inf>>sigma3d[n_of_current_material]>>sigma0[n_of_current_material];
	}
	inf.close();
	inf.clear();

	dpr3d = new double[n_materials];
	if(dpr3d == 0)
		Memory_allocation_error("dpr3d", "Vec_Prep_Data::Read_prep_data");

	dpr0 = new double[n_materials];
	if(dpr0 == 0)
		Memory_allocation_error("dpr0", "Vec_Prep_Data::Read_prep_data");

	inf.open("mu3d");
	if(!inf)
		Cannot_open_file("dpr3d", "Vec_Prep_Data::Read_prep_data");

	for (i=0; i<n_materials; i++)
	{
		int n_of_current_material;
		inf>>n_of_current_material;
		n_of_current_material--;
		inf>>dpr3d[n_of_current_material]>>dpr0[n_of_current_material];
		dpr3d[n_of_current_material]*=DPR_0;
		dpr0[n_of_current_material]*=DPR_0;
	}
	inf.close();
	inf.clear();

	if(strlen(pointres_fname)!=0)
	{
		if((fp=fopen("xyzVectorB", "r"))==0)
			Cannot_open_file("xyzVectorB", "Vec_Prep_Data::Read_mesh_for_nonstat_problem");

		fscanf(fp, "%ld", &n_pointresB);

		pointresB = new double[n_pointresB][3];
		if(pointresB == 0)
			Memory_allocation_error("pointresB", "Vec_Prep_Data::Read_mesh_for_nonstat_problem");

		for(i=0; i<n_pointresB; i++)
			fscanf(fp, "%lf %lf %lf", &pointresB[i][0], &pointresB[i][1], &pointresB[i][2]);
		fclose(fp);

		if((fp=fopen("xyzVectorE", "r"))==0)
			Cannot_open_file("xyzVectorE", "Vec_Prep_Data::Read_mesh_for_nonstat_problem");

		fscanf(fp, "%ld", &n_pointresE);

		pointresE = new double[n_pointresE][3];
		if(pointresE == 0)
			Memory_allocation_error("pointresE", "Vec_Prep_Data::Read_mesh_for_nonstat_problem");

		for(i=0; i<n_pointresE; i++)
			fscanf(fp, "%lf %lf %lf", &pointresE[i][0], &pointresE[i][1], &pointresE[i][2]);
		fclose(fp);
	}

	R.Read_inftry("inftry.dat", &kuzlov, &kpar, &kt1);

	nver = new int[kpar][14];
	if(nver == 0)
		Memory_allocation_error("nver", "Vec_Prep_Data::Read_mesh_for_nonstat_problem");

	R.Read_Bin_File_Of_Long("nver.dat", (int*)nver, kpar, 14);
	for(i=0; i<kpar; i++)
		for(j=0; j<14; j++)
			nver[i][j]--;

	nvkat = new int[kpar];
	if(nvkat == 0)
		Memory_allocation_error("nvkat", "Vec_Prep_Data::Read_mesh_for_nonstat_problem");

	R.Read_Bin_File_Of_Long("nvkat.dat", nvkat, kpar, 1);
	for(i=0; i<kpar; i++)
		nvkat[i]--;

	xyz = new double[kuzlov][3];
	if(xyz == 0)
		Memory_allocation_error("xyz", "Vec_Prep_Data::Read_mesh_for_nonstat_problem");

	R.Read_Bin_File_Of_Double("xyz.dat", (double*)xyz, kuzlov, 3);

	l13d = new int[kt1];
	if(l13d == 0)
		Memory_allocation_error("l13d", "Vec_Prep_Data::Read_mesh_for_nonstat_problem");

	R.Read_Bin_File_Of_Long("l13d.dat", l13d, kt1, 1);
	for(i=0; i<kt1; i++)
		l13d[i]--;


	return 0;
}
int Vec_Prep_Data::Read_3dmeshregular(int interval)
{
	int i, j, t, r;
	int temp;
	int int_whole, int_ost;
	int nx, ny, nz;
	ifstream inf;

	inf.open("3dmeshregular");
	if (!inf)
		Cannot_open_file("3dmeshregular", "Vec_Prep_Data::Read_3dmeshregular");

	inf>>nx;

	N_X=nx;
	Xcrd.resize(N_X);


	int_whole = nx/(1+interval);
	if(int_whole*(1+interval)-nx==0)
	{
		int_ost = 0;
	}
	else
	{	
		int_ost = 1;	
	}
	n_mesh_regular_x = int_whole + int_ost;

	mesh_regular_x = new double[n_mesh_regular_x];
	if(mesh_regular_x == 0)
		Memory_allocation_error("mesh_regular_x", "Vec_Prep_Data::Read_3dmeshregular");

	r=0;
	for(i=0; i<int_whole; i++)
	{
		t = i*(interval+1);
		inf>>temp>>Xcrd[r];
		mesh_regular_x[i]=Xcrd[r];
		r++;
		for(j=0; j<interval; j++)
		{
			inf>>temp>>Xcrd[r];
			r++;
			t++;
		}
	}

	if(int_ost!=0)
	{
		for(i=t+1; i<nx-1; i++){
			inf>>temp>>Xcrd[r];
			r++;
		}
		inf>>temp>>Xcrd[r];
		mesh_regular_x[n_mesh_regular_x-1]=Xcrd[r];
		r++;
	}

	inf>>ny;

	N_Y=ny;
	Ycrd.resize(N_Y);

	int_whole = ny/(1+interval);
	if(int_whole*(1+interval)-ny==0)
	{
		int_ost = 0;
	}
	else
	{	
		int_ost = 1;	
	}
	n_mesh_regular_y = int_whole + int_ost;

	mesh_regular_y = new double[n_mesh_regular_y];
	if(mesh_regular_y == 0)
		Memory_allocation_error("mesh_regular_y", "Vec_Prep_Data::Read_3dmeshregular");

	r=0;
	for(i=0; i<int_whole; i++)
	{
		t = i*(interval+1);
		inf>>temp>>Ycrd[r];
		mesh_regular_y[i]=Ycrd[r];
		r++;
		for(j=0; j<interval; j++)
		{
			inf>>temp>>Ycrd[r];
			r++;
			t++;
		}
	}

	if(int_ost!=0)
	{
		for(i=t+1; i<ny-1; i++){
			inf>>temp>>Ycrd[r];
			r++;
		}
		inf>>temp>>Ycrd[r];
		mesh_regular_y[n_mesh_regular_y-1]=Ycrd[r];
		r++;
	}

	inf>>nz;

	N_Z=nz;
	Zcrd.resize(N_Z);
	for(r=0;r<N_Z;r++){
		inf>>temp>>Zcrd[r];
	}

	inf>>nobj;
	if(inf.good())
	{
		Xobj[0].resize(nobj);
		Xobj[1].resize(nobj);
		Yobj[0].resize(nobj);
		Yobj[1].resize(nobj);
		Zobj[0].resize(nobj);
		Zobj[1].resize(nobj);
		Mobj.resize(nobj);
		for(i=0;i<nobj;i++)
		{
			inf>>temp;Xobj[0][i]=Xcrd[temp-1];
			inf>>temp;Xobj[1][i]=Xcrd[temp-1];
			inf>>temp;Yobj[0][i]=Ycrd[temp-1];
			inf>>temp;Yobj[1][i]=Ycrd[temp-1];
			inf>>temp;Zobj[0][i]=Zcrd[temp-1];
			inf>>temp;Zobj[1][i]=Zcrd[temp-1];
			inf>>Mobj[i];
		}
		inf.close();
	}
	else
	{
		nobj=0;
	}
	inf.clear();

	return 0;
}
int Vec_Prep_Data::Read_infite0()
{
	ifstream inf("infite.0");
	int i;
	inf.ignore(1000, '=');
	inf >> ntime;
	time=new double[ntime];
	inf.ignore(1000, '\n');
	inf.ignore(1000, '\n');
	inf.ignore(1000, '\n');
	inf.ignore(1000, '\n');
	for (i=0; i<ntime; i++)
	{
		if (!inf.good()) return 1;
		inf >> time[i];
		inf.ignore(1000, ';');
	}
	inf.close();
	return 0;
}
int Vec_Prep_Data::LoadReceivers(char *fname, int& _n_pointres, double (*(&_pointres))[3])
{
	FILE *fp=NULL;
	int i;

	if (_pointres) {delete [] _pointres; _pointres=NULL;}
	_n_pointres = 0;

	if((fp=fopen(fname, "r"))==0)
		Cannot_open_file(fname, "Vec_Prep_Data::LoadReceivers");

	fscanf(fp, "%ld", &_n_pointres);

	_pointres = new double[_n_pointres][3];
	if(_pointres == 0)
		Memory_allocation_error("_pointres", "Vec_Prep_Data::LoadReceivers");

	for(i=0; i<_n_pointres; i++)
		fscanf(fp, "%lf %lf %lf", &_pointres[i][0], &_pointres[i][1], &_pointres[i][2]);
	fclose(fp);

	return 0;
}
