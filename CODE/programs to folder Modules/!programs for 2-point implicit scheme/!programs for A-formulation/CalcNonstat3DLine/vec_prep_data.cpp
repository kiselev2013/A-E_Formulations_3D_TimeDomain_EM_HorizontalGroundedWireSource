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
	usin = NULL;
	ucos = NULL;
	z_1d = NULL;
	sigma_1d = NULL;
	layers_1d = NULL;
	time = NULL;
	tasktype = 0;
	dpr3d = NULL;
	dpr0 = NULL;
	fdirect=false;
	xyzt=NULL;
	nSrsRazb=0;
}
Vec_Prep_Data::~Vec_Prep_Data()
{
	if(mu0) {delete [] mu0; mu0=NULL;}
	if(xyz) {delete [] xyz; xyz=NULL;}
	if(nver) {delete [] nver; nver=NULL;}
	if(mu3d) {delete [] mu3d; mu3d=NULL;}
	if(l13d) {delete [] l13d; l13d=NULL;}
	if(usin) {delete [] usin; usin=NULL;}
	if(ucos) {delete [] ucos; ucos=NULL;}
	if(z_1d) {delete [] z_1d; z_1d=NULL;}
	if(time) {delete [] time; time=NULL;}
	if(nvkat) {delete [] nvkat; nvkat=NULL;}
	if(sigma0) {delete [] sigma0; sigma0=NULL;}
	if(sigma3d) {delete [] sigma3d; sigma3d=NULL;}
	if(pointresB) {delete [] pointresB; pointresB=NULL;}
	if(pointresE) {delete [] pointresE; pointresE=NULL;}
	if(sigma_1d) {delete [] sigma_1d; sigma_1d=NULL;}
	if(layers_1d) {delete [] layers_1d; layers_1d=NULL;}
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


	fdirect=false;
	inf.open("fdirect");
	if(inf)
	{
		inf>>fdirect;
		inf.close();
	}
	inf.clear();


	return 0;
}
int Vec_Prep_Data::ReadPrepDataHarmLoop(char *pointres_fname)
{
	In_Out R;
	FILE *fp=NULL;
	double temp1, temp2;
	int i, j;
	int max_material;
	int tmp;

	In_Out r;

	int n_reg;
	ifstream inf;


	r.Read_Double_From_Txt_File("nu", &nu);


	if((fp=fopen("sig3d", "r"))==0)
		Cannot_open_file("sig3d", "Vec_Prep_Data::Read_mesh_for_nonstat_problem");

	max_material = 0;
	while(!feof(fp))
	{
		fscanf(fp, "%ld %lf %lf", &tmp, &temp1, &temp2);
		if(tmp > max_material)
			max_material = tmp;
	}
	n_materials = max_material;
	fclose(fp);

	mu3d = new double[n_materials];
	if(mu3d == 0)
		Memory_allocation_error("mu3d", "Vec_Prep_Data::Read_mesh_for_nonstat_problem");

	mu0 = new double[n_materials];
	if(mu0 == 0)
		Memory_allocation_error("mu0", "Vec_Prep_Data::Read_mesh_for_nonstat_problem");

	if((fp=fopen("mu3d", "r"))==0)
		Cannot_open_file("mu3d", "Vec_Prep_Data::Read_prep_data");

	for (i=0; i<n_materials; i++)
	{
		int n_of_current_material;
		fscanf(fp, "%ld", &n_of_current_material);
		n_of_current_material--;
		fscanf(fp, "%lf %lf",
			&mu3d[n_of_current_material], &mu0[n_of_current_material]);
		mu3d[n_of_current_material]*=MU_0;
		mu0[n_of_current_material]*=MU_0;
	}
	fclose(fp);


	sigma3d = new double[n_materials];
	if(sigma3d == 0)
		Memory_allocation_error("sigma3d", "Vec_Prep_Data::Read_mesh_for_nonstat_problem");

	sigma0 = new double[n_materials];
	if(sigma0 == 0)
		Memory_allocation_error("sigma0", "Vec_Prep_Data::Read_mesh_for_nonstat_problem");

	if((fp=fopen("sig3d", "r"))==0)
		Cannot_open_file("sig3d", "Vec_Prep_Data::Read_mesh_for_nonstat_problem");

	for (i=0; i<n_materials; i++)
	{
		int n_of_current_material;
		fscanf(fp, "%ld", &n_of_current_material);
		n_of_current_material--;
		fscanf(fp, "%lf %lf",
			&sigma3d[n_of_current_material], &sigma0[n_of_current_material]);
	}
	fclose(fp);

	dpr3d = new double[n_materials];
	if(dpr3d == 0)
		Memory_allocation_error("dpr3d", "Vec_Prep_Data::Read_prep_data");

	dpr0 = new double[n_materials];
	if(dpr0 == 0)
		Memory_allocation_error("dpr0", "Vec_Prep_Data::Read_prep_data");

	if((fp=fopen("dpr3d", "r"))==0)
		Cannot_open_file("dpr3d", "Vec_Prep_Data::Read_prep_data");

	for (i=0; i<n_materials; i++)
	{
		int n_of_current_material;
		fscanf(fp, "%ld", &n_of_current_material);
		n_of_current_material--;
		fscanf(fp, "%lf %lf",
			&dpr3d[n_of_current_material], &dpr0[n_of_current_material]);
		dpr3d[n_of_current_material]*=DPR_0;
		dpr0[n_of_current_material]*=DPR_0;
	}
	fclose(fp);

	if(strlen(pointres_fname)!=0)
	{
		if((fp=fopen(pointres_fname, "r"))==0)
			Cannot_open_file(pointres_fname, "Vec_Prep_Data::Read_mesh_for_nonstat_problem");

		fscanf(fp, "%ld", &n_pointresB);

		pointresB = new double[n_pointresB][3];
		if(pointresB == 0)
			Memory_allocation_error("pointresB", "Vec_Prep_Data::Read_mesh_for_nonstat_problem");

		for(i=0; i<n_pointresB; i++)
			fscanf(fp, "%lf %lf %lf", &pointresB[i][0], &pointresB[i][1], &pointresB[i][2]);
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
int Vec_Prep_Data::Read_prep_data()
{
	In_Out R;
	FILE *fp=NULL;
	double temp, temp1, temp2;
	int i, j;
	int max_material;
	int tmp;

	int n_reg;
	ifstream inf;


	if (tasktype==0) //  
	{
		if((fp=fopen("config", "r"))==0)
			Cannot_open_file("config", "Vec_Prep_Data::Read_prep_data");

		fscanf(fp, "%lf", &temp);
		fscanf(fp, "%ld", &maxiter);
		fscanf(fp, "%lf", &eps);
		fclose(fp);
	}
	else
	{
		maxiter=10000;
		eps=1e-6;
	}

	if((fp=fopen("sig3d", "r"))==0)
		Cannot_open_file("sig3d", "Vec_Prep_Data::Read_prep_data");

	max_material = 0;
	while(!feof(fp))
	{
		fscanf(fp, "%ld %lf %lf", &tmp, &temp1, &temp2);
		if(tmp > max_material)
			max_material = tmp;
	}
	n_materials = max_material;
	fclose(fp);

	mu3d = new double[n_materials];
	if(mu3d == 0)
		Memory_allocation_error("mu3d", "Vec_Prep_Data::Read_prep_data");

	mu0 = new double[n_materials];
	if(mu0 == 0)
		Memory_allocation_error("mu0", "Vec_Prep_Data::Read_prep_data");

	if((fp=fopen("mu3d", "r"))==0)
		Cannot_open_file("mu3d", "Vec_Prep_Data::Read_prep_data");

	for (i=0; i<n_materials; i++)
	{
		int n_of_current_material;
		fscanf(fp, "%ld", &n_of_current_material);
		n_of_current_material--;
		fscanf(fp, "%lf %lf",
			&mu3d[n_of_current_material], &mu0[n_of_current_material]);
		mu3d[n_of_current_material]*=MU_0;
		mu0[n_of_current_material]*=MU_0;
	}
	fclose(fp);


	sigma3d = new double[n_materials];
	if(sigma3d == 0)
		Memory_allocation_error("sigma3d", "Vec_Prep_Data::Read_prep_data");

	sigma0 = new double[n_materials];
	if(sigma0 == 0)
		Memory_allocation_error("sigma0", "Vec_Prep_Data::Read_prep_data");	

	if((fp=fopen("sig3d", "r"))==0)
		Cannot_open_file("sig3d", "Vec_Prep_Data::Read_prep_data");

	for (i=0; i<n_materials; i++)
	{
		int n_of_current_material;
		fscanf(fp, "%ld", &n_of_current_material);
		n_of_current_material--;
		fscanf(fp, "%lf %lf",
			&sigma3d[n_of_current_material], &sigma0[n_of_current_material]);
	}
	fclose(fp);

	dpr3d = new double[n_materials];
	if(dpr3d == 0)
		Memory_allocation_error("dpr3d", "Vec_Prep_Data::Read_prep_data");

	dpr0 = new double[n_materials];
	if(dpr0 == 0)
		Memory_allocation_error("dpr0", "Vec_Prep_Data::Read_prep_data");

	if((fp=fopen("dpr3d", "r"))==0)
		Cannot_open_file("dpr3d", "Vec_Prep_Data::Read_prep_data");

	for (i=0; i<n_materials; i++)
	{
		int n_of_current_material;
		fscanf(fp, "%ld", &n_of_current_material);
		n_of_current_material--;
		fscanf(fp, "%lf %lf",
			&dpr3d[n_of_current_material], &dpr0[n_of_current_material]);
		dpr3d[n_of_current_material]*=DPR_0;
		dpr0[n_of_current_material]*=DPR_0;
	}
	fclose(fp);


	if (tasktype==0)
	{
		LoadReceivers("pointres", n_pointresB, pointresB);
		LoadReceivers("pointres", n_pointresE, pointresE);
	}






	R.Read_inftry("inftry.dat", &kuzlov, &kpar, &kt1);

	nver = new int[kpar][14];
	if(nver == 0)
		Memory_allocation_error("nver", "Vec_Prep_Data::Read_prep_data");

	R.Read_Bin_File_Of_Long("nver.dat", (int*)nver, kpar, 14);
	for(i=0; i<kpar; i++)
		for(j=0; j<14; j++)
			nver[i][j]--;

	nvkat = new int[kpar];
	if(nvkat == 0)
		Memory_allocation_error("nvkat", "Vec_Prep_Data::Read_prep_data");

	R.Read_Bin_File_Of_Long("nvkat.dat", nvkat, kpar, 1);
	for(i=0; i<kpar; i++)
		nvkat[i]--;

	xyz = new double[kuzlov][3];
	if(xyz == 0)
		Memory_allocation_error("xyz", "Vec_Prep_Data::Read_prep_data");

	R.Read_Bin_File_Of_Double("xyz.dat", (double*)xyz, kuzlov, 3);

	l13d = new int[kt1];
	if(l13d == 0)
		Memory_allocation_error("l13d", "Vec_Prep_Data::Read_prep_data");

	R.Read_Bin_File_Of_Long("l13d.dat", l13d, kt1, 1);
	for(i=0; i<kt1; i++)
		l13d[i]--;


	if (tasktype==0) //  
	{
		if((fp=fopen("sreda1d.ay","r"))==0)
			Cannot_open_file("sreda1d.ay", "Vec_Prep_Data::Read_prep_data");

		fscanf(fp, "%ld", &n_layers_1d);

		layers_1d = new double[n_layers_1d];
		if(layers_1d == 0)
			Memory_allocation_error("layers_1d", "Vec_Prep_Data::Read_prep_data");

		sigma_1d = new double[n_layers_1d];
		if(sigma_1d == 0)
			Memory_allocation_error("sigma_1d", "Vec_Prep_Data::Read_prep_data");

		for(i=0; i<n_layers_1d; i++)
		{
			fscanf(fp, "%lf", &layers_1d[i]);
			fscanf(fp, "%lf", &sigma_1d[i]);
			fscanf(fp, "%lf", &temp); // mu 
		}

		fclose(fp);
	}
	else
	{
		In_Out r;
		
		r.Read_Double_From_Txt_File("nu", &nu);

		r.Read_Long_From_Txt_File("norvect", &norvect);
	}

	return 0;
}
int Vec_Prep_Data::UnloadxyzVectorE0ForLine(const int &n, int (*ed)[25], int (*edges)[2])
{
	int i, j, k;
	bool f;
	const double eps = 1e-6;
	EnLine En;
	
	EnForLine.resize(n);
	
	for (i=0; i<kpar; i++)
	{
		if(fabs(sigma0[nvkat[i]] - sigma3d[nvkat[i]]) > eps)
		{
			for (j=0; j<12; j++)
			{
				vector<EnLine>& Enij=EnForLine[ed[i][j]];
				f=false;
				for(k=0; (k<(int)Enij.size())&&!f; k++)
					f=(Enij[k].mtr==nvkat[i]);
				if (!f)
				{
					En.mtr=nvkat[i];
					Enij.push_back(En);
				}
			}
		}
	}

	int xyzVectorE0_n=0;
	for (i=0; i<n; i++)
		xyzVectorE0_n+=(int)EnForLine[i].size();
		
	ofstream outf1, outf2;

	outf1.open("xyzVectorE0");
	outf2.open("xyzVectorE0n");

	outf1<<xyzVectorE0_n<<endl;
	outf2<<xyzVectorE0_n<<endl;

	for (i=0; i<n; i++)
	{
		double pntEn[3];
		for (j=0; j<3; j++)
			pntEn[j] = 0.5*(xyz[edges[i][0]][j] + xyz[edges[i][1]][j]);
		for (j=0; j<(int)EnForLine[i].size(); j++)
		{
			outf1<<pntEn[0]<<' '<<pntEn[1]<<' '<<pntEn[2]<<endl;
			outf2<<EnForLine[i][j].mtr+1<<endl;
		}
	}

	outf1.close();
	outf2.close();
	
	return 0;
}
int Vec_Prep_Data::LoadVectorE0ForLine()
{
	int i, j;
	ifstream inf1, inf2;
	inf1.open("e_s.dat", ios::binary);
	inf2.open("e_c.dat", ios::binary);
	if (!inf1)
		return 1;
	if (!inf2)
	{
		inf1.close();
		return 1;
	}
	for (i=0; i<(int)EnForLine.size(); i++)
		for (j=0; j<(int)EnForLine[i].size(); j++)
		{
			EnLine& EnL=EnForLine[i][j];
			inf1>EnL.es[0]>EnL.es[1]>EnL.es[2];
			inf2>EnL.ec[0]>EnL.ec[1]>EnL.ec[2];
		}
	inf1.close();
	inf2.close();
	return 0;
}
int Vec_Prep_Data::Read_mtz_1d()
{
	In_Out r;
	double alpha_double;

	r.Read_Double_From_Txt_File("alfa", &alpha_double);
	if (fabs(1.0 - alpha_double)<0.01)
	{
		alfa = 1; 
	} 
	else
	{
		alfa = 0;
	}

	r.Read_Double_From_Txt_File("nu", &nu);

	r.Read_Long_From_Txt_File("norvect", &norvect);

	r.Read_Long_From_Txt_File("setka1DEy", &n_1d);

	z_1d = new double[n_1d];
	if(z_1d == 0)
		Memory_allocation_error("z_1d", "Vec_Prep_Data::Read_mtz_1d");

	usin = new double[n_1d];
	if(usin == 0)
		Memory_allocation_error("usin", "Vec_Prep_Data::Read_mtz_1d");

	ucos = new double[n_1d];
	if(ucos == 0)
		Memory_allocation_error("ucos", "Vec_Prep_Data::Read_mtz_1d");

	return r.Read_1d_data(n_1d, z_1d, usin, ucos);
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

	inf.close();
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

int Vec_Prep_Data::GetNearestElement(loc_source &ls)
{
	double scurr,sbest,m_s_best[8];
	int i,elem,ib,jb,i_best,n_best,m_i_best[8];
	PointXYZ Pmin,Pmax,Pcur,Hex[8],p,loc;
	if(!kpar){return -1;}
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
	if(kpar<n_best){n_best=kpar;}

	for(elem=0;elem<kpar;elem++)
	{
		Pcur.x=xyz[nver[elem][0]][0];
		Pcur.y=xyz[nver[elem][0]][1];
		Pcur.z=xyz[nver[elem][0]][2];
		Hex[0]=Pmax=Pmin=Pcur;
		for(i=1;i<8;i++)
		{
			Pcur.x=xyz[nver[elem][i]][0];
			Pcur.y=xyz[nver[elem][i]][1];
			Pcur.z=xyz[nver[elem][i]][2];
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

	i_best=-1;
	sbest=1e+30;
	for(ib=0;ib<n_best;ib++)
	{
		elem=m_i_best[ib];
		if(elem==-1){break;}
		for(i=0;i<8;i++)
		{
			Hex[i].x=xyz[nver[elem][i]][0];
			Hex[i].y=xyz[nver[elem][i]][1];
			Hex[i].z=xyz[nver[elem][i]][2];
		}
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

int Vec_Prep_Data::LoadAndInitCurrent(int npls)
{
	int i,j,k,l,elem,nn;
	ifstream inf;
	ofstream outf;
	double tmpd;
	const double eps_loc_crd=1e-2;
	const double bnd_loc_min=0.0-eps_loc_crd;
	const double bnd_loc_max=1.0+eps_loc_crd;

	currentfunction.resize(ntime);

	inf.open("currentfunction");
	if(!inf)
	{
		logfile<<"Error in open file "<<"currentfunction"<<endl;
		cout<<"Error in open file "<<"sours"<<endl;
		return 1;
	}
	for(i=0;i<ntime;i++)
	{
		inf>>tmpd>>currentfunction[i];
	}
	inf.close();
	inf.clear();

	for(i=1;i<ntime;i++)
	{
		if(currentfunction[i])
		{
			break;
		}
	}

	if(i==ntime)
	{
		return 0;
	}

	nSrsRazb=500;
	inf.open("nSrsRazb");
	if(inf)
	{
		inf>>nSrsRazb;
		inf.close();
	}
	inf.clear();

	srs.resize(npls);

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
	}
	inf.close();
	inf.clear();

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

	for(elem=0;elem<kpar;elem++)
	{
		PointXYZ Pmin,Pmax,Pcur,Hex[8];	// Pmin,Pmax -  

		Pcur.x=xyz[nver[elem][0]][0];
		Pcur.y=xyz[nver[elem][0]][1];
		Pcur.z=xyz[nver[elem][0]][2];
		Hex[0]=Pmax=Pmin=Pcur;
		for(i=1;i<8;i++)
		{
			Pcur.x=xyz[nver[elem][i]][0];
			Pcur.y=xyz[nver[elem][i]][1];
			Pcur.z=xyz[nver[elem][i]][2];
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
				elem=GetNearestElement(ls);
				if(elem!=-1)
				{
					PointXYZ p,dir,pp,pm,Hex[8];
					const double delta=0.1;
					for(l=0;l<8;l++)
					{
						Hex[l].x=xyz[nver[elem][l]][0];
						Hex[l].y=xyz[nver[elem][l]][1];
						Hex[l].z=xyz[nver[elem][l]][2];
					}
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
	}

	outf.open("sours_dip_calc");
	for(i=0;i<npls;i++)
	{
		k=(int)srs[i].size();
		for(j=0;j<k;j++)
		{
			loc_source &ls=srs[i][j];
			outf<<ls.x<<' '<<ls.y<<' '<<ls.z<<' ';
			outf<<ls.ksi<<' '<<ls.eta<<' '<<ls.dzeta<<' ';
			outf<<ls.elem+1<<' '<<nvkat[ls.elem]+1<<' '<<'\n';
		}
	}
	outf.close();
	outf.clear();

	return 0;
}
