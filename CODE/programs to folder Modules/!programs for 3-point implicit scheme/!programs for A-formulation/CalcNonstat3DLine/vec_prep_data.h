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
 *  This file contains the headers for reading and storing edge mesh in 3D VFEM
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin                         
 *  Novosibirsk State Technical University,                 
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia       
 *  vdv_wk@mail.ru                                          
 *  Version 2.0 January 16, 2023                            
*/

#pragma once
#include "T_Mapping.h"
#include "PointXYZ.h"
#include "CheckInHex.h"

struct loc_source
{
	int elem;					//  
	double x,y,z;				//     
	double dx,dy,dz;			//     
	double ksi,eta,dzeta;		//     
	double d_ksi,d_eta,d_dzeta;	//    
	double len;					//   
	loc_source(){elem=-1;}
};

struct LineSource
{
	PointXYZ A,B;
};

// The class contains procedures for reading initial data on a non-stationary 3D VFEM problem
class Vec_Prep_Data
{
public:

	struct EnLine
	{
		int mtr;
		double es[3];
		double ec[3];
	};

	Vec_Prep_Data();
	~Vec_Prep_Data();

	int Read_prep_data(); //      
	int ReadPrepDataHarmLoop(char *pointres_fname); //       
	int Read_mtz_1d(); //  usin.dat, ucos.dat, alfa, nu (  )

	int Read_mesh_for_nonstat_problem(char *pointres_fname); //    
	int Read_infite0();  //  -     infite.0

	int nSrsRazb;
	vector<LineSource> lsrs;
	vector<vector<loc_source>> srs;	//    
	vector<double> currentfunction;

	int LoadAndInitCurrent(int npls);
	int GetNearestElement(loc_source &ls);

	int Read_3dmeshregular(int interval); 

	int LoadReceivers(char *fname, int& _n_pointres, double (*(&_pointres))[3]);

	int maxiter; //      ( config)
	double eps;   // epsilon   ( config)
	int n_materials; //   
	double *mu3d;     //(mu      (    mu3d))
	double *mu0;      //(mu   (    mu3d))
	int n_pointresB;  //    B
	int n_pointresE;  //    E
	double (*pointresB)[3]; //  
	double (*pointresE)[3]; //  
	double *sigma3d;       //(sigma      (   sig3d))
	double *sigma0;        //(sigma   (   sig3d))
	double *dpr3d;       //(dpr      (   dpr3d))
	double *dpr0;        //(dpr   (   dpr3d))
	int kuzlov;      //   (,     )
	int kpar;        //    
	int kt1;         //     
	int *l13d;       //     
	int (*nver)[14]; //     (13- )
	int *nvkat;       //    -
	double (*xyz)[3];  //  
	double (*xyzt)[3];  //  
	int n_layers_1d;  // -   ()   ( sreda1d.ay)
	double *layers_1d; //    ( sreda1d.ay)
	double *sigma_1d;  // sigma0   ( sreda1d.ay)

	int UnloadxyzVectorE0ForLine(const int &n, int (*ed)[25], int (*edges)[2]);
	int LoadVectorE0ForLine();
	vector< vector<EnLine> > EnForLine;

	int N_X,N_Y,N_Z;
	vector<int> regular;
	vector<double> Xcrd,Ycrd,Zcrd;

	int n_mesh_regular_x;  //    x  3dmeshregular
	int n_mesh_regular_y;  //    y  3dmeshregular
	double *mesh_regular_x; // x-  3dmeshregular
	double *mesh_regular_y; // y-  3dmeshregular

	int ntime;     //   
	double *time;   //  

 	int norvect; //    GMRES
 	double nu;    // 
 	int alfa;    //  :  x - alpha=1;  y - alpha=0; J=(alpha*Jx,(1-alpha)*Jy, 0)
 	double *usin; //    (sin-)
 	double *ucos; //    (sin-)
 	double *z_1d; //  
 	int n_1d;    //      

	T_Mapping_Vec *tmap;
	int tasktype;

	int npr, nfreq, fdirect;
};
