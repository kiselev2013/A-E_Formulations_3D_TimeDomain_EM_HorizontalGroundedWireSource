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
 *  This file contains the headers for converting the time grid from a variable step to a grid with piecewise constant steps
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin                         
 *  Novosibirsk State Technical University,                 
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia       
 *  vdv_wk@mail.ru                                          
 *  Version 2.0 January 16, 2023                            
*/

#pragma once
class TimeStep
{
public:
	vector<double> nonRegTimes;

	vector<double> intBounds;
	vector<double> hInt;
	vector<int> nSubInt;
	vector<int> ig_times;
	vector<double> gg_times;

	int Read_infite0(char *fname);
	int Replace_infite0();
	void GenRegMesh(double coefTime);
	int WriteMeshes(char *f_old, char *f_new);
	void Calc3LayersCoef(double t, double t_1, double t_2, double *gamma);
	double SmoothParabolic(double x, double x0, double x1, double x2, double x3, double y0, double y1, double y2, double y3);
	double Parabola(int n, double x, double x0, double x1, double x2, double y0, double y1, double y2);
	double dA_dt_3(double t, double u_j, double u_j1, double u_j2,
		double t_j, double t_j1, double t_j2);
	double dA_dt_2(double u_j, double u_j1, double t_j, double t_j1);

	int WriteTimes(char *fname);
	int WriteTimes2(char *fname);


	int WriteDoubleRegMesh(char *fname);
	int WriteDoubleRegMesh2(char *fname);
	int WriteDoubleNonRegMesh(char *fname);
	int WriteDoubleNonRegMesh2(char *fname);



	TimeStep(void);
	~TimeStep(void);
};
