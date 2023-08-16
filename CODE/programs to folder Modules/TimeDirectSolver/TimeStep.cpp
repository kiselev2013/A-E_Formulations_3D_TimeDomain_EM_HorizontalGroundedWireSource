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
 *  This file contains the code for converting the time grid from a variable step to a grid with piecewise constant steps
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin                         
 *  Novosibirsk State Technical University,                 
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia       
 *  vdv_wk@mail.ru                                          
 *  Version 2.0 January 16, 2023                            
*/

#include "stdafx.h"
#include "TimeStep.h"

extern ofstream logfile;
TimeStep::TimeStep(void)
{
}
TimeStep::~TimeStep(void)
{
}
int TimeStep::Replace_infite0()
{
	

	return 0;
}
int TimeStep::Read_infite0(char *fname)
{
	FILE *fp=NULL;
	int i;
	char buffer[20];
	int ntime;


	if((fp=fopen(fname, "r"))==0)
		return 1;

	while(!feof(fp))
	{
		fscanf(fp, "%s", buffer);

		if(strcmp(buffer,"ktime=")==0)
		{
			fscanf(fp, "%ld", &ntime);
			nonRegTimes.resize(ntime);
			continue;
		}

		if(strcmp(buffer,"T")==0)
		{
			for(i=0; i<4; i++) //   "I M E :"
				fscanf(fp, "%s", buffer);

			for (i=0; i<ntime; i++)
			{
				fscanf(fp, "%lf", &nonRegTimes[i]); //    
				fscanf(fp, "%s", buffer); //  ";"
			}
			break;
		}
	}

	fclose(fp);
	return 0;
}
void TimeStep::GenRegMesh(double coefTime)
{
	int i, j, k;
	double stepOld;
	double step;




	stepOld = nonRegTimes[1] - nonRegTimes[0];
	intBounds.push_back(nonRegTimes[0]);
	hInt.push_back(stepOld);

	k = 2;
	while (k < nonRegTimes.size())
	{
		step = nonRegTimes[k] - nonRegTimes[k-1];

		if (step > coefTime*stepOld  || k+1 == nonRegTimes.size())
		{
			nSubInt.push_back((nonRegTimes[k] - intBounds[intBounds.size()-1])/stepOld + 1);

			intBounds.push_back(intBounds[intBounds.size()-1] + nSubInt[nSubInt.size()-1]*stepOld);

			stepOld *= coefTime;
			hInt.push_back(stepOld);
		}
		k++;
	}
	
	ig_times.push_back(0);
	gg_times.push_back(nonRegTimes[0]);

	for (i=0; i<intBounds.size()-1; i++)
	{
		for (j=0; j<nSubInt[i]; j++)
		{
			double time;
			time = intBounds[i] + (j+1)*hInt[i];
			gg_times.push_back(time);
		}
		ig_times.push_back(gg_times.size());
	}


	while(true)
	{
		if(gg_times[gg_times.size()-2] > nonRegTimes[nonRegTimes.size()-1])
		{
			if (nSubInt[nSubInt.size()-1] == 1)
			{
				nSubInt.pop_back();
				gg_times.pop_back();
				ig_times.pop_back();
			}
			else
			{
				nSubInt[nSubInt.size()-1]--;
				gg_times.pop_back();
				ig_times[ig_times.size()-1]--;
			}
		}
		else
		{
			break;
		}
	}
}
int TimeStep::WriteTimes(char *fname)
{
	int i,j,nt;
	ofstream f(fname);

	nt=(int)gg_times.size();



	f<<scientific<<setprecision(7);
	f<<"    ktime=  "<<nt<<";  ntstop=  "<<nt<<";   kiter=    1;   ntime=    1;   niter=    1;"<<'\n';
    f<<"   kprogm=    1;  kpropi=    1;  kpropt=   -1;  kitrel=  250;              ;"<<'\n';
    f<<"       u0=   0.00000e+000"<<'\n';
	f<<" T I M E :"<<'\n';
	for(i=0;i<nt;i++)
	{
		for(j=0;j<5;j++)
		{
			if(!(i+j<nt)){break;}
			f<<" "<<gg_times[i+j]<<';';
		}
		f<<'\n';
		i+=4;
	}
	f.close();
	f.clear();

	return 0;
}
int TimeStep::WriteTimes2(char *fname)
{
	int i;
	ofstream f(fname);

	f << scientific << gg_times[0] << " 1" << endl;

	for (i=1; i<gg_times.size(); i++)
	{
		f << scientific << gg_times[i] << " 0" << endl;
	}

	f.close();

	return 0;
}
int TimeStep::WriteDoubleNonRegMesh(char *fname)
{
	int i;
	ofstream f(fname);

	f << nonRegTimes.size()*2-1 << endl;

	f << scientific << nonRegTimes[0] << ";" << endl;

	for (i=1; i<nonRegTimes.size(); i++)
	{
		f << scientific << 0.5*(nonRegTimes[i]+nonRegTimes[i-1]) << ";" << endl;
		f << scientific << nonRegTimes[i] << ";" << endl;
	}

	f.close();

	return 0;
}
int TimeStep::WriteDoubleRegMesh(char *fname)
{
	int i;
	ofstream f(fname);

	f << gg_times.size()*2-1 << endl;

	f << scientific << gg_times[0] << ";" << endl;

	for (i=1; i<gg_times.size(); i++)
	{
		f << scientific << 0.5*(gg_times[i]+gg_times[i-1]) << ";" << endl;
		f << scientific << gg_times[i] << ";" << endl;
	}

	f.close();

	return 0;
}
int TimeStep::WriteDoubleRegMesh2(char *fname)
{
	int i;
	ofstream f(fname);

	f << scientific << gg_times[0] << " 1" << endl;

	for (i=1; i<gg_times.size(); i++)
	{
		f << scientific << 0.5*(gg_times[i]+gg_times[i-1]) << " 0" << endl;
		f << scientific << gg_times[i] << " 0" << endl;
	}

	f.close();

	return 0;
}
int TimeStep::WriteDoubleNonRegMesh2(char *fname)
{
	int i;
	ofstream f(fname);

	f << scientific << nonRegTimes[0] << " 1" << endl;

	for (i=1; i<nonRegTimes.size(); i++)
	{
		f << scientific << 0.5*(nonRegTimes[i]+nonRegTimes[i-1]) << " 0" << endl;
		f << scientific << nonRegTimes[i] << " 0" << endl;
	}

	f.close();

	return 0;
}
int TimeStep::WriteMeshes(char *f_old, char *f_new)
{
	int i;

	ofstream fold(f_old);
	ofstream fnew(f_new);


	fold << endl << endl;
	for (i=0; i<nonRegTimes.size(); i++)
	{
		fold << nonRegTimes[i] << "\t1\n";
	}


	fnew << endl << endl;
	for (i=0; i<gg_times.size(); i++)
	{
		fnew << gg_times[i] << "\t1.02\n";
	}

	fold.close();
	fnew.close();

	return 0;
}
void TimeStep::Calc3LayersCoef(double t, double t_1, double t_2, double *gamma)
{
	double dt[3];

	dt[0] = t - t_2;
	dt[1] = t_1 - t_2;
	dt[2] = t - t_1;

	gamma[0] = (dt[0] + dt[2]) / (dt[0] * dt[2]);
	gamma[2] = -dt[2] / (dt[0]*dt[1]);
	gamma[1] =  dt[0] / (dt[1]*dt[2]);
}
double TimeStep::Parabola(int n, double x, double x0, double x1, double x2, double y0, double y1, double y2)
{
	double y=-1e+30;

	switch (n)
	{
	case 0:
		y = 2.0*(x - 0.5)*(x - 1);
		break;

	case 1:
		y = -4.0*x*(x - 1);
		break;

	case 2:
		y = 2.0*x*(x - 0.5);
		break;

	default:
		logfile << "error\n";
		cout << "error\n";
		exit(1);
	}

	return y;
}
double TimeStep::SmoothParabolic(double x, double x0, double x1, double x2, double x3, double y0, double y1, double y2, double y3)
{
	double f;
	double f1;
	double f2;
	double a1;
	double a2;

	if (x < x0 || x > x3)
	{
		logfile << "error\n";
		cout << "error\n";
		exit(1);
	}

	if (x >= x0 && x <= x1)
	{
		f1 = Parabola(0, x, x0, x1, x2, y0, y1, y2) +
			 Parabola(1, x, x0, x1, x2, y0, y1, y2) +
			 Parabola(2, x, x0, x1, x2, y0, y1, y2);

		return f1;
	}

	if (x >= x2 && x <= x3)
	{
		f2 = Parabola(0, x, x1, x2, x3, y1, y2, y3) +
			 Parabola(1, x, x1, x2, x3, y1, y2, y3) +
			 Parabola(2, x, x1, x2, x3, y1, y2, y3);

		return f2;
	}

	if (x > x1 && x < x2)
	{
		f1 = Parabola(0, x, x0, x1, x2, y0, y1, y2) +
			 Parabola(1, x, x0, x1, x2, y0, y1, y2) +
			 Parabola(2, x, x0, x1, x2, y0, y1, y2);

		f2 = Parabola(0, x, x1, x2, x3, y1, y2, y3) +
			 Parabola(1, x, x1, x2, x3, y1, y2, y3) +
			 Parabola(2, x, x1, x2, x3, y1, y2, y3);

		a1 = (x2 - x)/(x2 - x1); 
		a2 = 1.0 - a1; 

		f = a1*f1 + a2*f2;
	}

	logfile << "error\n";
	cout << "error\n";
	exit(1);

	return f;
}
double TimeStep::dA_dt_3(double t, double u_j, double u_j1, double u_j2,
									double t_j, double t_j1, double t_j2)
{
	double du_dt;

	double dt = t_j - t_j2; 
	double dt0 = t_j - t_j1;
	double dt1 = t_j1 - t_j2;



	double mt0,mt1,mt2;

	mt0 = (dt + dt0)/(dt*dt0);
	mt1 = -dt/(dt1*dt0);
	mt2 = dt0/(dt*dt1);

	du_dt=mt0*u_j+mt1*u_j1+mt2*u_j2;

	return du_dt;
}
double TimeStep::dA_dt_2(double u_j, double u_j1, double t_j, double t_j1)
{
	double du_dt;

	du_dt = (u_j - u_j1)/(t_j - t_j1);

	return du_dt;
}

