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
 *  This file contains code for rename 3D nonstationary task results
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 January 16, 2023
*/

#include "stdafx.h"
#include "utils.h"
#include "utils_direct_solver.h"
#include "picket.h"

int idSth,idFrq,idGen;


int CreateProcessForEXE(char *cmdline, char *workdir)
{
	int retp;
	STARTUPINFOA si;
	PROCESS_INFORMATION pi;
	ZeroMemory(&si, sizeof(si));
	si.cb = sizeof(si);
	ZeroMemory(&pi, sizeof(pi));
	if (!(retp=CreateProcessA(NULL,  // No module name (use command line). 
		(LPSTR)(const char*)cmdline,// Command line. 
		NULL,				// Process handle not inheritable. 
		NULL,				// Thread handle not inheritable. 
		FALSE,				// Set handle inheritance to FALSE. 
		0/*CREATE_NO_WINDOW*/,	// No creation flags. 
		NULL,				// Use parent's environment block. 
		workdir,			// Use parent's starting directory. 
		&si,				// Pointer to STARTUPINFO structure.
		&pi)))				// Pointer to PROCESS_INFORMATION structure.
	{
		sprintf(str,"Can't create process for %s, error code %d",cmdline,GetLastError());
		cout<<str<<endl;
		logfile<<str<<endl;
		return retp;
	}
	WaitForSingleObject(pi.hProcess, INFINITE);
	GetExitCodeProcess(pi.hProcess, (LPDWORD)&retp);
	CloseHandle(pi.hProcess);
	CloseHandle(pi.hThread);
	return retp;
}

void run(char *cmd)
{
	int retp;
	retp=CreateProcessForEXE(cmd,NULL);
	if(retp)
	{
		cout<<"Error: "<<cmd<<" returned "<<retp<<endl;
		logfile<<"Error: "<<cmd<<" returned "<<retp<<endl;
		CloseProgrammWithError(retp);
	}
}

int CreateProcessForEXENoWait(char *cmdline, char *workdir, PROCESS_INFORMATION& pi)
{
	int retp;
	STARTUPINFOA si;
	ZeroMemory(&si, sizeof(si));
	si.cb = sizeof(si);
	ZeroMemory(&pi, sizeof(pi));
	cout << "Start " << cmdline;
	if (workdir){ cout << " in " << workdir; }
	cout << endl;
	retp = CreateProcessA(NULL,  // No module name (use command line). 
		(LPSTR)(const char*)cmdline,// Command line. 
		NULL,				// Process handle not inheritable. 
		NULL,				// Thread handle not inheritable. 
		FALSE,				// Set handle inheritance to FALSE. 
		CREATE_NO_WINDOW,	// No creation flags. 
		NULL,				// Use parent's environment block. 
		workdir,				// Use parent's starting directory. 
		&si,				// Pointer to STARTUPINFO structure.
		&pi);				// Pointer to PROCESS_INFORMATION structure.
	return retp;
}

bool CheckStop(void)
{
	bool fstop;
	_getcwd(str, 1023);
	_chdir(RunPath);
	fstop = isFileExist("stop");
	if (fstop)
	{
		logfile << "File stop exist!" << endl;
	}
	_chdir(str);
	return fstop;
}

void WriteStopFile()
{
	ofstream ofa;
	ofa.open("stop");
	ofa << '\n';
	ofa.close();
	ofa.clear();
}

struct LineAB
{
	double Ax,Ay,Az;
	double Bx,By,Bz;
	bool isGel(){return (fabs(Az-Bz)<=line_eps);}
	bool isVel(){return (fabs(Az-Bz)>line_eps);}
	bool isXY(){return (fabs(Ax-Bx)>line_eps && fabs(Ay-By)>line_eps);}
	bool isX(){return (fabs(Ay-By)<=line_eps);}
	bool isY(){return (fabs(Ax-Bx)<=line_eps);}
};

void _diff_t2(double *du_j,double *u_j,double *u_j1,int n,double *time,int tnum)
{
	int i;
	double invdt;
	invdt=1.0/(time[tnum]-time[tnum-1]);
	for(i=0;i<n;i++)du_j[i]=-((u_j[i]-u_j1[i])*invdt);
}

void _diff_t3(double *du_j,double *u_j,double *u_j1,double *u_j2,int n,double *time,int tnum)
{
	int i;
	double dt,dt0,dt1;
	double mt0,mt1,mt2;

	dt =  time[tnum]  - time[tnum-2];
	dt0 = time[tnum]  - time[tnum-1];
	dt1 = time[tnum-1] - time[tnum-2];

	mt0 = (dt + dt0)/(dt*dt0);
	mt1 = -dt/(dt1*dt0);
	mt2 = dt0/(dt*dt1);
	
	for(i=0;i<n;i++)du_j[i]=-(mt0*u_j[i]+mt1*u_j1[i]+mt2*u_j2[i]);
}

void _diff_t3(double *du_j,double *u_j,double *u_j1,double *u_j2,int n,double t,double t1,double t2)
{
	int i;
	double dt,dt0,dt1;
	double mt0,mt1,mt2;

	dt =  t  - t2;
	dt0 = t  - t1;
	dt1 = t1 - t2;

	mt0 = (dt + dt0)/(dt*dt0);
	mt1 = -dt/(dt1*dt0);
	mt2 = dt0/(dt*dt1);
	
	for(i=0;i<n;i++)du_j[i]=-(mt0*u_j[i]+mt1*u_j1[i]+mt2*u_j2[i]);
}

int main(int argc,char **argv)
{
	int retp,i,j,k,l,m;
	char buff[256];

	ifstream inf;
	ofstream ofp;

	int ipls,npls,nall,na,nh,ns;

	int it,nt;
	double *times;
	int nRec,nTim;
	vector<Picket> fanm,fsum;
	double ftime,ltime;
	vector<int> vNrec;
	char LegendSummB[3][128],LegendSummF[3][128],LegendSummE[3][128];
	char Bname[3][16],Fname[3][16],Ename[3][16];

	logfile.open("LogSumNonstat2D3D");

	_getcwd(RunPath, 1023);
	cout << "Starting programm in " << RunPath << '\n';
	logfile << "Starting programm in " << RunPath << '\n';

	strcpy(LegendSummB[0],"t (ms)        Bx (mTl)     Bx (mTl)     Bx (mTl)");
	strcpy(LegendSummB[1],"t (ms)        By (mTl)     By (mTl)     By (mTl)");
	strcpy(LegendSummB[2],"t (ms)        Bz (mTl)     Bz (mTl)     Bz (mTl)");

	strcpy(LegendSummF[0],"t (ms)      Emfx (mV)   Emfx (mV)   Emfx (mV)");
	strcpy(LegendSummF[1],"t (ms)      Emfy (mV)   Emfy (mV)   Emfy (mV)");
	strcpy(LegendSummF[2],"t (ms)      Emfz (mV)   Emfz (mV)   Emfz (mV)");

	strcpy(LegendSummE[0],"t (ms)        Ex (V/m)     Ex (V/m)     Ex (V/m)");
	strcpy(LegendSummE[1],"t (ms)        Ey (V/m)     Ey (V/m)     Ey (V/m)");
	strcpy(LegendSummE[2],"t (ms)        Ez (V/m)     Ez (V/m)     Ez (V/m)");

	strcpy(Bname[0],"Bxall");
	strcpy(Bname[1],"Byall");
	strcpy(Bname[2],"Bzall");

	strcpy(Fname[0],"Edsxall");
	strcpy(Fname[1],"Edsyall");
	strcpy(Fname[2],"Edsall");

	strcpy(Ename[0],"Exall");
	strcpy(Ename[1],"Eyall");
	strcpy(Ename[2],"Ezall");


	retp=GetNumberOfPlaces(npls);
	if(retp)
	{
		cout << "Function GetNumberOfPlaces retrurned " << retp << '\n';
		logfile << "Function GetNumberOfPlaces retrurned " << retp << '\n';
		return 1;
	}

	vNrec.resize(npls);

	inf.open("timeintervalforprint");
	if(!inf)
	{
		logfile<<"No file "<<"timeintervalforprint"<<'\n';
		return 1;
	}
	inf>>ftime>>ltime;
	inf.close();
	inf.clear();

	retp=read_time_mesh(nt,times);
	if(retp)
	{
		cout << "Function read_time_mesh retrurned " << retp << '\n';
		logfile << "Function read_time_mesh retrurned " << retp << '\n';
		return 1;
	}

	nTim=0;
	for(it=0;it<nt;it++)
	{
		const double& tj=times[it];
		if(tj>=ftime && tj<=ltime)
		{
			nTim++;
		}
	}

	inf.open("xyzVectorB");
	if(!inf)
	{
		logfile<<"Error in open file "<<"xyzVectorB"<<endl;
		cout<<"Error in open file "<<"xyzVectorB"<<endl;
		return 1;
	}
	inf>>nRec;
	inf.close();
	inf.clear();

	inf.open("recvsb");
	if(!inf)
	{
		logfile<<"Error in open file "<<"recvsb"<<endl;
		cout<<"Error in open file "<<"recvsb"<<endl;
		return 1;
	}
	for(i=0;i<npls;i++){inf>>vNrec[i];}
	inf.close();
	inf.clear();

	for(ipls=0;ipls<npls;ipls++)
	{
		nall=vNrec[ipls];
		sprintf(str,"ball_anom.%d",ipls+1);
		ReadPickets(nall,nTim,fanm,str,true);
		for(j=0;j<3;j++)
		{
			fsum.clear();
			fsum.resize(nall);
			for(i=0;i<nall;i++)
			{
				fsum[i].Clear();
				fsum[i].Allocate(nTim);
				m=0;
				for(it=0;it<nt;it++)
				{
					const double& tj=times[it];
					if(tj>=ftime && tj<=ltime)
					{
						fsum[i].t[m]=tj*1000.0;
						m++;
					}
				}
				fsum[i].CopyPicketInfo(fanm[i]);
				fsum[i].SetLegend(LegendSummB[j]);
				fsum[i].CopyValuesOnlyAnomal(fanm[i].v[j]);
			}
			sprintf(str,"%s.%d",Bname[j],ipls+1);
			WritePickets(0,nall,nTim,fsum,str);
		}

		sprintf(str,"edsall_anom.%d",ipls+1);
		ReadPickets(nall,nTim,fanm,str,false);
		for(j=0;j<3;j++)
		{
			for(i=0;i<nall;i++)
			{
				fsum[i].CopyPicketInfo(fanm[i]);
				fsum[i].SetLegend(LegendSummF[j]);
				fsum[i].CopyValuesOnlyAnomal(fanm[i].v[j]);
			}
			sprintf(str,"%s.%d",Fname[j],ipls+1);
			WritePickets(0,nall,nTim,fsum,str);
		}
	}

	inf.open("xyzVectorE");
	if(!inf)
	{
		logfile<<"Error in open file "<<"xyzVectorE"<<endl;
		cout<<"Error in open file "<<"xyzVectorE"<<endl;
		return 1;
	}
	inf>>nRec;
	inf.close();
	inf.clear();

	inf.open("recvse");
	if(!inf)
	{
		logfile<<"Error in open file "<<"recvse"<<endl;
		cout<<"Error in open file "<<"recvse"<<endl;
		return 1;
	}
	for(i=0;i<npls;i++){inf>>vNrec[i];}
	inf.close();
	inf.clear();

	for(ipls=0;ipls<npls;ipls++)
	{
		nall=vNrec[ipls];
		sprintf(str,"eall_anom.%d",ipls+1);
		ReadPickets(nall,nTim,fanm,str,true);
		for(j=0;j<3;j++)
		{
			fsum.clear();
			fsum.resize(nall);
			for(i=0;i<nall;i++)
			{
				fsum[i].Clear();
				fsum[i].Allocate(nTim);
				m=0;
				for(it=0;it<nt;it++)
				{
					const double& tj=times[it];
					if(tj>=ftime && tj<=ltime)
					{
						fsum[i].t[m]=tj*1000.0;
						m++;
					}
				}
				fsum[i].CopyPicketInfo(fanm[i]);
				fsum[i].SetLegend(LegendSummE[j]);
				fsum[i].CopyValuesOnlyAnomal(fanm[i].v[j]);
			}
			sprintf(str,"%s.%d",Ename[j],ipls+1);
			WritePickets(0,nall,nTim,fsum,str);
		}
	}

	close_logfile();

	return 0;
}
