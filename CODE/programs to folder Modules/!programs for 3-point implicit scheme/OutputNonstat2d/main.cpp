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
 *  This file contains code for outputting nonstationary 2D task results (including primary field)
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
#include "TimeWatcher.h"
#include "MemFile.h"
#include "task2d.h"

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
	string sss;
	_getcwd(str, 1023);
	_chdir(RunPath);
	fstop = isFileExist("stop");
	if (fstop)
	{
		logfile << "File stop exist!" << endl;
	}
	sss=RunPath;
	if(!fstop && (sss.find("main")!=string::npos || sss.find("task")!=string::npos))
	{
		fstop = isFileExist("..\\stop");
		if (fstop)
		{
			logfile << "File stop exist!" << endl;
		}
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

void StopTask(PROCESS_INFORMATION &pi, int &retp)
{
	WriteStopFile();
	do{
		Sleep(1000);
		GetExitCodeProcess(pi.hProcess, (LPDWORD)&retp);
		if(retp != STILL_ACTIVE)break;
	}while(true);
}

int runnowait(char *path,PROCESS_INFORMATION &pi)
{
	if (!(CreateProcessForEXENoWait(path,".",pi)))
	{
		cout << "Can't create process" << '\n';
		logfile << "Can't create process" << '\n';
		return 1;
	}
	return 0;
}

int main(int argc,char **argv)
{
	int retp,i,j,k,l,m;
	char buff[256],PName[256];

	ifstream inf;
	ofstream ofp;

	int ipls,npls,nall,na,nh,ns;

	int it,nt;
	double *times;
	int nRec,nTim,npntE0;
	vector<int> nGels,nVels;
	vector<Picket> f_ax,f_ar,f_hp,f_sm;
	double ftime,ltime;
	vector<float> e_ax,e_ar,e_hp,e_sm;
	vector<int> vNrec;
	vector<vector<double>> Acurr,Hcurr;
	int gall,vall;
	bool fExists3D;
	vector<PointXYZ> f_ax_stat,f_ar_stat,f_hp_stat,f_sm_stat;

	logfile.open("LogOutputNonstat2d");

	_getcwd(RunPath, 1023);
	cout << "Starting programm in " << RunPath << '\n';
	logfile << "Starting programm in " << RunPath << '\n';

	AddTimeSpot("Begin Output2d");

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

	inf.open("xyzVectorE0");
	if(inf)
	{
		fExists3D=true;
		inf>>npntE0;
		inf.close();
	}
	else
	{
		fExists3D=false;
		npntE0=0;
	}
	inf.clear();

	retp=read_time_mesh(nt,times);
	if(retp)
	{
		cout << "Function read_time_mesh retrurned " << retp << '\n';
		logfile << "Function read_time_mesh retrurned " << retp << '\n';
		return 1;
	}

	gall=0;
	inf.open("clcnplsa");
	if(!inf)
	{
		cout << "Error in open file " << "clcnplsh" << '\n';
		logfile << "Error in open file " << "clcnplsh" << '\n';
	}
	inf>>gall;
	inf.close();
	inf.clear();

	vall=0;
	inf.open("clcnplsh");
	if(!inf)
	{
		cout << "Error in open file " << "clcnplsh" << '\n';
		logfile << "Error in open file " << "clcnplsh" << '\n';
	}
	inf>>vall;
	inf.close();
	inf.clear();















	AddTimeSpot("Begin outputing and summing fields");

	retp=GetNumberOfPlaces(npls);
	if(retp)
	{
		cout << "Function GetNumberOfPlaces retrurned " << retp << '\n';
		logfile << "Function GetNumberOfPlaces retrurned " << retp << '\n';
		return 1;
	}

	nGels.resize(npls);
	nVels.resize(npls);
	vNrec.resize(npls+1);
	Acurr.resize(npls);
	Hcurr.resize(npls);

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

	inf.open("srsvala");
	if(!inf)
	{
		logfile<<"Error in open file "<<"srsvala"<<endl;
		cout<<"Error in open file "<<"srsvala"<<endl;
		return 1;
	}
	for(i=0;i<npls;i++)
	{
		Acurr[i].resize(nGels[i]);
		for(j=0;j<nGels[i];j++)
		{
			inf>>Acurr[i][j];
		}
	}
	inf.close();
	inf.clear();

	inf.open("srsvalh");
	if(!inf)
	{
		logfile<<"Error in open file "<<"srsvalh"<<endl;
		cout<<"Error in open file "<<"srsvalh"<<endl;
		return 1;
	}
	for(i=0;i<npls;i++)
	{
		Hcurr[i].resize(nGels[i]+nVels[i]);
	}
	for(i=0;i<npls;i++)
	{
		for(j=0;j<nGels[i];j++)
		{
			inf>>Hcurr[i][j];
		}
	}
	for(i=0;i<npls;i++)
	{
		for(j=0;j<nVels[i];j++)
		{
			inf>>Hcurr[i][nGels[i]+j];
		}
	}
	inf.close();
	inf.clear();

	inf.open("timeintervalforprint");
	if(!inf)
	{
		logfile<<"No file "<<"timeintervalforprint"<<'\n';
		return 1;
	}
	inf>>ftime>>ltime;
	inf.close();
	inf.clear();

	nTim=0;
	for(it=0;it<nt;it++)
	{
		const double& tj=times[it];
		if(tj>=ftime && tj<=ltime)
		{
			nTim++;
		}
	}

	task2d_Ax Task_Ax;
	task2d_Ar Task_Ar;
	task2d_Hfi Task_Hfi;

	if(gall)
	{
		Task_Ax.init();
		Task_Ar.init();
	}
	if(vall)
	{
		Task_Hfi.init();
	}


	inf.open("xyzVectorE0");
	if(inf)
	{
		inf>>nRec;
		inf.close();
	}
	else
	{
		nRec=0;
	}
	inf.clear();

	if(nRec)
	{
		na=0;
		for(i=0;i<npls;i++){na+=nGels[i];}
		na*=nRec;
		e_ax.resize(na);
		e_ar.resize(na);

		nh=0;
		for(i=0;i<npls;i++){nh+=(nGels[i]+nVels[i]);}
		nh*=nRec;
		e_hp.resize(nh);

		ns=npls*nRec;
		e_sm.clear();
		e_sm.resize(ns);

		for(it=0;it<nt;it++)
		{
			if(CheckStop())
			{
				cout<<"find stop file - exit"<<endl;
				logfile<<"find stop file - exit"<<endl;
				return 1;
			}

			for(i=0;i<ns;i++){e_sm[i]=0.0;}


			if(gall)
			{





				Task_Ax.output(it,e_ax);
				Task_Ar.output(it,e_ar);
			}

			if(vall)
			{



				Task_Hfi.output(it,e_hp);
			}

			i=l=0;
			
			for(ipls=0;ipls<npls;ipls++)
			{
				for(k=0;k<nGels[ipls];k++)
				{
					for(j=0;j<nRec;j++)
					{
						e_sm[j+ipls*nRec]+=(float)Acurr[ipls][k]*e_ax[i];
						e_sm[j+ipls*nRec]+=(float)Acurr[ipls][k]*e_ar[i];
						e_sm[j+ipls*nRec]+=(float)Hcurr[ipls][k]*e_hp[l];
						i++;
						l++;
					}
				}
			}
			for(ipls=0;ipls<npls;ipls++)
			{
				for(k=0;k<nVels[ipls];k++)
				{
					for(j=0;j<nRec;j++)
					{
						e_sm[j+ipls*nRec]+=(float)Hcurr[ipls][nGels[ipls]+k]*e_hp[l];
						l++;
					}
				}
			}

			sprintf(buff,"m%s.enor.%04d",PName,it);
			WriteDataToMem(&(e_sm.front()),ns,buff);

			cout<<"Outputing time layer "<<it<<" for E0 done"<<'\n';
		}
	}
	else
	{
		for(it=0;it<nt;it++)
		{
			if(gall)
			{
				Task_Ax.output(it,e_ax);
				Task_Ar.output(it,e_ar);
			}
			if(vall)
			{
				Task_Hfi.output(it,e_hp);
			}
		}
	}



	if(gall)
	{
		Task_Ax.output_rec();
		Task_Ax.clear();

		Task_Ar.output_rec();
		Task_Ar.clear();
	}

	if(vall)
	{
		Task_Hfi.output_rec();
		Task_Hfi.clear();
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

	nall=0;
	for(i=0;i<npls;i++){nall+=nGels[i]*vNrec[i];}
	if(nall)
	{
		ReadPickets(nall,nTim,f_ax,"ball.ax",true);
		ReadPickets(nall,nTim,f_ar,"ball.ar",true);
		ReadPicketsStat(nall,f_ax_stat,"ball.ax.stat");
		ReadPicketsStat(nall,f_ar_stat,"ball.ar.stat");
	}

	nall=0;
	for(i=0;i<npls;i++){nall+=(nGels[i]+nVels[i])*vNrec[i];}
	if(nall)
	{
		ReadPickets(nall,nTim,f_hp,"ball.hphi",true);
		ReadPicketsStat(nall,f_hp_stat,"ball.hphi.stat");
	}

	m=vNrec[0];
	vNrec[0]=0;
	for(i=0;i<npls;i++)
	{
		k=vNrec[i+1];
		vNrec[i+1]=vNrec[i]+m;
		m=k;
	}

	nall=nRec;
	f_sm.clear();
	f_sm.resize(nall);
	f_sm_stat.clear();
	f_sm_stat.resize(nall);
	for(i=0;i<nall;i++)
	{
		f_sm[i].Clear();
		f_sm[i].Allocate(nTim);
		f_sm_stat[i].Clear();
		m=0;
		for(it=0;it<nt;it++)
		{
			const double& tj=times[it];
			if(tj>=ftime && tj<=ltime)
			{
				f_sm[i].t[m]=tj*1000.0;
				m++;
			}
		}
	}

	i=l=0;


	for(ipls=0;ipls<npls;ipls++)
	{
		for(k=0;k<nGels[ipls];k++)
		{
			for(j=vNrec[ipls];j<vNrec[ipls+1];j++)
			{
				f_sm[j].AddPicket(f_ax[i],Acurr[ipls][k]);
				f_sm[j].AddPicket(f_ar[i],Acurr[ipls][k]);
				f_sm[j].AddPicket(f_hp[l],Hcurr[ipls][k]);
				f_sm_stat[j].AddPoint(f_ax_stat[i],Acurr[ipls][k]);
				f_sm_stat[j].AddPoint(f_ar_stat[i],Acurr[ipls][k]);
				f_sm_stat[j].AddPoint(f_hp_stat[l],Hcurr[ipls][k]);
				i++;
				l++;
			}
		}
	}

	for(ipls=0;ipls<npls;ipls++)
	{
		for(k=0;k<nVels[ipls];k++)
		{
			for(j=vNrec[ipls];j<vNrec[ipls+1];j++)
			{
				f_sm[j].AddPicket(f_hp[l],Hcurr[ipls][nGels[ipls]+k]);
				f_sm_stat[j].AddPoint(f_hp_stat[l],Hcurr[ipls][nGels[ipls]+k]);
				l++;
			}
		}
	}

	for(ipls=0;ipls<npls;ipls++)
	{
		sprintf(buff,"ball.%d",ipls+1);
		WritePickets(vNrec[ipls],vNrec[ipls+1],nTim,f_sm,buff);
		sprintf(buff,"b2dstat.%d",ipls+1);
		WritePicketsStat(vNrec[ipls],vNrec[ipls+1],f_sm_stat,buff);
	}



	/**********************************************************/

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

	nall=0;
	for(i=0;i<npls;i++){nall+=nGels[i]*vNrec[i];}
	if(nall)
	{
		ReadPickets(nall,nTim,f_ax,"edsall.ax",true);
		ReadPickets(nall,nTim,f_ar,"edsall.ar",true);
	}

	nall=0;
	for(i=0;i<npls;i++){nall+=(nGels[i]+nVels[i])*vNrec[i];}
	if(nall)
	{
		ReadPickets(nall,nTim,f_hp,"edsall.hphi",true);
	}

	m=vNrec[0];
	vNrec[0]=0;
	for(i=0;i<npls;i++)
	{
		k=vNrec[i+1];
		vNrec[i+1]=vNrec[i]+m;
		m=k;
	}

	nall=nRec;
	f_sm.clear();
	f_sm.resize(nall);
	for(i=0;i<nall;i++)
	{
		f_sm[i].Clear();
		f_sm[i].Allocate(nTim);
		m=0;
		for(it=0;it<nt;it++)
		{
			const double& tj=times[it];
			if(tj>=ftime && tj<=ltime)
			{
				f_sm[i].t[m]=tj*1000.0;
				m++;
			}
		}
	}

	i=l=0;


	for(ipls=0;ipls<npls;ipls++)
	{
		for(k=0;k<nGels[ipls];k++)
		{
			for(j=vNrec[ipls];j<vNrec[ipls+1];j++)
			{
				f_sm[j].AddPicket(f_ax[i],Acurr[ipls][k]);
				f_sm[j].AddPicket(f_ar[i],Acurr[ipls][k]);
				f_sm[j].AddPicket(f_hp[l],Hcurr[ipls][k]);
				i++;
				l++;
			}
		}
	}

	for(ipls=0;ipls<npls;ipls++)
	{
		for(k=0;k<nVels[ipls];k++)
		{
			for(j=vNrec[ipls];j<vNrec[ipls+1];j++)
			{
				f_sm[j].AddPicket(f_hp[l],Hcurr[ipls][nGels[ipls]+k]);
				l++;
			}
		}
	}


	for(ipls=0;ipls<npls;ipls++)
	{
		sprintf(buff,"edsall_norm.%d",ipls+1);
		WritePickets(vNrec[ipls],vNrec[ipls+1],nTim,f_sm,buff);
	}


	/**********************************************************/


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

	nall=0;
	for(i=0;i<npls;i++){nall+=nGels[i]*vNrec[i];}
	if(nall)
	{
		ReadPickets(nall,nTim,f_ax,"eall.ax",true);
		ReadPickets(nall,nTim,f_ar,"eall.ar",true);
		f_ax_stat.resize(nall);
		f_ar_stat.resize(nall);
		for(i=0;i<nall;i++)
		{
			PointXYZ &pkt_ax=f_ax_stat[i];
			PointXYZ &pkt_ar=f_ar_stat[i];
			pkt_ax.x=pkt_ax.y=pkt_ax.z=0.0;
			pkt_ar.x=pkt_ar.y=pkt_ar.z=0.0;
		}
	}

	nall=0;
	for(i=0;i<npls;i++){nall+=(nGels[i]+nVels[i])*vNrec[i];}
	if(nall)
	{
		ReadPickets(nall,nTim,f_hp,"eall.hphi",true);
		ReadPicketsStat(nall,f_hp_stat,"eall.hphi.stat");
	}

	m=vNrec[0];
	vNrec[0]=0;
	for(i=0;i<npls;i++)
	{
		k=vNrec[i+1];
		vNrec[i+1]=vNrec[i]+m;
		m=k;
	}

	nall=nRec;
	f_sm.clear();
	f_sm.resize(nall);
	f_sm_stat.clear();
	f_sm_stat.resize(nall);
	for(i=0;i<nall;i++)
	{
		f_sm[i].Clear();
		f_sm[i].Allocate(nTim);
		f_sm_stat[i].Clear();
		m=0;
		for(it=0;it<nt;it++)
		{
			const double& tj=times[it];
			if(tj>=ftime && tj<=ltime)
			{
				f_sm[i].t[m]=tj*1000.0;
				m++;
			}
		}
	}

	i=l=0;



	for(ipls=0;ipls<npls;ipls++)
	{
		for(k=0;k<nGels[ipls];k++)
		{
			for(j=vNrec[ipls];j<vNrec[ipls+1];j++)
			{
				f_sm[j].AddPicket(f_ax[i],Acurr[ipls][k]);
				f_sm[j].AddPicket(f_ar[i],Acurr[ipls][k]);
				f_sm[j].AddPicket(f_hp[l],Hcurr[ipls][k]);
				f_sm_stat[j].AddPoint(f_ax_stat[i],Acurr[ipls][k]);
				f_sm_stat[j].AddPoint(f_ar_stat[i],Acurr[ipls][k]);
				f_sm_stat[j].AddPoint(f_hp_stat[l],Hcurr[ipls][k]);
				i++;
				l++;
			}
		}
	}

	for(ipls=0;ipls<npls;ipls++)
	{
		for(k=0;k<nVels[ipls];k++)
		{
			for(j=vNrec[ipls];j<vNrec[ipls+1];j++)
			{
				f_sm[j].AddPicket(f_hp[l],Hcurr[ipls][nGels[ipls]+k]);
				f_sm_stat[j].AddPoint(f_hp_stat[l],Hcurr[ipls][nGels[ipls]+k]);
				l++;
			}
		}
	}

	for(ipls=0;ipls<npls;ipls++)
	{
		sprintf(buff,"eall.%d",ipls+1);
		WritePickets(vNrec[ipls],vNrec[ipls+1],nTim,f_sm,buff);
		sprintf(buff,"e2dstat.%d",ipls+1);
		WritePicketsStat(vNrec[ipls],vNrec[ipls+1],f_sm_stat,buff);
	}


	AddTimeSpot("Finish Output2d");

	close_logfile();

	return 0;
}
