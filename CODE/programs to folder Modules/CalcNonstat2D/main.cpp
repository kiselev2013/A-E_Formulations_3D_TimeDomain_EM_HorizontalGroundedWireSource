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
 *  This file contains code for calculating nonstationary 2D task
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
#include "TimeWatcher.h"

char PathToRoot[1024],PathToExe[1024];

void ConCatPath(const char *str1,const char *str2,char *res)
{
	int i,j,k,n,i1,i2,l1,l2;
	l1=strlen(str1);
	l2=strlen(str2);
	if(!l1)
	{
		strcpy(res,str2);
	}
	else if(!l2)
	{
		strcpy(res,str1);
	}
	else
	{
		for(i=l1-1;i>=0;i--)
		{
			if(str1[i]!=' ' && str1[i]!='\t')
			{
				break;
			}
		}
		n=i+1;
		if(n)
		{
			if(str1[n-1]=='\\' || str1[n-1]=='/')
			{
				n--;
			}
		}
		j=0;
		if(str2[0]=='\\' || str2[0]=='/')
		{
			j=1;
		}

		if(n+(l2-j)>=1024)
		{
			logfile << "Error: long path" << '\n';
			exit(1);
		}

		for(k=0;k<n;k++)
		{
			res[k]=str1[k];
			if(res[k]=='/')
			{
				res[k]='\\';
			}
		}
		res[k]='\\';
		k++;
		for(;j<l2;j++)
		{
			res[k]=str2[j];
			k++;
		}
		res[k]='\0';
	}
}


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

int main(int argc,char **argv)
{
	ifstream inf;
	int gall,vall;
	
	logfile.open("LogCalcNonstat2d");

	inf.open("PathToRoot");
	if(!inf)
	{
		cout<<"Can't open file PathToRoot"<<endl;
		logfile<<"Can't open file PathToRoot"<<endl;
		return 1;
	}
	inf.getline(PathToRoot,1024);
	inf.close();
	inf.clear();

	AddTimeSpot("Begin CalcNonstat2d, Start SrsMesh2D.exe");

	ConCatPath(PathToRoot,"..\\Modules\\SrsMesh2D.exe",PathToExe);
	run(PathToExe);

	AddTimeSpot("Start calculation line");

	inf.open("clcnplsa");
	if(!inf)
	{
		cout << "Error in open file " << "clcnplsh" << '\n';
		logfile << "Error in open file " << "clcnplsh" << '\n';
	}
	inf>>gall;
	inf.close();
	inf.clear();

	if(gall)
	{
		AddTimeSpot("Begin CalcNonstat2dAx");
		ConCatPath(PathToRoot,"..\\Modules\\CalcNonstat2dAx.exe 1",PathToExe);
		run(PathToExe);
		AddTimeSpot("Begin CalcNonstat2dAr");
		ConCatPath(PathToRoot,"..\\Modules\\CalcNonstat2dAr.exe 1",PathToExe);
		run(PathToExe);
	}

	inf.open("clcnplsh");
	if(!inf)
	{
		cout << "Error in open file " << "clcnplsh" << '\n';
		logfile << "Error in open file " << "clcnplsh" << '\n';
	}
	inf>>vall;
	inf.close();
	inf.clear();

	if(vall)
	{
		AddTimeSpot("Begin CalcNonstat2dHphi");
		ConCatPath(PathToRoot,"..\\Modules\\CalcNonstat2dHphi.exe 1",PathToExe);
		run(PathToExe);
	}

	AddTimeSpot("Finish CalcNonstat2d");

	close_logfile();

	return 0;
}
