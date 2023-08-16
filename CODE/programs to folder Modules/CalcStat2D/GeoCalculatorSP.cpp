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
 *  This file contains the code for starting calculation of the primary field for stationary task
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin                        
 *  Novosibirsk State Technical University,                
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia      
 *  Corresponding author: vdv_wk@mail.ru                   
 *  Version 2.0 January 16, 2023                           
*/

#include "stdafx.h"
#include "TaskSP.h"

extern int CalcSP();
extern ofstream logfile;

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

bool isFileExists(const char *fname)
{
	ifstream inf;
	bool flag;
	flag=false;
	inf.open(fname);
	if(inf)
	{
		flag=true;
		inf.close();
	}
	inf.clear();
	return flag;
}

int CreateProcessForEXE(char *cmdline, char *workdir)
{
	char str[256];
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
		exit(retp);
	}
}

void Calc1AB()
{
	ifstream inf;

	logfile.open("LogCalcStat2D");

	inf.open("PathToRoot");
	if(!inf)
	{
		cout<<"Can't open file PathToRoot"<<endl;
		logfile<<"Can't open file PathToRoot"<<endl;
		exit(1);
	}
	inf.getline(PathToRoot,1024);
	inf.close();
	inf.clear();

	system("copy currentval4 currentval");
	system("copy inf2tr.dat4 inf2tr.dat");
	system("copy nvtr.dat4 nvtr.dat");
	system("copy nvkat2d.dat4 nvkat2d.dat");
	system("copy l1.dat4 l1.dat");
	system("copy rz.dat4 rz.dat");
	system("copy rz.txt4 rz.txt");
	system("copy r.dat4 r.dat");
	system("copy z.dat4 z.dat");
	system("copy tsize.dat4 tsize.dat");

	CalcSP();

	logfile.close();
	logfile.clear();
}

int main()
{
	setlocale( LC_ALL, "English" );
	Calc1AB();
	return 0;
}
