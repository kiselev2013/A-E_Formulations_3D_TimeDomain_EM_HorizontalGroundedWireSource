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
 *  This file contains the code for starting calculation of the grounded radial source part of primary field for nonstationary two-layer scheme task
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

int CalculateAllDecades=0;

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


bool CheckStop(void)
{
	bool fstop;
	ifstream ifstop;
	fstop=false;
	
	ifstop.open("stop");
	if(ifstop)
	{
		fstop=true;
		ifstop.close();
	}
	ifstop.clear();

	if(!ifstop)
	{
		ifstop.open("..\\stop");
		if(ifstop)
		{
			fstop=true;
			ifstop.close();
		}
		ifstop.clear();
	}

	return fstop;
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

int main(int argc, char **argv)
{
	ifstream inf;
	ofstream ofp;

	logfile.open("log2dNonStatHphi");

	if(argc>1)
	{
		CalculateAllDecades=1;
	}

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

	_mkdir("Hfi");

	system("copy mtr3d2d Hfi");
	system("copy sigma Hfi");
	system("copy sigmaZ Hfi");
	system("copy mu Hfi");
	system("copy dpr Hfi");
	system("copy infite.0 Hfi");
	system("copy currentfunction Hfi");
	system("copy deltafunction Hfi");
	system("copy timeintervalforprint Hfi");

	system("copy currentval1 Hfi\\currentval");
	system("copy inf2tr.dat1 Hfi\\inf2tr.dat");
	system("copy nvtr.dat1 Hfi\\nvtr.dat");
	system("copy nvkat2d.dat1 Hfi\\nvkat2d.dat");
	system("copy l1.dat1 Hfi\\l1.dat");
	system("copy rz.dat1 Hfi\\rz.dat");
	system("copy rz.txt1 Hfi\\rz.txt");
	system("copy r.dat1 Hfi\\r.dat");
	system("copy z.dat1 Hfi\\z.dat");
	system("copy tsize.dat1 Hfi\\tsize.dat");

	ConCatPath(PathToRoot,"..\\Modules\\BuildKUForHphi.exe",PathToExe);
	run(PathToExe);

	if (CalcSP()!=0)
	{
		logfile<<"Error!!!"<<endl;
		logfile.close();
		exit(1);
	}

	return 0;
}
