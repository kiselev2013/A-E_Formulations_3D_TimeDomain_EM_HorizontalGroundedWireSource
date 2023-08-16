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
 *  Main calculation module
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin                                
 *  Novosibirsk State Technical University,                        
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia              
 *  Corresponding author: vdv_wk@mail.ru                           
 *  Version 2.0 January 16, 2023                                   
*/

#include "stdafx.h"
#include "Picket.h"
#include "PointXYZ.h"
#include "TimeWatcher.h"

ofstream logfile;

extern int FindIntervalInDoubleVec(const vector<double> &vec,const double &elem);
extern double GetValue(const int i,const double& _t,const vector<double> &t,const vector<double> &v);

extern int ConCatEdsall(char *res,char *add);

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

bool isFileExists(char *fname)
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

void CloseProgramm(int code)
{
	logfile<<"Closing programm with code "<<code<<endl;
	logfile.close();
	exit(code);
}

bool CheckStop(void)
{
	bool fstop;
	ifstream ifstop;
	fstop=false;
	ifstop.open("stop");
	if(ifstop)
	{
		logfile<<"File stop exist!"<<endl;
		fstop=true;
		ifstop.close();
	}
	ifstop.clear();
	return fstop;
}

void WriteStopFile()
{
	ofstream ofa;
	ofa.open("stop");
	ofa<<'\n';
	ofa.close();
	ofa.clear();
}

int CreateProcessForEXENoWait(char *cmdline, char *workdir, PROCESS_INFORMATION& pi)
{
	int retp;
	STARTUPINFOA si;
	ZeroMemory(&si, sizeof(si));
	si.cb = sizeof(si);
	ZeroMemory(&pi, sizeof(pi));
	cout<<"Start "<<cmdline;
	if(workdir){cout<<" in "<<workdir;}
	cout<<endl;
	retp=CreateProcessA(NULL,  // No module name (use command line). 
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

int CreateProcessForEXE(char *cmdline, char *workdir)
{
	int retp;
	char str[256];
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

int run(char *cmd)
{
	int retp;
	retp=CreateProcessForEXE(cmd,NULL);
	if(retp)
	{
		cout<<"Error: "<<cmd<<" returned "<<retp<<endl;
		logfile<<"Error: "<<cmd<<" returned "<<retp<<endl;
		return 1;
	}
	return 0;
}

int run_no_wait(char *cmd)
{
	PROCESS_INFORMATION pi;
	CreateProcessForEXENoWait(cmd,NULL,pi);
	return 0;
}

void ReadFileAll(char *fname,vector<vector<PointXYZ>> &SumField,int nRec,int ktime,char *fstr)
{
	int iRec,itime;
	ifstream ifa;
	double t;

	SumField.resize(nRec);
	for(iRec=0;iRec<nRec;iRec++)
	{
		SumField[iRec].resize(ktime);
	}

	ifa.open(fname);
	if(!ifa){
		logfile<<"File "<<fname<<" not opend"<<'\n';
		CloseProgramm(1);
	}
	for(iRec=0;iRec<nRec;iRec++){
		ifa.getline(fstr,1000,'\n');
		ifa.getline(fstr,1000,'\n');
		ifa.getline(fstr,1000,'\n');
		ifa.getline(fstr,1000,'\n');
		ifa.getline(fstr,1000,'\n');
		ifa.getline(fstr,1000,'\n');
		for(itime=0;itime<ktime;itime++){
			ifa>>t>>SumField[iRec][itime].x>>SumField[iRec][itime].y>>SumField[iRec][itime].z;
		}
		ifa.ignore(1000,'\n');
	}
	ifa.close();
	ifa.clear();
}


bool FileExist(char *fname)
{
	bool flag;
	ifstream ifa;
	flag=false;
	ifa.open(fname);
	if(ifa)
	{
		flag=true;
		ifa.close();
	}
	ifa.clear();
	return flag;
}

void AddMainFieldAndSave(char *fname,vector<vector<PointXYZ>> &SumField,int nRec,int ktime,
						 char *fstr,string &np,string &pf,int ia,int ipls)
{
	int iRec,itime;
	double t,f;
	ifstream ifa;
	ofstream ofa;

	sprintf(fstr,"task_%d\\%s_%d",ia,fname,ipls);
	ifa.open(fstr);
	sprintf(fstr,"%s%s_%s_%d",fname,np.c_str(),pf.c_str(),ia);
	ofa.open(fstr);
	for(iRec=0;iRec<nRec;iRec++)
	{
		ifa.getline(fstr,1000,'\n');
		ofa<<fstr<<'\n';
		ifa.getline(fstr,1000,'\n');
		ofa<<fstr<<'\n';
		ifa.getline(fstr,1000,'\n');
		ofa<<fstr<<'\n';
		ifa.getline(fstr,1000,'\n');
		ofa<<fstr<<'\n';
		ifa.getline(fstr,1000,'\n');
		ofa<<fstr<<'\n';
		ifa.getline(fstr,1000,'\n');
		ofa<<fstr<<'\n';
		for(itime=0;itime<ktime;itime++)
		{
			ifa>>t>>f>>f>>f;
			ofa<<t<<'\t';
			ofa<<SumField[iRec][itime].z<<'\t';
			ofa<<f<<'\t';
			ofa<<f+SumField[iRec][itime].z<<'\n';
		}
		ifa.ignore(1000,'\n');
	}
	ifa.close();
	ifa.clear();
	ofa.close();
	ofa.clear();
}

int read_time_mesh(int &ntimes,vector<double> &times)
{
	int i;
	ifstream inf;
	inf.open("infite.0");
	if(!inf){return 1;}
	inf.ignore(1000, '=');
	inf>>ntimes;
	times.resize(ntimes);
	inf.ignore(1000, '\n');
	inf.ignore(1000, '\n');
	inf.ignore(1000, '\n');
	inf.ignore(1000, '\n');
	for(i=0;i<ntimes;i++)
	{
		if(!inf.good()){return 1;}
		inf>>times[i];
		inf.ignore(1000, ';');
	}
	inf.close();
	inf.clear();
	return 0;
}

int read_time_mesh_size(int &ntimes)
{
	ifstream inf;
	inf.open("infite.0");
	if(!inf){return 1;}
	inf.ignore(1000, '=');
	inf>>ntimes;
	inf.close();
	inf.clear();
	return 0;
}

int FileSize(const char *fname)
{
	struct _stat fi;
	_stat(fname,&fi);
    return fi.st_size;
}

void CurrentToDouble(char *buf,string &np,string &pf)
{
	sprintf(buf,"eall%d_%d",np.c_str(),pf.c_str());
	CopyFileA("eall",buf,false);
	remove("eall");
	sprintf(buf,"exall%d_%d",np.c_str(),pf.c_str());
	CopyFileA("exall",buf,false);
	remove("exall");
	sprintf(buf,"eyall%d_%d",np.c_str(),pf.c_str());
	CopyFileA("eyall",buf,false);
	remove("eyall");
	sprintf(buf,"ezall%d_%d",np.c_str(),pf.c_str());
	CopyFileA("ezall",buf,false);
	remove("ezall");
}

void PlaceToDouble(char *buf,int ipls_glob,string &np,string &pf)
{
	sprintf(buf,"rename dV_2d.000.%d dV_2d.000_%s_%s",ipls_glob+1,np.c_str(),pf.c_str());
	system(buf);

	sprintf(buf,"rename ddV_2d.000.%d ddV_2d.000_%s_%s",ipls_glob+1,np.c_str(),pf.c_str());
	system(buf);

	sprintf(buf,"rename ball.%d ball%s_%s",ipls_glob+1,np.c_str(),pf.c_str());
	system(buf);

	sprintf(buf,"rename eall.%d eall%s_%s",ipls_glob+1,np.c_str(),pf.c_str());
	system(buf);

	sprintf(buf,"rename dV.000.%d dV.000_%s_%s",ipls_glob+1,np.c_str(),pf.c_str());
	system(buf);

	sprintf(buf,"rename ddV.000.%d ddV.000_%s_%s",ipls_glob+1,np.c_str(),pf.c_str());
	system(buf);

	sprintf(buf,"rename exall.%d exall%s_%s",ipls_glob+1,np.c_str(),pf.c_str());
	system(buf);

	sprintf(buf,"rename eyall.%d eyall%s_%s",ipls_glob+1,np.c_str(),pf.c_str());
	system(buf);

	sprintf(buf,"rename ezall.%d ezall%s_%s",ipls_glob+1,np.c_str(),pf.c_str());
	system(buf);

	sprintf(buf,"rename bxall.%d bxall%s_%s",ipls_glob+1,np.c_str(),pf.c_str());
	system(buf);

	sprintf(buf,"rename byall.%d byall%s_%s",ipls_glob+1,np.c_str(),pf.c_str());
	system(buf);

	sprintf(buf,"rename bzall.%d bzall%s_%s",ipls_glob+1,np.c_str(),pf.c_str());
	system(buf);

	sprintf(buf,"rename edsxall.%d edsxall%s_%s",ipls_glob+1,np.c_str(),pf.c_str());
	system(buf);

	sprintf(buf,"rename edsyall.%d edsyall%s_%s",ipls_glob+1,np.c_str(),pf.c_str());
	system(buf);

	sprintf(buf,"rename edsall.%d edsall%s_%s",ipls_glob+1,np.c_str(),pf.c_str());
	system(buf);

	sprintf(buf,"rename b3dstat.%d b3dstat_%s_%s",ipls_glob+1,np.c_str(),pf.c_str());
	system(buf);

	sprintf(buf,"rename e3dstat.%d e3dstat_%s_%s",ipls_glob+1,np.c_str(),pf.c_str());
	system(buf);
}



void ReadTimeData(int &ktime,int &ntimes,vector<double> &times,double &ftime,double &ltime)
{
	int i,retp;
	ifstream inf;

	inf.open("timeintervalforprint");
	if(!inf)
	{
		cout<<"No file "<<"timeintervalforprint"<<'\n';
		logfile<<"No file "<<"timeintervalforprint"<<'\n';
		CloseProgramm(1);
	}
	inf>>ftime>>ltime;
	inf.close();
	inf.clear();

	retp=read_time_mesh(ntimes,times);
	if(retp)
	{
		cout << "Function read_time_mesh retrurned " << retp << '\n';
		logfile << "Function read_time_mesh retrurned " << retp << '\n';
		CloseProgramm(1);
	}

	ktime=0;
	for(i=0;i<ntimes;i++)
	{
		const double& tj=times[i];
		if(tj>=ftime && tj<=ltime)
		{
			ktime++;
		}
	}
}

void ReadPlaceData(int &npls_glob,int &nprf,vector<int> &ProfNum,vector<int> &ProfBegin,vector<int> &ProfEnd,
				   vector<int> &RecvPlsIgB,vector<int> &RecvPlsIgE,vector<int> &RecvPlsIgL)
{
	int i,j,p1,p2;
	int npls,ngrp;
	ifstream inf;

	npls=0;
	inf.open("group");
	if(!inf)
	{
		logfile<<"Error in open file "<<"group"<<endl;
		cout<<"Error in open file "<<"group"<<endl;
		exit(1);
	}
	inf>>ngrp;
	ProfNum.resize(ngrp);
	ProfBegin.resize(ngrp);
	ProfEnd.resize(ngrp);
	for(i=0;i<ngrp;i++)
	{
		inf>>j>>p1>>p2;
		ProfNum[i]=j-1;
		ProfBegin[i]=p1-1;
		ProfEnd[i]=p2-1;
		npls+=p2-p1+1;
	}
	inf.close();
	inf.clear();

	npls_glob=npls;
	nprf=ngrp;

	RecvPlsIgB.resize(npls+1);
	RecvPlsIgE.resize(npls+1);
	RecvPlsIgL.resize(npls+1);

	inf.open("recvsb");
	if(!inf)
	{
		logfile<<"Error in open file "<<"recvsb"<<endl;
		cout<<"Error in open file "<<"recvsb"<<endl;
		exit(1);
	}
	RecvPlsIgB[0]=0;
	for(i=0;i<npls;i++){inf>>RecvPlsIgB[i+1];}
	inf.close();
	inf.clear();

	inf.open("recvse");
	if(!inf)
	{
		logfile<<"Error in open file "<<"recvse"<<endl;
		cout<<"Error in open file "<<"recvse"<<endl;
		exit(1);
	}
	RecvPlsIgE[0]=0;
	for(i=0;i<npls;i++){inf>>RecvPlsIgE[i+1];}
	inf.close();
	inf.clear();

	inf.open("lin");
	if(!inf)
	{
		logfile<<"Error in open file "<<"lin"<<endl;
		cout<<"Error in open file "<<"lin"<<endl;
		exit(1);
	}
	RecvPlsIgL[0]=0;
	for(i=0;i<npls;i++){inf>>RecvPlsIgL[i+1];}
	inf.close();
	inf.clear();

	for(i=0;i<npls;i++)
	{
		RecvPlsIgB[i+1]+=RecvPlsIgB[i];
		RecvPlsIgE[i+1]+=RecvPlsIgE[i];
		RecvPlsIgL[i+1]+=RecvPlsIgL[i];
	}
}

void AddBRecsToERecs()
{
	int i,j,p1,p2;
	int npls,ngrp,nRecB,nRecE;
	ifstream inf;
	ofstream ofp;
	vector<int> RecvPlsB,RecvPlsE;
	vector<PointXYZ> RecB,RecE;

	npls=0;
	inf.open("group");
	if(!inf)
	{
		logfile<<"Error in open file "<<"group"<<endl;
		cout<<"Error in open file "<<"group"<<endl;
		exit(1);
	}
	inf>>ngrp;
	for(i=0;i<ngrp;i++)
	{
		inf>>j>>p1>>p2;
		npls+=p2-p1+1;
	}
	inf.close();
	inf.clear();

	RecvPlsB.resize(npls);
	RecvPlsE.resize(npls);

	inf.open("recvsb");
	if(!inf)
	{
		logfile<<"Error in open file "<<"recvsb"<<endl;
		cout<<"Error in open file "<<"recvsb"<<endl;
		exit(1);
	}
	j=0;
	for(i=0;i<npls;i++)
	{
		inf>>RecvPlsB[i];
		j+=RecvPlsB[i];
	}
	inf.close();
	inf.clear();

	if(!j)
	{
		return;
	}


	inf.open("recvse");
	if(!inf)
	{
		logfile<<"Error in open file "<<"recvse"<<endl;
		cout<<"Error in open file "<<"recvse"<<endl;
		exit(1);
	}
	for(i=0;i<npls;i++){inf>>RecvPlsE[i];}
	inf.close();
	inf.clear();

	inf.open("xyzVectorB");
	if(!inf)
	{
		logfile<<"Error in open file "<<"xyzVectorB"<<endl;
		cout<<"Error in open file "<<"xyzVectorB"<<endl;
		exit(1);
	}
	inf>>nRecB;
	RecB.resize(nRecB);
	for(i=0;i<nRecB;i++){inf>>RecB[i].x>>RecB[i].y>>RecB[i].z;}
	inf.close();
	inf.clear();

	inf.open("xyzVectorE");
	if(!inf)
	{
		logfile<<"Error in open file "<<"xyzVectorE"<<endl;
		cout<<"Error in open file "<<"xyzVectorE"<<endl;
		exit(1);
	}
	inf>>nRecE;
	RecE.resize(nRecE);
	for(i=0;i<nRecE;i++){inf>>RecE[i].x>>RecE[i].y>>RecE[i].z;}
	inf.close();
	inf.clear();

	ofp.open("xyzVectorE");
	ofp<<nRecB+nRecE<<'\n';
	nRecB=nRecE=0;
	for(j=0;j<npls;j++)
	{
		for(i=nRecB;i<(nRecB+RecvPlsB[j]);i++){ofp<<RecB[i].x<<' '<<RecB[i].y<<' '<<RecB[i].z<<'\n';}
		nRecB+=RecvPlsB[j];
		for(i=nRecE;i<(nRecE+RecvPlsE[j]);i++){ofp<<RecE[i].x<<' '<<RecE[i].y<<' '<<RecE[i].z<<'\n';}
		nRecE+=RecvPlsE[j];
	}
	ofp.close();
	ofp.clear();

	ofp.open("recvse");
	for(i=0;i<npls;i++){ofp<<RecvPlsB[i]+RecvPlsE[i]<<'\n';}
	ofp.close();
	ofp.clear();
}

void BuildErec()
{
	system("copy recvsb recvsbe");
	system("copy xyzVectorB xyzVectorErec");
}

void DeleteBReceivers(bool &flag)
{
	int i,j,p1,p2;
	int npls,ngrp,nRecB,nRecE;
	ifstream inf;
	ofstream ofp;
	vector<int> RecvPlsB,RecvPlsE;
	vector<PointXYZ> RecB,RecE;

	npls=0;
	inf.open("group");
	if(!inf)
	{
		logfile<<"Error in open file "<<"group"<<endl;
		cout<<"Error in open file "<<"group"<<endl;
		exit(1);
	}
	inf>>ngrp;
	for(i=0;i<ngrp;i++)
	{
		inf>>j>>p1>>p2;
		npls+=p2-p1+1;
	}
	inf.close();
	inf.clear();

	RecvPlsB.resize(npls);
	RecvPlsE.resize(npls);

	inf.open("recvsb");
	if(!inf)
	{
		logfile<<"Error in open file "<<"recvsb"<<endl;
		cout<<"Error in open file "<<"recvsb"<<endl;
		exit(1);
	}
	j=0;
	for(i=0;i<npls;i++)
	{
		inf>>RecvPlsB[i];
		j+=RecvPlsB[i];
	}
	inf.close();
	inf.clear();

	flag=false;

	if(!j)
	{
		return;
	}

	flag=true;
	system("copy recvse !recvse");
	system("copy xyzVectorE !xyzVectorE");
	system("copy recvsb !recvsb");
	system("copy xyzVectorB !xyzVectorB");

	inf.open("recvse");
	if(!inf)
	{
		logfile<<"Error in open file "<<"recvse"<<endl;
		cout<<"Error in open file "<<"recvse"<<endl;
		exit(1);
	}
	for(i=0;i<npls;i++){inf>>RecvPlsE[i];}
	inf.close();
	inf.clear();

	inf.open("xyzVectorB");
	if(!inf)
	{
		logfile<<"Error in open file "<<"xyzVectorB"<<endl;
		cout<<"Error in open file "<<"xyzVectorB"<<endl;
		exit(1);
	}
	inf>>nRecB;
	RecB.resize(nRecB);
	for(i=0;i<nRecB;i++){inf>>RecB[i].x>>RecB[i].y>>RecB[i].z;}
	inf.close();
	inf.clear();

	inf.open("xyzVectorE");
	if(!inf)
	{
		logfile<<"Error in open file "<<"xyzVectorE"<<endl;
		cout<<"Error in open file "<<"xyzVectorE"<<endl;
		exit(1);
	}
	inf>>nRecE;
	RecE.resize(nRecE);
	for(i=0;i<nRecE;i++){inf>>RecE[i].x>>RecE[i].y>>RecE[i].z;}
	inf.close();
	inf.clear();

	ofp.open("xyzVectorB");
	ofp<<0<<'\n';
	ofp.close();
	ofp.clear();

	ofp.open("recvsb");
	for(i=0;i<npls;i++){ofp<<0<<'\n';}
	ofp.close();
	ofp.clear();

	ofp.open("xyzVectorE");
	ofp<<nRecE-nRecB<<'\n';
	nRecE=0;
	for(j=0;j<npls;j++)
	{
		for(i=nRecE+RecvPlsB[j];i<(nRecE+RecvPlsE[j]);i++){ofp<<RecE[i].x<<' '<<RecE[i].y<<' '<<RecE[i].z<<'\n';}
		nRecE+=RecvPlsE[j];
	}
	ofp.close();
	ofp.clear();

	ofp.open("recvse");
	for(i=0;i<npls;i++){ofp<<(RecvPlsE[i]-RecvPlsB[i])<<'\n';}
	ofp.close();
	ofp.clear();
}

void DivideEReceivers()
{
	int i,j,k,m,p1,p2;
	int npls,ngrp,nRecEL,nRecE,nLin;
	ifstream inf;
	ofstream ofp;
	vector<int> RecvPlsE,RecvPlsEL;
	vector<PointXYZ> RecE,RecL;
	int NUMBEROFMNPOINTS=20;

	inf.open("numberofmnpoints");
	if(!inf)
	{
		logfile<<"Error in open file "<<"numberofmnpoints"<<endl;
		cout<<"Error in open file "<<"numberofmnpoints"<<endl;
		exit(1);
	}
	inf>>NUMBEROFMNPOINTS;
	inf.close();
	inf.clear();

	npls=0;
	inf.open("group");
	if(!inf)
	{
		logfile<<"Error in open file "<<"group"<<endl;
		cout<<"Error in open file "<<"group"<<endl;
		exit(1);
	}
	inf>>ngrp;
	for(i=0;i<ngrp;i++)
	{
		inf>>j>>p1>>p2;
		npls+=p2-p1+1;
	}
	inf.close();
	inf.clear();

	RecvPlsE.resize(npls);
	RecvPlsEL.resize(npls);

	inf.open("xyzmn");
	if(!inf)
	{
		logfile<<"Error in open file "<<"xyzmn"<<endl;
		cout<<"Error in open file "<<"xyzmn"<<endl;
		exit(1);
	}
	inf>>nLin;
	inf.close();
	inf.clear();

	nRecEL=0;
	inf.open("lin");
	if(!inf)
	{
		logfile<<"Error in open file "<<"lin"<<endl;
		cout<<"Error in open file "<<"lin"<<endl;
		exit(1);
	}
	for(i=0;i<npls;i++)
	{
		inf>>RecvPlsEL[i];
		RecvPlsEL[i]*=NUMBEROFMNPOINTS;
		nRecEL+=RecvPlsEL[i];
	}
	inf.close();
	inf.clear();

	inf.open("recvse");
	if(!inf)
	{
		logfile<<"Error in open file "<<"recvse"<<endl;
		cout<<"Error in open file "<<"recvse"<<endl;
		exit(1);
	}
	for(i=0;i<npls;i++){inf>>RecvPlsE[i];}
	inf.close();
	inf.clear();

	inf.open("xyzVectorE");
	if(!inf)
	{
		logfile<<"Error in open file "<<"xyzVectorE"<<endl;
		cout<<"Error in open file "<<"xyzVectorE"<<endl;
		exit(1);
	}
	inf>>nRecE;
	RecE.resize(nRecE);
	for(i=0;i<nRecE;i++){inf>>RecE[i].x>>RecE[i].y>>RecE[i].z;}
	inf.close();
	inf.clear();

	ofp.open("xyzVectorErec");
	ofp<<scientific<<setprecision(16);
	ofp<<nRecE-nRecEL<<'\n';
	k=0;
	for(j=0;j<npls;j++)
	{
		m=RecvPlsE[j]-RecvPlsEL[j];
		for(i=k;i<(k+m);i++){ofp<<RecE[i].x<<' '<<RecE[i].y<<' '<<RecE[i].z<<'\n';}
		k+=RecvPlsE[j];
	}
	ofp.close();
	ofp.clear();

	ofp.open("recvser");
	for(i=0;i<npls;i++){ofp<<(RecvPlsE[i]-RecvPlsEL[i])<<'\n';}
	ofp.close();
	ofp.clear();

	ofp.open("xyzVectorElin");
	ofp<<scientific<<setprecision(16);
	ofp<<nRecEL<<'\n';
	k=0;
	for(j=0;j<npls;j++)
	{
		m=RecvPlsE[j]-RecvPlsEL[j];
		for(i=k+m;i<(k+RecvPlsE[j]);i++){ofp<<RecE[i].x<<' '<<RecE[i].y<<' '<<RecE[i].z<<'\n';}
		k+=RecvPlsE[j];
	}
	ofp.close();
	ofp.clear();

	ofp.open("recvsel");
	for(i=0;i<npls;i++){ofp<<RecvPlsEL[i]<<'\n';}
	ofp.close();
	ofp.clear();
}

void ReturnBReceivers(bool flag)
{
	if(flag)
	{
		system("del /q recvse");
		system("del /q recvsb");
		system("del /q xyzVectorE");
		system("del /q xyzVectorB");

		system("copy !recvse recvse");
		system("copy !recvsb recvsb");
		system("copy !xyzVectorE xyzVectorE");
		system("copy !xyzVectorB xyzVectorB");

		system("del /q !recvse");
		system("del /q !recvsb");
		system("del /q !xyzVectorE");
		system("del /q !xyzVectorB");
	}
}

void AddMeshInfo()
{
	ifstream inf;
	int kpar,nodes_all,nodes_term,edges_all,edges_term;

	inf.open("inftry.dat");
	if(!inf){
		logfile<<"Error in open file inftry.dat"<<endl;
		exit(1);
	}
	inf.ignore(1000, '\n');
	inf.ignore(1000, '=');	inf>>nodes_all;
	inf.ignore(1000, '=');  inf>>kpar;
	inf.close();
	inf.clear();

	inf.open("tsize3d.dat");
	if(!inf){
		logfile<<"Error in open file tsize3d.dat"<<endl;
		exit(1);
	}
	inf>>nodes_term;
	inf.close();
	inf.clear();

	inf.open("tsize3d_.dat");
	if(!inf){
		logfile<<"Error in open file tsize3d_.dat"<<endl;
		exit(1);
	}
	inf>>edges_term;
	inf>>edges_all;
	inf.close();
	inf.clear();

	TimeWatcher::getInstance().AddMeshInfo(kpar,nodes_all,nodes_term,edges_all,edges_term);
}

wchar_t *convertCharArrayToLPCWSTR(const char* charArray)
{
	wchar_t* wString = new wchar_t[4096];
	MultiByteToWideChar(CP_ACP, 0, charArray, -1, wString, 4096);
	return wString;
}

int OpenMap(HANDLE &hMapFile,LPCTSTR &pBuf,int MemSize,char *EName)
{
	hMapFile = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, MemSize, convertCharArrayToLPCWSTR(EName));
	if (hMapFile == NULL || hMapFile == INVALID_HANDLE_VALUE)
	{
		printf("        (%d).\n",GetLastError());
		return 1;
	}
	pBuf = (LPTSTR)MapViewOfFile(hMapFile, FILE_MAP_ALL_ACCESS, 0, 0, MemSize);
	if (pBuf == NULL)
	{
		printf("     (%d).\n",GetLastError());
		return 1;
	}
	return 0;
}

//void GetPName(char *buf,char *PName)
//{
//	int i,j,k;
//	_getcwd(PName,255);
//	k=(int)strlen(PName);
//	j=0;
//	for(i=0;i<k;i++)
//	{
//		if(PName[i]=='\\')
//		{
//			j++;
//			if(j==2)
//			{
//				break;
//			}
//		}
//	}
//	if(j<2)
//	{
//		logfile<<"Error in Pname"<<endl;
//		exit(1);
//	}
//	j=0;
//	k++;
//	for(i=i+1;i<k;i++)
//	{
//		if(PName[i]=='\\' || PName[i]=='\0')
//		{
//			buf[j]='\0';
//			break;
//		}
//		else
//		{
//			buf[j]=PName[i];
//			j++;
//		}
//	}
//	strcpy(PName,buf);
//}

int GetNumberOfPlaces(int &npls)
{
	ifstream inf;
	int i,j,ngrp,p1,p2;
	npls=0;
	inf.open("group");
	if(!inf)
	{
		cout<<"Error in open file "<<"group"<<endl;
		return 1;
	}
	inf>>ngrp;
	for(i=0;i<ngrp;i++)
	{
		inf>>j>>p1>>p2;
		npls+=p2-p1+1;
	}
	inf.close();
	inf.clear();
	return 0;
}

void SumStatBE(int npls)
{
	char fname1[256],fname2[256],fres[256];
	int i,k,ipls;
	double t1,t2;
	ifstream inf,inf1,inf2;
	ofstream ofp;
	vector<int> recvsb,recvse;

	recvsb.resize(npls);
	recvse.resize(npls);

	inf.open("recvsb");
	if(!inf)
	{
		logfile<<"Error in open file "<<"recvsb"<<endl;
		exit(1);
	}
	for(ipls=0;ipls<npls;ipls++){inf>>recvsb[ipls];}
	inf.close();
	inf.clear();

	inf.open("recvse");
	if(!inf)
	{
		logfile<<"Error in open file "<<"recvse"<<endl;
		exit(1);
	}
	for(ipls=0;ipls<npls;ipls++){inf>>recvse[ipls];}
	inf.close();
	inf.clear();

	for(ipls=0;ipls<npls;ipls++)
	{
		sprintf(fname1,"b2dstat.%d",ipls+1);
		sprintf(fname2,"b3dstat_anom.%d",ipls+1);
		sprintf(fres,"b3dstat.%d",ipls+1);
		inf1.open(fname1);
		inf2.open(fname2);
		ofp.open(fres);
		if(inf1 && inf2)
		{
			for(i=0;i<recvsb[ipls];i++)
			{
				inf1>>t1;
				inf2>>t2;
				ofp<<(t1+t2)<<' ';
				inf1>>t1;
				inf2>>t2;
				ofp<<(t1+t2)<<' ';
				inf1>>t1;
				inf2>>t2;
				ofp<<(t1+t2)<<'\n';
			}
		}
		if(inf1){inf1.close();}
		if(inf2){inf2.close();}
		inf1.clear();
		inf2.clear();
		ofp.close();
		ofp.clear();
	}

	for(ipls=0;ipls<npls;ipls++)
	{
		sprintf(fname1,"e2dstat.%d",ipls+1);
		sprintf(fname2,"e3dstat_anom.%d",ipls+1);
		sprintf(fres,"e3dstat.%d",ipls+1);
		inf1.open(fname1);
		inf2.open(fname2);
		ofp.open(fres);
		if(inf1 && inf2)
		{
			for(i=0;i<recvsb[ipls];i++)
			{
				inf1>>t1;
				inf2>>t2;
				ofp<<(t1+t2)<<' ';
				inf1>>t1;
				inf2>>t2;
				ofp<<(t1+t2)<<' ';
				inf1>>t1;
				inf2>>t2;
				ofp<<(t1+t2)<<'\n';
			}
			k=recvse[ipls]-recvsb[ipls];
			for(i=0;i<k;i++)
			{
				inf1>>t1>>t1>>t1;
			}
		}
		if(inf1){inf1.close();}
		if(inf2){inf2.close();}
		inf1.clear();
		inf2.clear();
		ofp.close();
		ofp.clear();
	}
}

void WriteNullBE(int npls)
{
	char fres[256];
	int i,ipls;
	ifstream inf;
	ofstream ofp;
	vector<int> recvsb;

	recvsb.resize(npls);

	inf.open("recvsb");
	if(!inf)
	{
		logfile<<"Error in open file "<<"recvsb"<<endl;
		exit(1);
	}
	for(ipls=0;ipls<npls;ipls++){inf>>recvsb[ipls];}
	inf.close();
	inf.clear();

	for(ipls=0;ipls<npls;ipls++)
	{
		sprintf(fres,"b3dstat.%d",ipls+1);
		ofp.open(fres);
		for(i=0;i<recvsb[ipls];i++)
		{
			ofp<<0.0<<' '<<0.0<<' '<<0.0<<'\n';
		}
		ofp.close();
		ofp.clear();
	}

	for(ipls=0;ipls<npls;ipls++)
	{
		sprintf(fres,"e3dstat.%d",ipls+1);
		ofp.open(fres);
		for(i=0;i<recvsb[ipls];i++)
		{
			ofp<<0.0<<' '<<0.0<<' '<<0.0<<'\n';
		}
		ofp.close();
		ofp.clear();
	}
}

void WriteNullDV(int npls)
{
	char fres[256];
	int i,ipls;
	ifstream inf;
	ofstream ofp;
	vector<int> recvsb;

	recvsb.resize(npls);

	inf.open("lin");
	if(!inf)
	{
		logfile<<"Error in open file "<<"lin"<<endl;
		exit(1);
	}
	for(ipls=0;ipls<npls;ipls++){inf>>recvsb[ipls];}
	inf.close();
	inf.clear();

	for(ipls=0;ipls<npls;ipls++)
	{
		sprintf(fres,"dV_2d.000.%d",ipls+1);
		ofp.open(fres);
		for(i=0;i<recvsb[ipls];i++){ofp<<0.0<<'\n';}
		ofp.close();
		ofp.clear();

		sprintf(fres,"ddV_2d.000.%d",ipls+1);
		ofp.open(fres);
		for(i=0;i<recvsb[ipls];i++){ofp<<0.0<<'\n';}
		ofp.close();
		ofp.clear();

		sprintf(fres,"dV.000.%d",ipls+1);
		ofp.open(fres);
		for(i=0;i<recvsb[ipls];i++){ofp<<0.0<<'\n';}
		ofp.close();
		ofp.clear();

		sprintf(fres,"ddV.000.%d",ipls+1);
		ofp.open(fres);
		for(i=0;i<recvsb[ipls];i++){ofp<<0.0<<'\n';}
		ofp.close();
		ofp.clear();
	}
}

void WriteNullA0NodeFile(int npls)
{
	int ipls,i,n;
	ifstream inf;
	ofstream ofp;
	double tmpd;

	tmpd=0.0;

	inf.open("inftry.dat");
	if(!inf)
	{
		logfile<<"Error in open file "<<"inftry.dat"<<endl;
		exit(1);
	}
	inf.ignore(1000,'\n');
	inf.ignore(1000,'=');
	inf>>n;
	inf.close();
	inf.clear();

	ofp.open("a0.dat",ios::binary);
	for(ipls=0;ipls<npls;ipls++)
	{
		for(i=0;i<n;i++)
		{
			ofp.write((char *)&(tmpd),size_d);
			ofp.write((char *)&(tmpd),size_d);
			ofp.write((char *)&(tmpd),size_d);
		}
	}
	ofp.close();
	ofp.clear();
}

void GenerateDecInfo(int ntime,vector<double> &time)
{
	int iDec,nDec,i,j;
	double htpre,htcur;
	int max_dec_size;
	ofstream ofp;
	vector<int> DecIg;

	htpre=time[1]-time[0];
	DecIg.push_back(0);
	for(i=2;i<ntime;i++)
	{
		htcur=time[i]-time[i-1];
		if(fabs(1.0-htpre/htcur)>1e-3)
		{
			DecIg.push_back(i);
			htpre=htcur;
		}
	}
	DecIg.push_back(ntime);
	
	nDec=(int)DecIg.size()-1;

	ofp.open("nDec");
	ofp<<nDec<<'\n';
	ofp.close();
	ofp.clear();

	ofp.open("DecIg");
	ofp<<0<<'\n';
	for(iDec=0;iDec<nDec;iDec++){ofp<<DecIg[iDec+1]<<'\n';}
	ofp.close();
	ofp.clear();

	max_dec_size=0;
	for(iDec=0;iDec<nDec;iDec++)
	{
		j=DecIg[iDec+1]-DecIg[iDec];
		if(j>max_dec_size)
		{
			max_dec_size=j;
		}
	}
	ofp.open("max_dec_size");
	ofp<<max_dec_size<<'\n';
	ofp.close();
	ofp.clear();
}

void CalcMainTask(int &MemSize,int &npntE0,HANDLE *(&hMapFile),LPCTSTR *(&pBuf),int ntimes,vector<double> &times,bool fNeedBuildMesh, bool fNeedCalcStatTask,int fdirect)
{
	ifstream inf;
	ofstream ofp;
	char buf[256],PName[256];
	int it,nt,retp,ipls,npls;

	retp=GetNumberOfPlaces(npls);
	if(retp)
	{
		cout << "Function GetNumberOfPlaces retrurned " << retp << '\n';
		logfile << "Function GetNumberOfPlaces retrurned " << retp << '\n';
		exit(1);
	}

	if(fNeedCalcStatTask)
	{
		_mkdir("stationary");
		_chdir("stationary");
		system("copy ..\\PathToRoot PathToRoot");
		system("copy ..\\ImpStatCur ImpStatCur");
		system("copy ..\\geoprep.dat geoprep.dat");
		system("copy ..\\z_sig_2d z_sig_2d");
		system("copy ..\\lc.txt lc.txt");
		system("copy ..\\numberofabdipoles numberofabdipoles");
		system("copy ..\\numberofmnpoints numberofmnpoints");
		system("copy ..\\numberofsqdipoles numberofsqdipoles");
		system("copy ..\\SplineMesh SplineMesh");
		system("copy ..\\imesh imesh");
		system("copy ..\\BaseMeshLaunch BaseMeshLaunch");
		system("copy ..\\Contours Contours");
		system("copy ..\\ContoursMaterials ContoursMaterials");
		system("copy ..\\ContoursLatheralPoints ContoursLatheralPoints");
		system("copy ..\\MeshMatGeoDiscrete MeshMatGeoDiscrete");
		system("copy ..\\curvedTemplate.txt curvedTemplate.txt");
		system("copy ..\\template.txt template.txt");
		system("copy ..\\objectslevels objectslevels");
		system("copy ..\\layer_to_relief layer_to_relief");
		system("copy ..\\relief.* .");
		system("copy ..\\xyzVectorB xyzVectorB");
		system("copy ..\\xyzVectorE xyzVectorE");
		system("copy ..\\xyzVectorErec xyzVectorErec");
		system("copy ..\\xyzVectorElin xyzVectorElin");
		system("copy ..\\xyzmn xyzmn");
		system("copy ..\\sours sours");
		system("copy ..\\group group");
		system("copy ..\\recvsb recvsb");
		system("copy ..\\recvse recvse");
		system("copy ..\\recvser recvser");
		system("copy ..\\recvsel recvsel");
		system("copy ..\\lin lin");
		system("copy ..\\mobjects mobjects");
		system("copy ..\\mlayers mlayers");
		system("copy ..\\mobjectsAdditional mobjectsAdditional");
		system("copy ..\\mlayersAdditional mlayersAdditional");
		system("copy ..\\mobjects_arm mobjects_arm");
		system("copy ..\\inv_flags inv_flags");
		system("copy ..\\imp imp");
		system("copy ..\\Geoprepoptions Geoprepoptions");
		system("copy ..\\gpsettings gpsettings");	
		system("copy ..\\marine_settingsSD marine_settingsTD");
		system("copy ..\\dim_task dim_task");
		system("copy ..\\scalingZ scalingZ");
		system("copy ..\\spline_settings spline_settings");
		system("copy ..\\nthreads.txt nthreads.txt");
		_chdir("..");
	}

	if(fNeedBuildMesh)
	{
		sprintf(buf,"begin meshing");
		TimeWatcher::getInstance().AddTimeSpot(clock(),buf);

		ofp.open("TaskMainType");
		ofp<<2<<'\n';
		ofp.close();
		ofp.clear();

		ConCatPath(PathToRoot,"..\\Modules\\p.bat",PathToExe);
		system(PathToExe);
		
		if(FileExist("geoprep_orig.dat"))
		{
			system("del geoprep.dat");
			system("rename geoprep_orig.dat geoprep.dat");
		}

		sprintf(buf,"begin stationary task");
		TimeWatcher::getInstance().AddTimeSpot(clock(),buf);
	}

	if(fNeedCalcStatTask)
	{
		_chdir("stationary");

		if(fNeedBuildMesh)
		{

			ofp.open("TaskMainType");
			ofp<<2<<'\n';
			ofp.close();
			ofp.clear();

			ConCatPath(PathToRoot,"..\\Modules\\p.bat",PathToExe);
			system(PathToExe);
			
			if(FileExist("geoprep_orig.dat"))
			{
				system("del geoprep.dat");
				system("rename geoprep_orig.dat geoprep.dat");
			}

		}

		if(FileExist("AddFields"))
		{
			BuildErec();
		}

		if(!fdirect)
		{
			sprintf(buf,"Begin calc stat 2d");
			TimeWatcher::getInstance().AddTimeSpot(clock(),buf);

			ConCatPath(PathToRoot,"..\\Modules\\SrsMesh2D.exe",PathToExe);
			system(PathToExe);

			ConCatPath(PathToRoot,"..\\Modules\\CalcStat2D.exe",PathToExe);
			system(PathToExe);

			sprintf(buf,"Begin output stat 2d");
			TimeWatcher::getInstance().AddTimeSpot(clock(),buf);

			ConCatPath(PathToRoot,"..\\Modules\\OutputStat2D.exe",PathToExe);
			system(PathToExe);

			ofp.open("forA0");
			ofp<<1<<'\n';
			ofp.close();
			ofp.clear();

			sprintf(buf,"Begin BuildNormVFor3D");
			TimeWatcher::getInstance().AddTimeSpot(clock(),buf);

			ConCatPath(PathToRoot,"..\\Modules\\BuildNormVFor3D.exe",PathToExe);
			system(PathToExe);
		}

		sprintf(buf,"Begin calc stat 3D");
		TimeWatcher::getInstance().AddTimeSpot(clock(),buf);

		if(!fdirect)
		{
			ConCatPath(PathToRoot,"..\\Modules\\CalcStat3D_V.exe",PathToExe);
		}
		else
		{
			ConCatPath(PathToRoot,"..\\Modules\\CalcStat3D_3T_V.exe",PathToExe);
		}
		system(PathToExe);

		sprintf(buf,"Begin SumStat2D3D");
		TimeWatcher::getInstance().AddTimeSpot(clock(),buf);
		
		if(!fdirect)
		{
			ConCatPath(PathToRoot,"..\\Modules\\SumStat2D3D.exe",PathToExe);
			system(PathToExe);
		}
		else
		{
			for(ipls=0;ipls<npls;ipls++)
			{
				sprintf(buf,"copy dV_anom.000.%d dV.000.%d",ipls+1,ipls+1);
				system(buf);
				sprintf(buf,"copy ddV_anom.000.%d ddV.000.%d",ipls+1,ipls+1);
				system(buf);
			}
		}

		system("copy dV* ..");
		system("copy ddV* ..");

		system("copy e3dstat_anom* ..");
		system("copy b3dstat_anom* ..");

		_chdir("..");
	}

	if(!fdirect)
	{
		if(fNeedBuildMesh)
		{
			ConCatPath(PathToRoot,"..\\Modules\\SrsMesh2D.exe",PathToExe);
			system(PathToExe);
		}

		if(fNeedCalcStatTask)
		{
			sprintf(buf,"Begin calc stat 2d");
			TimeWatcher::getInstance().AddTimeSpot(clock(),buf);

			ConCatPath(PathToRoot,"..\\Modules\\CalcStat2D.exe",PathToExe);
			system(PathToExe);

			sprintf(buf,"Begin output stat 2d");
			TimeWatcher::getInstance().AddTimeSpot(clock(),buf);

			ConCatPath(PathToRoot,"..\\Modules\\OutputStat2D.exe",PathToExe);
			system(PathToExe);

			ofp.open("forA0");
			ofp<<1<<'\n';
			ofp.close();
			ofp.clear();
		
			sprintf(buf,"Begin BuildNormVFor3D");
			TimeWatcher::getInstance().AddTimeSpot(clock(),buf);

			ConCatPath(PathToRoot,"..\\Modules\\BuildNormVFor3D.exe",PathToExe);
			system(PathToExe);
		}
	}

	if(fNeedCalcStatTask)
	{
		sprintf(buf,"Begin calc stat 3D");
		TimeWatcher::getInstance().AddTimeSpot(clock(),buf);
		
		if(!fdirect)
		{
			ConCatPath(PathToRoot,"..\\Modules\\CalcStat3D_A.exe",PathToExe);
		}
		else
		{
			ConCatPath(PathToRoot,"..\\Modules\\CalcStat3D_3T_A.exe",PathToExe);
		}
		system(PathToExe);
	}
	
	if(fNeedBuildMesh)
	{
		sprintf(buf,"Begin UnloadAnomal");
		TimeWatcher::getInstance().AddTimeSpot(clock(),buf);

		if(!fdirect)
		{
			ConCatPath(PathToRoot,"..\\Modules\\UnloadAnomalNonstat.exe",PathToExe);
			system(PathToExe);
		}
		else
		{
			ofp.open("xyzVectorE0");
			ofp<<0<<'\n';
			ofp.close();
			ofp.clear();
		}
	}

	if(!fNeedCalcStatTask)
	{
		if(!FileExist("a0.edge"))
		{
			WriteNullA0NodeFile(npls);
		}
	}

	nt=ntimes;

	//GetPName(buf,PName);
	sprintf(PName,"%u",(unsigned int)time(NULL));
	logfile<<"PName= "<<PName<<endl;

	ofp.open("PName");
	ofp<<PName<<'\n';
	ofp.close();
	ofp.clear();

	npntE0=0;
	inf.open("xyzVectorE0");
	if(inf)
	{
		inf>>npntE0;
		inf.close();
	}
	inf.clear();

	MemSize = npntE0 * npls * sizeof(float);

	if(MemSize)
	{
		hMapFile = new HANDLE[nt];
		pBuf = new LPCTSTR[nt];

		for(it=0;it<nt;it++)
		{
			sprintf(buf,"m%s.enor.%04d",PName,it);
			retp=OpenMap(hMapFile[it],pBuf[it],MemSize,buf);
			if(retp)
			{
				cout<<"Function OpenMap returned "<<retp<<endl;
				logfile<<"Function OpenMap returned "<<retp<<endl;
				exit(1);
			}
		}
	}

	if(!fdirect)
	{
		sprintf(buf,"Begin calc nonstat 2d");
		TimeWatcher::getInstance().AddTimeSpot(clock(),buf);

		ConCatPath(PathToRoot,"..\\Modules\\CalcNonstat2D.exe",PathToExe);
		system(PathToExe);

		sprintf(buf,"Begin output nonstat 2d");
		TimeWatcher::getInstance().AddTimeSpot(clock(),buf);

		if(MemSize)
		{
			ConCatPath(PathToRoot,"..\\Modules\\OutputNonstat2D.exe",PathToExe);
		}
		else
		{
			ConCatPath(PathToRoot,"..\\Modules\\OutputNonstat2DOnlyRec.exe",PathToExe);
		}
		system(PathToExe);
	}

	if(MemSize || fdirect)
	{
		sprintf(buf,"Begin calc nonstat 3D");
		TimeWatcher::getInstance().AddTimeSpot(clock(),buf);

		ConCatPath(PathToRoot,"..\\Modules\\CalcNonstat3DLine.exe",PathToExe);
		system(PathToExe);

		sprintf(buf,"Begin SumNonstat2D3D");
		TimeWatcher::getInstance().AddTimeSpot(clock(),buf);

		if(!fdirect)
		{
			ConCatPath(PathToRoot,"..\\Modules\\SumNonstat2D3D.exe",PathToExe);
		}
		else
		{
			ConCatPath(PathToRoot,"..\\Modules\\SumNonstat3DOnly.exe",PathToExe);
		}
		system(PathToExe);

	}
	else
	{
		sprintf(buf,"Begin SumNonstat2DOnly");
		TimeWatcher::getInstance().AddTimeSpot(clock(),buf);

		ConCatPath(PathToRoot,"..\\Modules\\SumNonstat2DOnly.exe",PathToExe);
		system(PathToExe);
	}

	if(fNeedCalcStatTask)
	{
		if(!fdirect)
		{
			SumStatBE(npls);
		}
		else
		{
			for(ipls=0;ipls<npls;ipls++)
			{
				sprintf(buf,"copy b3dstat_anom.%d b3dstat.%d",ipls+1,ipls+1);
				system(buf);
				sprintf(buf,"copy e3dstat_anom.%d e3dstat.%d",ipls+1,ipls+1);
				system(buf);
			}
		}
	}
	else
	{
		WriteNullBE(npls);
		WriteNullDV(npls);
	}
}

void WaitTasksFinish(int nTask,vector<PROCESS_INFORMATION> &vpi)
{
	int jp,retp;
	do{
		Sleep(1000);
		for(jp=0;jp<nTask;jp++)
		{
			GetExitCodeProcess(vpi[jp].hProcess,(LPDWORD)&retp);
			if(retp==STILL_ACTIVE){break;}
		}
	}while(jp<nTask);
}

void GetLastTime(double &LastTime)
{
	int i,retp;
	ifstream inf;
	double FirstTime,tmpd;

	LastTime=1e+30;

	inf.open("t");
	if(inf)
	{
		inf>>retp;
		for(i=0;i<retp;i++)
		{
			inf>>LastTime;
		}
		inf.close();
	}
	inf.clear();

	inf.open("TimeIntervalForPrint");
	if(inf)
	{
		inf>>FirstTime>>tmpd;
		inf.close();
		if(tmpd<LastTime)
		{
			LastTime=tmpd;
		}
	}
	inf.clear();

	inf.open("imp");
	if(inf)
	{
		inf>>tmpd;
		inf>>tmpd;
		inf>>i;
		inf.close();
		if(i && tmpd<LastTime)
		{
			LastTime=tmpd;
		}
	}
	inf.clear();

	LastTime*=1.01;
}

void BuildZSig2DByRelief(char *buf,int i,vector<int> &rid,double &zave,int ia)
{
	ifstream inf;
	ofstream ofp;
	int j,npnt;
	PointXYZ TmpPnt;
	vector<PointXYZ> vp;

	npnt=0;
	if(ia==-1)
		sprintf(buf,"relief.%d",rid[i]+1);
	else
		sprintf(buf,"relief.%d.%d",rid[i]+1,ia+1);
	inf.open(buf);
	if(!inf)
	{
		logfile<<"Error: open file "<<buf<<endl;
		cout<<"Error: open file "<<buf<<endl;
		exit(1);
	}
	vp.clear();
	while(!inf.eof())
	{
		inf>>TmpPnt.x;
		if(!inf.good() || inf.eof()){break;}
		inf>>TmpPnt.y>>TmpPnt.z;
		zave+=TmpPnt.z;
		npnt++;
		vp.push_back(TmpPnt);
	}
	inf.close();
	inf.clear();

	if(npnt)
	{
		zave=0.0;
		for(j=0;j<npnt;j++)
		{
			zave+=vp[j].z;
		}
		zave/=npnt;
	}
}

int RepairGeoprepDat(char *buf,int nlay,vector<double> &Zgs,vector<double> &ZgsCpy)
{
	int i,j,k,l;
	double tmpd;
	ifstream inf;
	ofstream ofp;

	inf.open("geoprep.dat");
	if(!inf)
	{
		logfile<<"Error: open file "<<"geoprep.dat"<<endl;
		cout<<"Error: open file "<<"geoprep.dat"<<endl;
		return 1;
	}

	ofp.open("geoprep.dat.ttt");

	inf.getline(buf,1023,'\n');		ofp<<buf<<'\n';
	inf.getline(buf,1023,'\n');		ofp<<buf<<'\n';
	inf.getline(buf,1023,'\n');		ofp<<buf<<'\n';
	inf.getline(buf,1023,'\n');		ofp<<buf<<'\n';
	inf.getline(buf,1023,'\n');		ofp<<buf<<'\n';
	inf.getline(buf,1023,'\n');		ofp<<buf<<'\n';
	inf.getline(buf,1023,'\n');		ofp<<buf<<'\n';
	inf.getline(buf,1023,'\n');		ofp<<buf<<'\n';
	inf.getline(buf,1023,'\n');		ofp<<buf<<'\n';
	inf.getline(buf,1023,'\n');		ofp<<buf<<'\n';
	inf.getline(buf,1023,'\n');		ofp<<buf<<'\n';
	inf.getline(buf,1023,'\n');		ofp<<buf<<'\n';
	inf.getline(buf,1023,'\n');		ofp<<buf<<'\n';	//  

	inf>>k;		ofp<<k;	// -
	inf.getline(buf,1023,'\n');		ofp<<buf<<'\n';
	for(i=0;i<k;i++)
	{
		inf.getline(buf,1023,'\n');		ofp<<buf<<'\n';
	}
	inf.getline(buf,1023,'\n');		ofp<<buf<<'\n';	// 
	inf>>k;		ofp<<k;	// -
	inf.getline(buf,1023,'\n');		ofp<<buf<<'\n';


	vector<string> cc;
	vector<double> aa,bb;
	aa.resize(k);
	bb.resize(k);
	cc.clear();

	for(i=0;i<k;i++)
	{
		bb[i]=1e+30;
		inf>>tmpd;
		aa[i]=tmpd;
		inf.getline(buf,1023,'\n');
		cc.push_back(buf);

		for(j=0;j<nlay;j++)
		{
			if(tmpd<-1e-3)
			{
				if(fabs(1.0-ZgsCpy[j]/tmpd)<1e-3)
				{
					break;
				}
			}
			else
			{
				if(fabs(tmpd-ZgsCpy[j])<1e-3)
				{
					break;
				}
			}
		}
		if(j<nlay)
		{
			bb[i]=Zgs[j];
		}	
	}

	for(i=0;i<k;i++)
	{
		if(bb[i]==1e+30)
		{
			for(j=i-1;j>=0;j--)
			{
				if(bb[j]!=1e+30)
				{
					break;
				}
			}
			for(l=i+1;l<k;l++)
			{
				if(bb[l]!=1e+30)
				{
					break;
				}
			}
			if(j<0 || l==k)
			{
				cout<<"Error: object not in layer"<<endl;
				logfile<<"Error: object not in layer"<<endl;
				exit(1);
			}
			bb[i]=(bb[j]*(aa[l]-aa[i])+bb[l]*(aa[i]-aa[j]))/(aa[l]-aa[j]);
		}
	}

	for(i=0;i<k;i++)
	{
		ofp<<bb[i];
		ofp<<cc[i]<<'\n';
	}

	while(!inf.eof())
	{
		inf.getline(buf,1023,'\n');
		if(!inf.good() || inf.eof()){break;}
		ofp<<buf<<'\n';
	}
	inf.close();
	inf.clear();
	ofp.close();
	ofp.clear();

	sprintf(buf,"rename geoprep.dat !geoprep.dat");
	system(buf);
	sprintf(buf,"rename geoprep.dat.ttt geoprep.dat");
	system(buf);

	return 0;
}

void build_mesh_add_to_up_without_last(vector<double> &TmpVec,double hf,double kdz,double z1,double z2,bool wf)
{
	double vz,hdz;
	hdz=hf;
	vz=z1+hdz/(2.0-wf);
	while(vz+0.5*hdz<z2)
	{
		TmpVec.push_back(vz);
		hdz*=kdz;
		vz+=hdz;
	}
}

void build_mesh_add_to_up(vector<double> &TmpVec,double hf,double kdz,double z1,double z2,bool wf)
{
	double vz,hdz;
	hdz=hf;
	vz=z1+hdz/(2.0-wf);
	while(vz+0.5*hdz<z2)
	{
		TmpVec.push_back(vz);
		hdz*=kdz;
		vz+=hdz;
	}
	if(vz<z2 || (fabs(z2)<1e-6 && fabs(vz-z2)<1e-9) || fabs(1.0-vz/z2)<1e-3)
	{
		TmpVec.push_back(vz);
	}
}

void UnloadInifte0(int ktime,vector<double> &time)
{
	ofstream ofp;
	int ntstop=ktime, kiter=1, ntime=1, niter=1, kprogm=1, kpropi=1, kpropt=-1,	kitrel=250;
	float u0=0.0;

	ofp.open("infite.0");
	ofp<<"    ktime="<<setw(5)<<ktime<<";  ntstop="<<setw(5)<<ntstop<<";   kiter="<<setw(5)<<kiter<<
		";   ntime="<<setw(5)<<ntime<<";   niter="<<setw(5)<<niter<<";"<<'\n';
	ofp<<"   kprogm="<<setw(5)<<kprogm<<";  kpropi="<<setw(5)<<kpropi<<";  kpropt="<<setw(5)<<kpropt<<
		";  kitrel="<<setw(5)<<kitrel<<";              ;"<<'\n';
	ofp<<"       u0="<<setiosflags(ios_base::scientific)<<setprecision(5)<<setw(15)<<u0<<'\n';
	ofp<<" T I M E :"<<'\n';
	ofp<<setiosflags(ios_base::scientific)<<setprecision(7)<<setiosflags(ios_base::scientific);
	for(int itime=0;itime<ktime;itime++)
	{
		ofp<<setw(15)<<(float)time[itime]<<";";
		if((itime+1)%5==0){ofp<<'\n';}
	}
	if(ktime%5!=0){ofp<<'\n';}
	ofp.close();
	ofp.clear();
}

int FindTimeForStep(int ntime,vector<double> &time,double hnxt)
{
	int itime;
	double hcrr;
	for(itime=1;itime<ntime;itime++)
	{
		hcrr=time[itime]-time[itime-1];
		if(hcrr>hnxt)
		{
			return itime;
		}
	}
	return ntime-1;
}

double CheckTimeMeshDirect(double tcffTM,int ntime,vector<double> &time,double tbeg,double tfst,double tend)
{
	int itime;
	double tcur,hcur,hnxt;
	tcur=tbeg+tfst;
	hcur=tfst;
	hnxt=tcffTM*hcur;
	itime=0;
	while(itime<ntime)
	{
		itime=FindTimeForStep(ntime,time,hnxt);
		while(tcur<time[itime])
		{
			tcur+=hcur;
		}
		if(itime<ntime-1)
		{
			hcur=hnxt;
			hnxt=tcffTM*hcur;
		}
		else
		{
			break;
		}
	}
	return tend/tcur;
}

void Calc1Loop(bool fNeedBuildMesh,bool fNeedCalcStatTask,int fdirect)
{
	try
	{
		char buf[1024],fstr[512],fname[256];
		ifstream inf;
		ofstream ofp;
		PointXYZ pE,O;
		string pf,np,pft,npt;
		int ktime;
		ifstream ifa;
		ofstream ofa,off;
		bool fstop,WithDDV,WithP1,front;
		int i,it,retp;
		int mn_sea,mn_sea_type;

		int ipls_glob,npls_glob,ipls,iprf,nprf;
		vector<int> ProfBegin,ProfEnd,ProfNum;
		vector<int> RecvPlsIgB,RecvPlsIgE,RecvPlsIgL;

		int ntimes;
		vector<double> times;
		double ftime,ltime;

		int iLin,nLin,nRec;
		vector<PointXYZ> A,B,R;

		int igr;

		int MemSize,npntE0;
		HANDLE *hMapFile;
		LPCTSTR *pBuf;

		logfile.open("logfile");	

		inf.open("igroup");
		if(!inf)
		{
			cout<<"File "<<"igroup"<<" not opend"<<'\n';
			logfile<<"File "<<"igroup"<<" not opend"<<'\n';
			CloseProgramm(1);
		}
		inf>>igr;
		inf.close();
		inf.clear();

		fstop=false;
		WithDDV=WithP1=false;
		
		inf.open("inv_flags");
		if(!inf)
		{
			logfile<<"File "<<"inv_flags"<<" not opend"<<'\n';
			exit(1);
		}
		inf>>mn_sea;
		inf>>WithDDV;
		inf>>mn_sea_type;
		inf.close();
		inf.clear();

		if(CheckStop()){CloseProgramm(1);}

		DivideEReceivers();

		ConCatPath(PathToRoot,"..\\Modules\\TimeDirectSolver.exe",PathToExe);
		retp=CreateProcessForEXE(PathToExe,NULL);
		if(retp)
		{
			logfile<<"Error: "<<"TimeDirectSolver.exe "<<"returned "<<retp<<endl;
			logfile.close();
			exit(retp);
		}

		ConCatPath(PathToRoot,"..\\Modules\\CorrectTimeMesh.exe",PathToExe);
		retp=CreateProcessForEXE(PathToExe,NULL);
		if(retp)
		{
			logfile<<"Error: "<<"CorrectTimeMesh.exe "<<"returned "<<retp<<endl;
			logfile.close();
			exit(retp);
		}

		ReadTimeData(ktime,ntimes,times,ftime,ltime);

		ofp.open("ntout");
		ofp<<ktime<<'\n';
		ofp.close();
		ofp.clear();

		CalcMainTask(MemSize,npntE0,hMapFile,pBuf,ntimes,times,fNeedBuildMesh,fNeedCalcStatTask,fdirect);

		sprintf(buf,"Finish main task");
		TimeWatcher::getInstance().AddTimeSpot(clock(),buf);

		ReadPlaceData(npls_glob,nprf,ProfNum,ProfBegin,ProfEnd,RecvPlsIgB,RecvPlsIgE,RecvPlsIgL);

		nLin=RecvPlsIgL[npls_glob];
		nRec=RecvPlsIgB[npls_glob];
		A.resize(nLin);
		B.resize(nLin);
		R.resize(nRec);

		inf.open("xyzmn");
		if(!inf)
		{
			logfile<<"File "<<"xyzmn"<<" not opend"<<'\n';
			exit(1);
		}
		inf>>i;
		if(i!=nLin)
		{
			logfile<<"Number of lines in file xyzmn doesn't correspond with file lin"<<'\n';
			exit(1);
		}
		for(iLin=0;iLin<nLin;iLin++)
		{
			inf>>A[iLin].x>>A[iLin].y>>A[iLin].z>>B[iLin].x>>B[iLin].y>>B[iLin].z;
		}
		inf.close();
		inf.clear();

		inf.open("xyzVectorB");
		if(!inf)
		{
			logfile<<"File "<<"xyzVectorB"<<" not opend"<<'\n';
			exit(1);
		}
		inf>>i;
		if(i!=nRec)
		{
			logfile<<"Number of lines in file xyzmn doesn't correspond with file lin"<<'\n';
			exit(1);
		}
		for(iLin=0;iLin<nRec;iLin++)
		{
			inf>>R[iLin].x>>R[iLin].y>>R[iLin].z;
		}
		inf.close();
		inf.clear();

		if(CheckStop()){CloseProgramm(1);}

		front=false;
		inf.open("impulse");
		if(inf)
		{
			front=true;
			inf.close();
		}
		inf.clear();

		int nRecE;
		vector<PointXYZ> Ag,Bg,P;
		vector<int> RecvPlsIgEL;

		Ag.resize(npls_glob);
		Bg.resize(npls_glob);
		RecvPlsIgEL.resize(npls_glob+1);

		inf.open("recvsel");
		if(!inf)
		{
			cout<<"Error in open file "<<"recvsel"<<endl;
			logfile<<"Error in open file "<<"recvsel"<<endl;
			exit(1);
		}
		RecvPlsIgEL[0]=0;
		for(i=0;i<npls_glob;i++)
		{
			inf>>RecvPlsIgEL[i+1];
		}
		inf.close();
		inf.clear();

		for(i=1;i<=npls_glob;i++)
		{
			RecvPlsIgEL[i]+=RecvPlsIgEL[i-1];
		}

		nRecE=RecvPlsIgEL[npls_glob];
		P.resize(nRecE);

		inf.open("xyzVectorElin");
		if(!inf)
		{
			cout<<"Error in open file "<<"xyzVectorElin"<<endl;
			logfile<<"Error in open file "<<"xyzVectorElin"<<endl;
			exit(1);
		}
		inf>>i;
		if(i!=nRecE)
		{
			cout<<"Error ret!=nRecE"<<endl;
			logfile<<"Error ret!=nRecE"<<endl;
			exit(1);
		}
		for(i=0;i<nRecE;i++)
		{
			inf>>P[i].x>>P[i].y>>P[i].z;
		}
		inf.close();
		inf.clear();

		inf.open("sours");
		if(!inf)
		{
			cout<<"Error in open file "<<"sours"<<endl;
			logfile<<"Error in open file "<<"sours"<<endl;
			exit(1);
		}
		for(i=0;i<npls_glob;i++)
		{
			inf>>Ag[i].x>>Ag[i].y>>Ag[i].z;
			inf>>Bg[i].x>>Bg[i].y>>Bg[i].z;
		}
		inf.close();
		inf.clear();

		iprf=-1;
		ipls=-1;
		for(ipls_glob=0;ipls_glob<npls_glob;ipls_glob++)
		{
			if(ipls==-1)
			{
				iprf++;
				ipls=ProfBegin[iprf];
			}
			sprintf(fname,"%d",ProfNum[iprf]);
			np=fname;
			sprintf(fname,"%d",ipls);
			pf=fname;
			
			PlaceToDouble(buf,ipls_glob,np,pf);
			
			sprintf(buf,"ab%d_%d",ipls+1,ProfNum[iprf]);
			ofp.open(buf);
			ofp<<scientific<<setprecision(16);
			ofp<<Ag[ipls_glob].x<<' '<<Ag[ipls_glob].y<<' '<<Ag[ipls_glob].z<<' '
				<<Bg[ipls_glob].x<<' '<<Bg[ipls_glob].y<<' '<<Bg[ipls_glob].z<<'\n';
			ofp.close();
			ofp.clear();

			sprintf(buf,"xyzmn%d_%d",ipls+1,ProfNum[iprf]);
			ofp.open(buf);
			ofp<<scientific<<setprecision(16);
			ofp<<RecvPlsIgL[ipls_glob+1]-RecvPlsIgL[ipls_glob]<<'\n';
			for(i=RecvPlsIgL[ipls_glob];i<RecvPlsIgL[ipls_glob+1];i++)
			{
				ofp<<A[i].x<<' '<<A[i].y<<' '<<A[i].z<<' '
					<<B[i].x<<' '<<B[i].y<<' '<<B[i].z<<'\n';
			}
			ofp.close();
			ofp.clear();

			sprintf(buf,"xyzVectorElin%d_%d",ipls+1,ProfNum[iprf]);
			ofp.open(buf);
			ofp<<scientific<<setprecision(16);
			ofp<<RecvPlsIgEL[ipls_glob+1]-RecvPlsIgEL[ipls_glob]<<'\n';
			for(i=RecvPlsIgEL[ipls_glob];i<RecvPlsIgEL[ipls_glob+1];i++)
			{
				ofp<<P[i].x<<' '<<P[i].y<<' '<<P[i].z<<'\n';
			}
			ofp.close();
			ofp.clear();

			ConCatPath(PathToRoot,"..\\Modules\\SummatorEDMN.exe",PathToExe);
			sprintf(buf,"%s %d %s %s",PathToExe,WithDDV,np.c_str(),pf.c_str());
			system(buf);


			ipls++;
			if(ipls-1==ProfEnd[iprf])
			{
				ipls=-1;
			}
		}

		if(MemSize)
		{
			for(it=0;it<ntimes;it++)
			{
				UnmapViewOfFile(pBuf[it]);
				CloseHandle(hMapFile[it]);
			}
			if(hMapFile){delete [] hMapFile;hMapFile=NULL;}
			if(pBuf){delete [] pBuf;pBuf=NULL;}
		}

		iprf=-1;
		ipls=-1;
		for(ipls_glob=0;ipls_glob<npls_glob;ipls_glob++)
		{
			if(ipls==-1)
			{
				iprf++;
				ipls=ProfBegin[iprf];
			}

			sprintf(fname,"%d",ProfNum[iprf]);
			np=fname;
			sprintf(fname,"%d",ipls);
			pf=fname;

			sprintf(fname,"finish_%s_%s",np.c_str(),pf.c_str());
			off.open(fname);
			off<<"done"<<'\n';
			off.close();
			off.clear();

			fstop=CheckStop();
			if(fstop){CloseProgramm(1);}

			ipls++;
			if(ipls-1==ProfEnd[iprf])
			{
				ipls=-1;
			}
		}
	}
	catch(...)
	{
		exit(1);
	}
	logfile.close();
	logfile.clear();
}

int main(int argc, char **argv)
{
	ifstream inf;
	ofstream ofp;
	int fdirect;
	bool fNeedBuildMesh,fNeedCalcStatTask;

	setlocale(LC_ALL, "English");
	fNeedBuildMesh=true;
	fNeedCalcStatTask=true;

	fdirect=false;
	inf.open("fdirect");
	if(inf)
	{
		inf>>fdirect;
		inf.close();
	}
	inf.clear();

	_getcwd(PathToRoot,1024);

	ofp.open("PathToRoot");
	ofp<<PathToRoot<<'\n';
	ofp.close();
	ofp.clear();

	Calc1Loop(fNeedBuildMesh,fNeedCalcStatTask,fdirect);

	return 0;
}
