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
 *  This file contains code for cutting inclined source line
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

int idSth,idFrq,idGen;

double line_eps=1e-3;

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

int GetNumberOfPlaces(int &npls)
{
	int i,k,p1,p2,nprof;
	ifstream inf;

	npls=0;
	inf.open("group");
	if(!inf)
	{
		logfile<<"Error in open file "<<"group"<<endl;
		cout<<"Error in open file "<<"group"<<endl;
		return 1;
	}
	inf>>nprof;
	for(i=0;i<nprof;i++)
	{
		inf>>k>>p1>>p2;
		npls+=p2-p1+1;
	}
	inf.close();
	inf.clear();

	return 0;
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

struct LineAB
{
	double Ax,Ay,Az;
	double Bx,By,Bz;
	bool isGel(){return (fabs(Az-Bz)<=line_eps);}
	bool isVel(){return (fabs(Az-Bz)>line_eps);}
	bool isXY(){return (fabs(Ax-Bx)>line_eps && fabs(Ay-By)>line_eps);}
	bool isX(){return (fabs(Ay-By)<=line_eps);}
	bool isY(){return (fabs(Ax-Bx)<=line_eps);}
	bool isGel2()
	{
		double dx=Ax-Bx;
		double dy=Ay-By;
		double dz=fabs(Az-Bz);
		double dr=sqrt(dx*dx+dy*dy);
		if(dr>line_eps)
		{
			return (dz/dr<line_eps);
		}
		else
		{
			return false;
		}
	}
	bool isVel2()
	{
		double dx=Ax-Bx;
		double dy=Ay-By;
		double dz=fabs(Az-Bz);
		double dr=sqrt(dx*dx+dy*dy);
		if(dz>line_eps)
		{
			return (dr/dz<line_eps);
		}
		else
		{
			return false;
		}
	}
};

int nr;
double eps;
double PlCff,PvCff,LnCff;

int CutLine(LineAB &Line,vector<LineAB> &LineGel,vector<LineAB> &LineVel)
{
	ifstream inf;
	ofstream ofpg;
	ofstream ofpv;
	int i,j,k;
	double Ax,Ay,Az,Bx,By,Bz;
	double hx,hy,hz;
	int nx,ny,nz;
	int kx,ky,kz;
	vector<double> X,Y,Z;


	if(Line.isGel())
	{
		LineGel.push_back(Line);
		return 0;
	}

	if(Line.isGel2())
	{
		Line.Az=Line.Bz=0.5*(Line.Az+Line.Bz);
		LineGel.push_back(Line);
		return 0;
	}

	if(Line.isVel() && Line.isX() && Line.isY())
	{
		LineVel.push_back(Line);
		return 0;
	}

	if(Line.isVel2())
	{
		Line.Ax=Line.Bx=0.5*(Line.Ax+Line.Bx);
		Line.Ay=Line.By=0.5*(Line.Ay+Line.By);
		LineVel.push_back(Line);
		return 0;
	}


	Ax=Line.Ax;
	Ay=Line.Ay;
	Az=Line.Az;
	
	Bx=Line.Bx;
	By=Line.By;
	Bz=Line.Bz;

	if(nr)
	{
		if(fabs(Ax-Bx)>eps){hx=(Bx-Ax)/nr;kx=nr;}else{hx=(Bx-Ax);kx=0;}
		if(fabs(Ay-By)>eps){hy=(By-Ay)/nr;ky=nr;}else{hy=(By-Ay);ky=0;}
		if(fabs(Az-Bz)>eps){hz=(Bz-Az)/nr;kz=nr;}else{hz=(Bz-Az);kz=0;}

		nx=ny=nz=nr+1;

		X.resize(nx);
		Y.resize(ny);
		Z.resize(nz);

		if(kx)
			for(i=0;i<nx;i++){X[i]=Ax+i*hx;}
		else
			for(i=0;i<nx;i++){X[i]=0.5*(Ax+Bx);}

		if(ky)
			for(j=0;j<ny;j++){Y[j]=Ay+j*hy;}
		else
			for(j=0;j<ny;j++){Y[j]=0.5*(Ay+By);}

		if(kz)
			for(k=0;k<nz;k++){Z[k]=Az+k*hz;}
		else
			for(k=0;k<nz;k++){Z[k]=Bz;/*0.5*(Az+Bz);*/}

		for(i=0;i<nr;i++)
		{
			LineAB LineI;
			LineI.Ax=X[i];
			LineI.Ay=Y[i];
			LineI.Az=Z[i+1]-0.5*hz;
			LineI.Bx=X[i+1];
			LineI.By=Y[i+1];
			LineI.Bz=Z[i+1]-0.5*hz;
			LineGel.push_back(LineI);
		}

		if(kz>0)
		{
			LineAB LineI;


			
			LineI.Ax=X[0];
			LineI.Ay=Y[0];
			LineI.Az=Z[0];
			LineI.Bx=X[0];
			LineI.By=Y[0];
			LineI.Bz=Z[1]-0.5*hz;

			LineVel.push_back(LineI);
		
			for(k=1;k<kz;k++)
			{
				LineI.Ax=X[k];
				LineI.Ay=Y[k];
				LineI.Az=Z[k]-0.5*hz;
				LineI.Bx=X[k];
				LineI.By=Y[k];
				LineI.Bz=Z[k+1]-0.5*hz;

				LineVel.push_back(LineI);
			}


			LineI.Ax=X[k];
			LineI.Ay=Y[k];
			LineI.Az=Z[k]-0.5*hz;
			LineI.Bx=X[k];
			LineI.By=Y[k];
			LineI.Bz=Z[k];

			LineVel.push_back(LineI);
		}
	}
	else
	{
		LineAB LineI;

		
		hx=Ax*(1.0-PlCff)+Bx*PlCff;
		hy=Ay*(1.0-PlCff)+By*PlCff;
		hz=Az*(1.0-PvCff)+Bz*PvCff;


		LineI.Ax=Ax;
		LineI.Ay=Ay;
		LineI.Az=0.5*(Az+Bz);
		LineI.Bx=Bx;
		LineI.By=By;
		LineI.Bz=0.5*(Az+Bz);

		LineGel.push_back(LineI);
	

		LineI.Ax=hx;
		LineI.Ay=hy;
		LineI.Az=hz+0.5*(Az-Bz)*LnCff;
		LineI.Bx=hx;
		LineI.By=hy;
		LineI.Bz=hz+0.5*(Bz-Az)*LnCff;

		LineVel.push_back(LineI);
	}

	return 0;
}

int main(int argc,char **argv)
{
	int retp,i,j,k,m,sz;
	char path[256],buff[256];
	int tmpi;
	double tmpd;

	ifstream inf;
	ofstream ofp;

	int ipls,npls,nall,gall,vall;

	vector<LineAB> GenAB;
	vector<vector<LineAB>> LineGel,LineVel;

	int fnsrc;
	vector<int> nsrc;

	logfile.open("SrsMesh2d");

	_getcwd(RunPath, 1023);
	cout << "Starting programm in " << RunPath << '\n';
	logfile << "Starting programm in " << RunPath << '\n';

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

	inf.open("lc.txt");
	if(!inf)
	{
		cout<<"Can't open file lc.txt"<<endl;
		return 1;
	}
	inf>>nr;
	inf>>eps;
	inf>>line_eps;
	if(!nr)
	{
		inf>>PlCff;
		inf>>PvCff;
		inf>>LnCff;
	}
	inf.close();
	inf.clear();

	retp=GetNumberOfPlaces(npls);
	if(retp)
	{
		cout << "Function GetNumberOfPlaces retrurned " << retp << '\n';
		logfile << "Function GetNumberOfPlaces retrurned " << retp << '\n';
		return 1;
	}

	GenAB.resize(npls);
	LineGel.resize(npls);
	LineVel.resize(npls);

	for(i=0;i<npls;i++)
	{
		LineGel[i].clear();
		LineVel[i].clear();
	}

	fnsrc=0;
	inf.open("nsrc");
	if(inf)
	{
		fnsrc=1;
		nsrc.resize(npls);
		for(i=0;i<npls;i++)
		{
			inf>>nsrc[i];
		}
		inf.close();
	}
	inf.clear();

	if(!fnsrc)
	{
		inf.open("sours");
		if(!inf)
		{
			logfile<<"Error in open file "<<"sours"<<endl;
			cout<<"Error in open file "<<"sours"<<endl;
			return 1;
		}
		for(i=0;i<npls;i++){inf>>GenAB[i].Ax>>GenAB[i].Ay>>GenAB[i].Az>>GenAB[i].Bx>>GenAB[i].By>>GenAB[i].Bz;}
		inf.close();
		inf.clear();

		for(i=0;i<npls;i++)
		{
			retp=CutLine(GenAB[i],LineGel[i],LineVel[i]);
			if(retp)
			{
				cout << "Function CutLine retrurned " << retp << '\n';
				logfile << "Function CutLine retrurned " << retp << '\n';
				return 1;
			}
		}
	}
	else
	{
		for(i=0;i<npls;i++)
		{
			GenAB.resize(nsrc[i]);
			sprintf(buff,"sours.%d",i+1);
			inf.open(buff);
			for(j=0;j<nsrc[i];j++)
			{
				inf>>GenAB[j].Ax>>GenAB[j].Ay>>GenAB[j].Az>>GenAB[j].Bx>>GenAB[j].By>>GenAB[j].Bz;
				retp=CutLine(GenAB[j],LineGel[i],LineVel[i]);
				if(retp)
				{
					cout << "Function CutLine retrurned " << retp << '\n';
					logfile << "Function CutLine retrurned " << retp << '\n';
					return 1;
				}
			}
			inf.close();
			inf.clear();
		}
	}

	ofp.open("srsclcgsz");
	for(i=0;i<npls;i++)
	{
		ofp<<(int)LineGel[i].size()<<'\n';
	}
	ofp.close();
	ofp.clear();

	ofp.open("srsclcvsz");
	for(i=0;i<npls;i++)
	{
		ofp<<(int)LineVel[i].size()<<'\n';
	}
	ofp.close();
	ofp.clear();

	ofp.open("srsclcsz");
	for(i=0;i<npls;i++)
	{
		ofp<<(int)(LineGel[i].size()+LineVel[i].size())<<'\n';
	}
	ofp.close();
	ofp.clear();

	ofp.open("srsvala");
	for(i=0;i<npls;i++)
	{
		sz=(int)LineGel[i].size();
		for(j=0;j<sz;j++)
		{
			ofp<<' '<<1;
		}
		ofp<<'\n';
	}
	ofp.close();
	ofp.clear();

	ofp.open("srsvalh");
	for(i=0;i<npls;i++)
	{
		sz=(int)LineGel[i].size();
		for(j=0;j<sz;j++)
		{
			ofp<<' '<<-1;
		}
		ofp<<'\n';
	}
	for(i=0;i<npls;i++)
	{
		sz=(int)LineVel[i].size();
		for(j=0;j<sz;j++)
		{
			ofp<<' '<<((LineVel[i][j].Az>LineVel[i][j].Bz)? 1 : -1);
		}
		ofp<<'\n';
	}
	ofp.close();
	ofp.clear();

	ofp.open("srsclca");
	ofp<<scientific<<setprecision(16);
	nall=0;
	for(i=0;i<npls;i++)
	{
		sz=(int)LineGel[i].size();
		nall+=sz;
		for(j=0;j<sz;j++)
		{
			LineAB &LineI=LineGel[i][j];
			ofp<<LineI.Ax<<' '<<LineI.Ay<<' '<<LineI.Az<<' '<<LineI.Bx<<' '<<LineI.By<<' '<<LineI.Bz<<'\n';
		}
	}
	ofp.close();
	ofp.clear();

	ofp.open("clcnplsa");
	ofp<<nall<<'\n';
	ofp.close();
	ofp.clear();

	gall=nall;

	ofp.open("srsclch");
	ofp<<scientific<<setprecision(16);
	nall=0;

	for(i=0;i<npls;i++)
	{
		sz=(int)LineGel[i].size();
		nall+=sz;
		for(j=0;j<sz;j++)
		{
			LineAB &LineI=LineGel[i][j];
			ofp<<LineI.Ax<<' '<<LineI.Ay<<' '<<LineI.Az<<' '<<LineI.Bx<<' '<<LineI.By<<' '<<LineI.Bz<<'\n';
		}
	}

	for(i=0;i<npls;i++)
	{
		sz=(int)LineVel[i].size();
		nall+=sz;
		for(j=0;j<sz;j++)
		{
			LineAB &LineI=LineVel[i][j];
			ofp<<LineI.Ax<<' '<<LineI.Ay<<' '<<LineI.Az<<' '<<LineI.Bx<<' '<<LineI.By<<' '<<LineI.Bz<<'\n';
		}
	}

	ofp.close();
	ofp.clear();

	ofp.open("clcnplsh");
	ofp<<nall<<'\n';
	ofp.close();
	ofp.clear();

	vall=nall;
	
	close_logfile();

	return 0;
}
