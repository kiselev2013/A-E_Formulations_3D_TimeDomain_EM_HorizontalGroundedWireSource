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
 *  This file contains code for summing 2D and 3D stationary tasks results
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
	int retp,i;

	ifstream inf;
	ofstream ofp;

	int ipls,npls,nall;

	int nRec;
	double fnrm,fanm,fsum;
	vector<int> vNrec;

	ifstream infa,infn;

	logfile.open("LogSumStat2D3D");

	_getcwd(RunPath, 1023);

	cout << "Starting programm in " << RunPath << '\n';
	logfile << "Starting programm in " << RunPath << '\n';

	retp=GetNumberOfPlaces(npls);
	if(retp)
	{
		cout << "Function GetNumberOfPlaces retrurned " << retp << '\n';
		logfile << "Function GetNumberOfPlaces retrurned " << retp << '\n';
		return 1;
	}

	vNrec.resize(npls);

	inf.open("xyzmn");
	if(!inf)
	{
		logfile<<"Error in open file "<<"xyzmn"<<endl;
		cout<<"Error in open file "<<"xyzmn"<<endl;
		return 1;
	}
	inf>>nRec;
	inf.close();
	inf.clear();

	inf.open("lin");
	if(!inf)
	{
		logfile<<"Error in open file "<<"lin"<<endl;
		cout<<"Error in open file "<<"lin"<<endl;
		return 1;
	}
	for(i=0;i<npls;i++){inf>>vNrec[i];}
	inf.close();
	inf.clear();

	for(ipls=0;ipls<npls;ipls++)
	{
		nall=vNrec[ipls];

		sprintf(str,"dV_2d.000.%d",ipls+1);
		infn.open(str);
		if(!infn)
		{
			logfile<<"Error in open file "<<str<<endl;
			cout<<"Error in open file "<<str<<endl;
			return 1;
		}
		sprintf(str,"dV_anom.000.%d",ipls+1);
		infa.open(str);
		if(!infa)
		{
			logfile<<"Error in open file "<<str<<endl;
			cout<<"Error in open file "<<str<<endl;
			return 1;
		}
		sprintf(str,"dV.000.%d",ipls+1);
		ofp.open(str);
		for(i=0;i<nall;i++)
		{
			infn>>fnrm;
			infa>>fanm;
			fsum=fnrm+fanm;
			ofp<<fsum<<'\n';
		}
		ofp.close();
		ofp.clear();

		infn.close();
		infn.clear();

		infa.close();
		infa.clear();
	}

	for(ipls=0;ipls<npls;ipls++)
	{
		nall=vNrec[ipls];

		sprintf(str,"ddV_2d.000.%d",ipls+1);
		infn.open(str);
		if(!infn)
		{
			logfile<<"Error in open file "<<str<<endl;
			cout<<"Error in open file "<<str<<endl;
			return 1;
		}
		sprintf(str,"ddV_anom.000.%d",ipls+1);
		infa.open(str);
		if(!infa)
		{
			logfile<<"Error in open file "<<str<<endl;
			cout<<"Error in open file "<<str<<endl;
			return 1;
		}
		sprintf(str,"ddV.000.%d",ipls+1);
		ofp.open(str);
		for(i=0;i<nall;i++)
		{
			infn>>fnrm;
			infa>>fanm;
			fsum=fnrm+fanm;
			ofp<<fsum<<'\n';
		}
		ofp.close();
		ofp.clear();

		infn.close();
		infn.clear();

		infa.close();
		infa.clear();
	}

	nRec=0;


	inf.open("xyzVectorErec");
	if(inf)
	{
		inf>>nRec;
		inf.close();
	}
	inf.clear();

	if(nRec)
	{
		inf.open("recvsbe");
		if(!inf)
		{
			logfile<<"Error in open file "<<"recvsbe"<<endl;
			cout<<"Error in open file "<<"recvsbe"<<endl;
			return 1;
		}
		for(i=0;i<npls;i++){inf>>vNrec[i];}
		inf.close();
		inf.clear();

		for(ipls=0;ipls<npls;ipls++)
		{
			nall=vNrec[ipls];

			sprintf(str,"e2dstat.%d",ipls+1);
			infn.open(str);
			if(!infn)
			{
				logfile<<"Error in open file "<<str<<endl;
				cout<<"Error in open file "<<str<<endl;
				return 1;
			}
			sprintf(str,"e3dstat_anom.%d",ipls+1);
			infa.open(str);
			if(!infa)
			{
				logfile<<"Error in open file "<<str<<endl;
				cout<<"Error in open file "<<str<<endl;
				return 1;
			}
			sprintf(str,"e3dstat.%d",ipls+1);
			ofp.open(str);
			for(i=0;i<nall;i++)
			{
				infn>>fnrm;
				infa>>fanm;
				fsum=fnrm+fanm;
				ofp<<fsum<<'\t';

				infn>>fnrm;
				infa>>fanm;
				fsum=fnrm+fanm;
				ofp<<fsum<<'\t';

				infn>>fnrm;
				infa>>fanm;
				fsum=fnrm+fanm;
				ofp<<fsum<<'\n';
			}
			ofp.close();
			ofp.clear();

			infn.close();
			infn.clear();

			infa.close();
			infa.clear();
		}
	}

	close_logfile();

	return 0;
}
