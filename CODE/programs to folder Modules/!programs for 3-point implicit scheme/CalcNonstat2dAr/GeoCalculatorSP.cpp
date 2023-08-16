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
 *  This file contains the code for starting calculation of the radial source part of primary field for nonstationary task
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

int main(int argc, char **argv)
{
	ifstream inf;
	ofstream ofp;
	

	logfile.open("log2dNonStatAr");

	if(argc>1)
	{
		CalculateAllDecades=1;
	}

	_mkdir("Ar");

	system("copy mtr3d2d Ar");
	system("copy sigma Ar");
	system("copy mu Ar");
	system("copy dpr Ar");
	system("copy infite.0 Ar");
	system("copy currentfunction Ar");
	system("copy deltafunction Ar");
	system("copy timeintervalforprint Ar");

	system("copy currentval2 Ar\\currentval");
	system("copy inf2tr.dat2 Ar\\inf2tr.dat");
	system("copy nvtr.dat2 Ar\\nvtr.dat");
	system("copy nvkat2d.dat2 Ar\\nvkat2d.dat");
	system("copy l1.dat2 Ar\\l1.dat");
	system("copy rz.dat2 Ar\\rz.dat");
	system("copy rz.txt2 Ar\\rz.txt");
	system("copy r.dat2 Ar\\r.dat");
	system("copy z.dat2 Ar\\z.dat");
	system("copy tsize.dat2 Ar\\tsize.dat");

	if (CalcSP()!=0)
	{
		logfile<<"Error!!!"<<endl;
		logfile.close();
		exit(1);
	}

	return 0;
}
