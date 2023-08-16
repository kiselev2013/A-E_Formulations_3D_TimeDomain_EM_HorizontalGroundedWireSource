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
 *  This file contains the code for starting calculation of the eleptic source part of primary field for nonstationary task
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
		

		logfile.open("log2dNonStatAx");


		if(argc>1)
		{
			CalculateAllDecades=1;
		}
		
		_mkdir("Ax");

		system("copy mtr3d2d Ax");
		system("copy sigma Ax");
		system("copy mu Ax");
		system("copy dpr Ax");
		system("copy infite.0 Ax");
		system("copy currentfunction Ax");
		system("copy deltafunction Ax");
		system("copy timeintervalforprint Ax");

		system("copy currentval3 Ax\\currentval");
		system("copy inf2tr.dat3 Ax\\inf2tr.dat");
		system("copy nvtr.dat3 Ax\\nvtr.dat");
		system("copy nvkat2d.dat3 Ax\\nvkat2d.dat");
		system("copy l1.dat3 Ax\\l1.dat");
		system("copy rz.dat3 Ax\\rz.dat");
		system("copy rz.txt3 Ax\\rz.txt");
		system("copy r.dat3 Ax\\r.dat");
		system("copy z.dat3 Ax\\z.dat");
		system("copy tsize.dat3 Ax\\tsize.dat");

		if (CalcSP()!=0)
		{
			logfile<<"Error!!!"<<endl;
			logfile.close();
			exit(1);
		}

	return 0;
}
