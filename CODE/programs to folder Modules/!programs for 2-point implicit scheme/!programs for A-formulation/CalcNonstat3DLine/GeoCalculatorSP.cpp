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
 *  This file contains main funvtion for calculating nonstationary 3D VFEM two-layer scheme tasks for A-formulation
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 January 16, 2023
*/

#include "stdafx.h"
#include "TaskSP.h"

extern int CalcSP();
extern ofstream logfile;

bool CheckStop(void)
{
	bool fstop;
	ifstream ifstop;
	fstop=false;
	ifstop.open("stop");
	if(ifstop){
		fstop=true;
		ifstop.close();
	}
	ifstop.clear();
	return fstop;
}

void Calc1AB(void)
{
	ifstream inf;


	logfile.open("log3dNonStat");

	if(CalcSP()!=0)
	{
		logfile<<"Error!!!"<<endl;
		logfile.close();
		return;
	}
}

int main(int argc, char **argv)
{
	Calc1AB();

	return 0;
}


