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
 *  This file contains the code for starting calculation of field V for stationary task
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

void Calc1AB()
{
	logfile.open("LogCalcStat3D");
	try
	{
		CalcSP();
	}
	catch(...)
	{
		exit(1);
	}
	logfile.close();
	logfile.clear();
}

int main()
{
	setlocale( LC_ALL, "English" );
	Calc1AB();
	return 0;
}
