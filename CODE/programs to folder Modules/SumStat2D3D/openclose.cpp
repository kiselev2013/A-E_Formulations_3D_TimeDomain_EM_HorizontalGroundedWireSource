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
 *  This file contains code for open-close file functions
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 January 16, 2023
*/

#include "stdafx.h"
#include "openclose.h"
#include "utils.h"

extern ofstream logfile;

int OpenFile(ifstream &inf,char *fname,bool fstop)
{
	inf.open(fname);
	if(!inf)
	{
		inf.clear();
		sprintf(str,"Can't open file %s",fname);
		cout<<str<<endl;
		logfile<<str<<endl;
		if(fstop)
		{
			CloseProgrammWithError(1);
		}
		else
		{
			return 1;
		}
	}
	return 0;
}

void CloseFile(ifstream &inf,char *fname,bool fstop)
{
	if(!inf.good() && fstop)
	{
		sprintf(str,"Error while reaing file %s",fname);
		cout<<str<<endl;
		logfile<<str<<endl;
		CloseProgrammWithError(1);
	}
	inf.close();
	inf.clear();
}
