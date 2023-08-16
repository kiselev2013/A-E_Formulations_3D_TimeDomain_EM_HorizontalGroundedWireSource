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
#include "open_close.h"

int OpenInputFile(ifstream &inf,char *fname,ios_base::open_mode mode)
{
	inf.open(fname,mode);
	if(inf)
	{
		return 0;
	}
	else
	{
		inf.clear();
		cout<<"Can't open file "<<fname<<endl;
		return 1;
	}
}

int OpenOutputFile(ofstream &ofp,char *fname,ios_base::open_mode mode)
{
	ofp.open(fname,mode);
	if(ofp)
	{
		return 0;
	}
	else
	{
		ofp.clear();
		cout<<"Can't open file "<<fname<<endl;
		return 1;
	}
}

int CloseInputFile(ifstream &inf,char *fname)
{
	if(inf.good())
	{
		inf.close();
		inf.clear();
		return 0;
	}
	else
	{
		inf.clear();
		cout<<"Error while reading file "<<fname<<endl;
		return 1;
	}
}

int CloseOutputFile(ofstream &ofp,char *fname)
{
	if(ofp.good())
	{
		ofp.close();
		ofp.clear();
		return 0;
	}
	else
	{
		ofp.clear();
		cout<<"Error while writing file "<<fname<<endl;
		return 1;
	}
}
