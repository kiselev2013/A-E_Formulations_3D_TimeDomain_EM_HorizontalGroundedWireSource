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
 *  This file contains code for utility functions
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

char str[1024];

ofstream logfile;

char RunPath[1024];

void close_logfile()
{
	logfile.close();
	logfile.clear();
}

void MoveFileFromDirToDir(char *file_name,char *dir_source,char *dir_target)
{
	sprintf(str,"move %s\\%s %s",dir_source,file_name,dir_target);
	system(str);
}

void CopyFileFromDirToDir(char *file_name,char *dir_source,char *dir_target)
{
	sprintf(str,"copy %s\\%s %s",dir_source,file_name,dir_target);
	system(str);
}

bool isFileExist(char *file_name)
{
	ifstream inf;
	bool flag;
	flag = false;
	inf.open(file_name);
	if (inf)
	{
		flag = true;
		inf.close();
	}
	inf.clear();
	return flag;
}

void CloseProgrammWithError(int rcode)
{
	ofstream ofp;
	_chdir(RunPath);
	ofp.open("stop");
	ofp << "1" << '\n';
	ofp.close();
	ofp.clear();
	close_logfile();
	exit(rcode);
}
