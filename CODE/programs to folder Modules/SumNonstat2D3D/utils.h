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
 *  This file contains headers for utility functions
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 January 16, 2023
*/

#pragma once

extern char str[1024];

extern ofstream logfile;

extern char RunPath[1024];

void close_logfile();

struct profile_position
{
	double Ox,Oy,Oz,angle;
};

struct grpls
{
	int beg,end;
};

void MoveFileFromDirToDir(char *file_name,char *dir_source,char *dir_target);

void CopyFileFromDirToDir(char *file_name,char *dir_source,char *dir_target);

bool isFileExist(char *file_name);

void CloseProgrammWithError(int rcode);
