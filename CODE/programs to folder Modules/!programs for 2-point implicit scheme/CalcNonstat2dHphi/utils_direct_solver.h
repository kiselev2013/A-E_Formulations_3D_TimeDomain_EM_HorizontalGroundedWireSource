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
 *  This file contains headers for utility functions for calculating with source grouping
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 January 16, 2023
*/

#pragma once

int GetNumberOfPlaces(int &npls,int &nplsa);
int FindTimeLayer(int nt,double *times,double val);
int read_time_mesh(int &ntime,double *(&time));
int ReadField2d(char *dpath,double *u,int n,int ind);
int ReadField2d(char *dpath,double *u,int n,char *fname);
int WriteField2d(char *dpath,double *u,int n,int ind);
int WriteField2d(char *dpath,double *u,int n,char *fname);
int ReadDecInfo(int &nDec,int &iDec,vector<int> &DecIg);
int FindTimeLayersForDec(int ntime,double *time,int iDec,vector<int> &DecIg,int &npre1,int &npre2,int &tbeg,int &tend);
void GenerateDecInfo(int ntime,double *time);
