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
 *  This file contains code for read-write mapping file functions
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 January 16, 2023
*/

#include "stdafx.h"
#include "MemFile.h"

extern ofstream logfile;

wchar_t *convertCharArrayToLPCWSTR(const char* charArray)
{
	wchar_t* wString = new wchar_t[4096];
	MultiByteToWideChar(CP_ACP, 0, charArray, -1, wString, 4096);
	return wString;
}

int OpenMemFile(int MemSize,char *buf,HANDLE &hMapFile,LPCTSTR &pBuf)
{
	if(MemSize)
	{
		hMapFile = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, MemSize, convertCharArrayToLPCWSTR(buf));
		if (hMapFile == NULL || hMapFile == INVALID_HANDLE_VALUE)
		{
			printf("        (%d).\n", GetLastError());
			return 1;
		}
		pBuf = (LPTSTR)MapViewOfFile(hMapFile, FILE_MAP_ALL_ACCESS, 0, 0, MemSize);
		if (pBuf == NULL)
		{
			printf("     (%d).\n", GetLastError());
			return 1;
		}
	}
	return 0;
}

int CloseMemFile(int MemSize,HANDLE &hMapFile,LPCTSTR &pBuf)
{
	if(MemSize)
	{
		UnmapViewOfFile(pBuf);
		CloseHandle(hMapFile);
	}
	return 0;
}

int WriteDataToMem(char *Buff,int MemSize,char *EName)
{
	HANDLE hMapFile;
	LPCTSTR pBuf;
	if(MemSize)
	{
		hMapFile = OpenFileMapping(FILE_MAP_WRITE,TRUE,convertCharArrayToLPCWSTR(EName));
		if (hMapFile == NULL || hMapFile == INVALID_HANDLE_VALUE)
		{
			printf("        (%d).\n", GetLastError());
			return 1;
		}
		pBuf = (LPTSTR)MapViewOfFile(hMapFile, FILE_MAP_WRITE, 0, 0, MemSize);
		if (pBuf == NULL)
		{
			printf("     (%d).\n", GetLastError());
			return 1;
		}
		CopyMemory((PVOID)pBuf, (char *)Buff, MemSize);
		UnmapViewOfFile(pBuf);
		CloseHandle(hMapFile);
	}
	return 0;
}

int WriteDataToMemConCat(char *Buff,int MemSize,char *EName,int AddOffset,char *AddBuff,int AddMemSize)
{
	HANDLE hMapFile;
	LPCTSTR pBuf;
	if(MemSize)
	{
		hMapFile = OpenFileMapping(FILE_MAP_ALL_ACCESS,TRUE,convertCharArrayToLPCWSTR(EName));
		if (hMapFile == NULL || hMapFile == INVALID_HANDLE_VALUE)
		{
			printf("        (%d).\n", GetLastError());
			return 1;
		}
		pBuf = (LPTSTR)MapViewOfFile(hMapFile, FILE_MAP_ALL_ACCESS, 0, 0, MemSize);
		if (pBuf == NULL)
		{
			printf("     (%d).\n", GetLastError());
			return 1;
		}

		CopyMemory(Buff, (PVOID)pBuf, MemSize);

		CopyMemory(Buff+AddOffset, AddBuff, AddMemSize);

		CopyMemory((PVOID)pBuf, (char *)Buff, MemSize);

		UnmapViewOfFile(pBuf);
		CloseHandle(hMapFile);
	}
	return 0;
}

int ReadDataFromMem(char *Buff,int MemSize,char *EName)
{
	HANDLE hMapFile;
	LPCTSTR pBuf;
	if(MemSize)
	{
		hMapFile = OpenFileMapping(FILE_MAP_READ,TRUE,convertCharArrayToLPCWSTR(EName));
		if (hMapFile == NULL || hMapFile == INVALID_HANDLE_VALUE)
		{
			printf("        (%d).\n", GetLastError());
			return 1;
		}
		pBuf = (LPTSTR)MapViewOfFile(hMapFile, FILE_MAP_READ, 0, 0, MemSize);
		if (pBuf == NULL)
		{
			printf("     (%d).\n", GetLastError());
			return 1;
		}
		CopyMemory((char *)Buff, (PVOID)pBuf, MemSize);
		UnmapViewOfFile(pBuf);
		CloseHandle(hMapFile);
	}
	return 0;
}

int WriteDataToMem(float *DataF,int Nrec,char *MemName)
{
	int MemSize=Nrec*sizeof(float);
	char *Buff=(char *)DataF;
	return WriteDataToMem(Buff,MemSize,MemName);
}

int ReadDataFromMem(float *DataF,int Nrec,char *MemName)
{
	int MemSize=Nrec*sizeof(float);
	char *Buff=(char *)DataF;
	return ReadDataFromMem(Buff,MemSize,MemName);
}

int WriteDataToMem(double *DataF,int Nrec,char *MemName)
{
	int MemSize=Nrec*sizeof(double);
	char *Buff=(char *)DataF;
	return WriteDataToMem(Buff,MemSize,MemName);
}

int ReadDataFromMem(double *DataF,int Nrec,char *MemName)
{
	int MemSize=Nrec*sizeof(double);
	char *Buff=(char *)DataF;
	return ReadDataFromMem(Buff,MemSize,MemName);
}

int WriteDataToMem(int *DataF,int Nrec,char *MemName)
{
	int MemSize=Nrec*sizeof(int);
	char *Buff=(char *)DataF;
	return WriteDataToMem(Buff,MemSize,MemName);
}

int ReadDataFromMem(int *DataF,int Nrec,char *MemName)
{
	int MemSize=Nrec*sizeof(int);
	char *Buff=(char *)DataF;
	return ReadDataFromMem(Buff,MemSize,MemName);
}

int WriteDataToMem(int &IntVal,char *MemName)
{
	int MemSize=1*sizeof(int);
	char *Buff=(char *)(&IntVal);
	return WriteDataToMem(Buff,MemSize,MemName);
}

int ReadDataFromMem(int &IntVal,char *MemName)
{
	int MemSize=1*sizeof(int);
	char *Buff=(char *)(&IntVal);
	return ReadDataFromMem(Buff,MemSize,MemName);
}

MemFile::MemFile(){MemName=NULL;}
void MemFile::SetFileName(char *fname)
{
	if(MemName){delete [] MemName;MemName=NULL;}
	MemName = new char[strlen(fname)+1];
	strcpy(MemName,fname);
}
int MemFile::OpenMemFileC(){return OpenMemFile(MemSize,MemName,hMapFile,pBuf);}
int MemFile::CloseMemFileC(){return CloseMemFile(MemSize,hMapFile,pBuf);}

MemFileFloatMass::MemFileFloatMass(){Nrec=0;DataF=NULL;}
void MemFileFloatMass::WriteDataToMemC(){CopyMemory((PVOID)pBuf,(char *)DataF, MemSize);}
void MemFileFloatMass::ReadDataFromMemC(){CopyMemory((char *)DataF, (PVOID)pBuf, MemSize);}
void MemFileFloatMass::GetMemmory(){DataF=new float[Nrec];}
void MemFileFloatMass::CopyDataToVector(vector<float> &vec)
{
	int i,n;
	n=(int)vec.size();
	if(Nrec!=n){vec.resize(Nrec);}
	for(i=0;i<Nrec;i++){vec[i]=DataF[i];}
}
void MemFileFloatMass::SetSize(int num){Nrec=num;MemSize=num*sizeof(float);}

MemFileDoubleMass::MemFileDoubleMass(){Nrec=0;DataF=NULL;}
void MemFileDoubleMass::WriteDataToMemC(){CopyMemory((PVOID)pBuf,(char *)DataF, MemSize);}
void MemFileDoubleMass::ReadDataFromMemC(){CopyMemory((char *)DataF, (PVOID)pBuf, MemSize);}
void MemFileDoubleMass::GetMemmory(){DataF=new double[Nrec];}
void MemFileDoubleMass::CopyDataToVector(vector<double> &vec)
{
	int i,n;
	n=(int)vec.size();
	if(Nrec!=n){vec.resize(Nrec);}
	for(i=0;i<Nrec;i++){vec[i]=DataF[i];}
}
void MemFileDoubleMass::SetSize(int num){Nrec=num;MemSize=num*sizeof(double);}

MemFileIntegerMass::MemFileIntegerMass(){Nrec=0;DataF=NULL;}
void MemFileIntegerMass::WriteDataToMemC(){CopyMemory((PVOID)pBuf,(char *)DataF, MemSize);}
void MemFileIntegerMass::ReadDataFromMemC(){CopyMemory((char *)DataF, (PVOID)pBuf, MemSize);}
void MemFileIntegerMass::GetMemmory(){DataF=new int[Nrec];}
void MemFileIntegerMass::CopyDataToVector(vector<int> &vec)
{
	int i,n;
	n=(int)vec.size();
	if(Nrec!=n){vec.resize(Nrec);}
	for(i=0;i<Nrec;i++){vec[i]=DataF[i];}
}
void MemFileIntegerMass::SetSize(int num){Nrec=num;MemSize=num*sizeof(int);}

MemFileIntegerValue::MemFileIntegerValue(){IntVal=0;}
void MemFileIntegerValue::WriteDataToMemC(){CopyMemory((PVOID)pBuf,(char *)&IntVal, MemSize);}
void MemFileIntegerValue::ReadDataFromMemC(){CopyMemory((char *)&IntVal, (PVOID)pBuf, MemSize);}
void MemFileIntegerValue::GetMemmory(){}
void MemFileIntegerValue::SetSize(){MemSize=sizeof(int);}
