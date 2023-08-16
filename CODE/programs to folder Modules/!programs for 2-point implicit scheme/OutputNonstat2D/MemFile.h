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
 *  This file contains headers for read-write mapping file functions
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 January 16, 2023
*/

#pragma once

wchar_t *convertCharArrayToLPCWSTR(const char* charArray);

int OpenMemFile(int MemSize,char *buf,HANDLE &hMapFile,LPCTSTR &pBuf);
int CloseMemFile(int MemSize,HANDLE &hMapFile,LPCTSTR &pBuf);

int WriteDataToMem(char *Buff,int MemSize,char *EName);
int WriteDataToMemConCat(char *Buff,int MemSize,char *EName,int AddOffset,char *AddBuff,int AddMemSize);
int ReadDataFromMem(char *Buff,int MemSize,char *EName);

int WriteDataToMem(float *DataF,int Nrec,char *MemName);
int ReadDataFromMem(float *DataF,int Nrec,char *MemName);

int WriteDataToMem(double *DataF,int Nrec,char *MemName);
int ReadDataFromMem(double *DataF,int Nrec,char *MemName);

int WriteDataToMem(int *DataF,int Nrec,char *MemName);
int ReadDataFromMem(int *DataF,int Nrec,char *MemName);

int WriteDataToMem(int &IntVal,char *MemName);
int ReadDataFromMem(int &IntVal,char *MemName);

struct MemFile
{
	HANDLE hMapFile;
	LPCTSTR pBuf;
	int MemSize;
	char *MemName;
	MemFile();
	void SetFileName(char *fname);
	int OpenMemFileC();
	int CloseMemFileC();
	virtual void GetMemmory()=0;
};

struct MemFileFloatMass : public MemFile
{
	int Nrec;
	float *DataF;
	MemFileFloatMass();
	void WriteDataToMemC();
	void ReadDataFromMemC();
	void GetMemmory();
	void CopyDataToVector(vector<float> &vec);
	void SetSize(int num);
};

struct MemFileDoubleMass : public MemFile
{
	int Nrec;
	double *DataF;
	MemFileDoubleMass();
	void WriteDataToMemC();
	void ReadDataFromMemC();
	void GetMemmory();
	void CopyDataToVector(vector<double> &vec);
	void SetSize(int num);
};

struct MemFileIntegerMass : public MemFile
{
	int Nrec;
	int *DataF;
	MemFileIntegerMass();
	void WriteDataToMemC();
	void ReadDataFromMemC();
	void GetMemmory();
	void CopyDataToVector(vector<int> &vec);
	void SetSize(int num);
};

struct MemFileIntegerValue : public MemFile
{
	int IntVal;
	MemFileIntegerValue();
	void WriteDataToMemC();
	void ReadDataFromMemC();
	void GetMemmory();
	void SetSize();
};
