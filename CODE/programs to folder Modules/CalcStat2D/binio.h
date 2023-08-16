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
 *  This file contains functions for binary I/O
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin                        
 *  Novosibirsk State Technical University,                
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia      
 *  Corresponding author: vdv_wk@mail.ru                   
 *  Version 2.0 January 16, 2023                           
*/

#pragma once
#include <iostream>
#include <fstream>

using namespace std;

namespace binio
{

using std::ostream;
using std::istream;

__forceinline ostream& operator < (ostream& file,const double& data)
{
	file.write((char*)&data,sizeof(data));
	return file;
}

__forceinline istream& operator > (istream& file,double&  data)
{
	file.read((char*)&data,sizeof(data));
	return file;
}

__forceinline ostream& operator < (ostream& file,const float& data)
{
	file.write((char*)&data,sizeof(data));
	return file;
}

__forceinline istream& operator > (istream& file,float&  data)
{
	file.read((char*)&data,sizeof(data));
	return file;
}

__forceinline ostream& operator < (ostream& file,const int& data)
{
	file.write((char*)&data,sizeof(data));
	return file;
}

__forceinline istream& operator > (istream& file,int&  data)
{
	file.read((char*)&data,sizeof(data));
	return file;
}


__forceinline ostream& operator < (ostream& file,const bool& data)
{
	file.write((char*)&data,sizeof(data));
	return file;
}

__forceinline istream& operator > (istream& file,bool&  data)
{
	file.read((char*)&data,sizeof(data));
	return file;
}

__forceinline ostream& operator < (ostream& file,const char& data)
{
	file.write((char*)&data,sizeof(data));
	return file;
}

__forceinline istream& operator > (istream& file,char&  data)
{
	file.read((char*)&data,sizeof(data));
	return file;
}

__forceinline ostream& operator < (ostream& file,const unsigned char& data)
{
	file.write((char*)&data,sizeof(data));
	return file;
}

__forceinline istream& operator > (istream& file,unsigned char&  data)
{
	file.read((char*)&data,sizeof(data));
	return file;
}

__forceinline ostream& operator < (ostream& file,const short& data)
{
	file.write((char*)&data,sizeof(data));
	return file;
}

__forceinline istream& operator > (istream& file,short&  data)
{
	file.read((char*)&data,sizeof(data));
	return file;
}

__forceinline ostream& operator < (ostream& file,const unsigned int & data)
{
	file.write((char*)&data,sizeof(data));
	return file;
}

__forceinline istream& operator > (istream& file,unsigned int &  data)
{
	file.read((char*)&data,sizeof(data));
	return file;
}

}
