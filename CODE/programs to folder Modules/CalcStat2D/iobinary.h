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

#include <fstream>

namespace MESHG
{

using namespace std;

namespace BinIO 
{
	template<class type>
	ifstream& operator>(ifstream& file, type &id)
	{
		file.read((char*)&id, sizeof(type));
		return file;
	}
	template<class type>
	ofstream& operator<(ofstream& file, const type &id)
	{
		file.write((char*)&id, sizeof(type));
		return file;
	}

}

};
