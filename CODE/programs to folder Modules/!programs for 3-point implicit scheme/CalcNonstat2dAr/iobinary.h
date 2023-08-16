/** \addtogroup Mesh_Generator */
/*! @{ */
#pragma once

#include <fstream>

namespace MESHG
{

using namespace std;

namespace BinIO 
{
	/*!     */
	template<class type>
	ifstream& operator>(ifstream& file, type &id)
	{
		file.read((char*)&id, sizeof(type));
		return file;
	}
	/*!     */
	template<class type>
	ofstream& operator<(ofstream& file, const type &id)
	{
		file.write((char*)&id, sizeof(type));
		return file;
	}

}

};
/*! @} */