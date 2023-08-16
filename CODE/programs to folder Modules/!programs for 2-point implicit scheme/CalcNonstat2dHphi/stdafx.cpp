#include "stdafx.h"

bool isFileExists(const char *fname)
{
	ifstream inf;
	bool flag;
	flag=false;
	inf.open(fname);
	if(inf)
	{
		flag=true;
		inf.close();
	}
	inf.clear();
	return flag;
}
