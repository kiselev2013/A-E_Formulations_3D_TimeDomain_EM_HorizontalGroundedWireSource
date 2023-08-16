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
 *  This file contains the code of auxiliary utilities for outputting with smoothing
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin                         
 *  Novosibirsk State Technical University,                 
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia       
 *  vdv_wk@mail.ru                                          
 *  Version 2.0 January 16, 2023                            
*/

#include "stdafx.h"
#include "T_Brick.h"
#include "OutputArbitrary.h"


extern ofstream logfile;

Output3dArbitrary::Output3dArbitrary(int withSpline3d, int withSpline2d, int zeroPlane,
									 int kuzlov, int kpar, int n_pointres,
									 double (*pointres)[3], double (*xyz)[3])
{
	this->withSpline3d = withSpline3d;	
	this->withSpline2d = withSpline2d;	
	this->zeroPlane = zeroPlane;
	this->kuzlov = kuzlov;
	this->kpar = kpar;
	this->n_pointres = n_pointres;
	this->pointres = pointres;
	this->xyz = xyz;

	this->nver = NULL;
	this->nvtr = NULL;

	PointresSort();
}
Output3dArbitrary::~Output3dArbitrary()
{
}
OutputNode3d::OutputNode3d(int withSpline3d, int withSpline2d, int zeroPlane,
						   int kuzlov, int kpar, int n_pointres, double (*pointres)[3],
						   double (*xyz)[3], int (*nvtr)[8]) : Output3dArbitrary(withSpline3d,
						   withSpline2d, zeroPlane, kuzlov, kpar, n_pointres, pointres, xyz)
{
	this->nvtr = nvtr;
	
	FindElemForReceivers();
}
OutputNode3d::~OutputNode3d()
{
}
int Output3dArbitrary::GetGlobalVertex(int nElem, int nLocVertex)
{
	if (nvtr!=NULL)
	{
		return nvtr[nElem][nLocVertex];
	}
	else
	{
		return nver[nElem][nLocVertex];
	}
}
void Output3dArbitrary::PointresSort()
{
	vector<PointRes2> p;
	int i;

	p.resize(n_pointres);

	for (i=0; i<n_pointres; i++)
	{
		p[i].num = i;
		p[i].point[0] = pointres[i][0];
		p[i].point[1] = pointres[i][1];
		p[i].point[2] = pointres[i][2];
	}
	sort(p.begin(), p.end(), PointRes_less_x());

	PointresXsorted.resize(n_pointres);
	for (i=0; i<n_pointres; i++)
	{
		PointresXsorted[i].i = p[i].num;
		PointresXsorted[i].d = p[i].point[0];
	}
	
	for (i=0; i<n_pointres; i++)
	{
		p[i].num = i;
		p[i].point[0] = pointres[i][0];
		p[i].point[1] = pointres[i][1];
		p[i].point[2] = pointres[i][2];
	}
	sort(p.begin(), p.end(), PointRes_less_y());
	PointresYsorted.resize(n_pointres);
	for (i=0; i<n_pointres; i++)
	{
		PointresYsorted[i].i = p[i].num;
		PointresYsorted[i].d = p[i].point[1];
	}

	for (i=0; i<n_pointres; i++)
	{
		p[i].num = i;
		p[i].point[0] = pointres[i][0];
		p[i].point[1] = pointres[i][1];
		p[i].point[2] = pointres[i][2];
	}
	sort(p.begin(), p.end(), PointRes_less_z());
	PointresZsorted.resize(n_pointres);
	for (i=0; i<n_pointres; i++)
	{
		PointresZsorted[i].i = p[i].num;
		PointresZsorted[i].d = p[i].point[2];
	}
}
int Output3dArbitrary::FindElemForReceivers()
{
	int i, j, sz, t;

	vector<long_double>::iterator beg_x, beg_y, beg_z, end_x, end_y, end_z; 

	long_double coordMin[3], coordMax[3]; 

	vector<long_double> pntInBrick, pntInBrickTmp;
	vector<long_double> tmp_x, tmp_y, tmp_z;


	PointresForElem.resize(kpar);

	elemForPoint.resize(n_pointres);
	for(i=0; i<n_pointres; i++)
		elemForPoint[i] = -1;		

	for(i=0; i<kpar; i++)
	{
		for (j=0; j<3; j++)
		{
			coordMin[j].d = xyz[GetGlobalVertex(i,0)][j];
			coordMax[j].d = xyz[GetGlobalVertex(i,7)][j];
		}

		beg_x = lower_bound(PointresXsorted.begin(), PointresXsorted.end(), coordMin[0], Long_double_less());
		if(beg_x == PointresXsorted.end())
			continue;

		beg_y = lower_bound(PointresYsorted.begin(), PointresYsorted.end(), coordMin[1], Long_double_less());
		if(beg_y == PointresYsorted.end())
			continue;

		beg_z = lower_bound(PointresZsorted.begin(), PointresZsorted.end(), coordMin[2], Long_double_less());
		if(beg_z == PointresZsorted.end())
			continue;

		end_x = upper_bound(PointresXsorted.begin(), PointresXsorted.end(), coordMax[0], Long_double_less());
		if(end_x == PointresXsorted.begin())
			continue;

		end_y = upper_bound(PointresYsorted.begin(), PointresYsorted.end(), coordMax[1], Long_double_less());
		if(end_y == PointresYsorted.begin())
			continue;

		end_z = upper_bound(PointresZsorted.begin(), PointresZsorted.end(), coordMax[2], Long_double_less());
		if(end_z == PointresZsorted.begin())
			continue;

		if(distance(beg_x, end_x) <= 0)
			continue;

		if(distance(beg_y, end_y) <= 0)
			continue;

		if(distance(beg_z, end_z) <= 0)
			continue;

		tmp_x.clear();
		copy(beg_x, end_x, back_inserter(tmp_x));
		sort(tmp_x.begin(), tmp_x.end(), Long_double_num_less());

		tmp_y.clear();
		copy(beg_y, end_y, back_inserter(tmp_y));
		sort(tmp_y.begin(), tmp_y.end(), Long_double_num_less());

		tmp_z.clear();
		copy(beg_z, end_z, back_inserter(tmp_z));
		sort(tmp_z.begin(), tmp_z.end(), Long_double_num_less());

		set_intersection(tmp_x.begin(), tmp_x.end(), tmp_y.begin(), tmp_y.end(), back_inserter(pntInBrickTmp), Long_double_num_less()); 
		set_intersection(tmp_z.begin(), tmp_z.end(), pntInBrickTmp.begin(), pntInBrickTmp.end(), back_inserter(pntInBrick), Long_double_num_less()); 

		sz = (int)pntInBrick.size();
		for (j=0; j<sz; j++)
		{
			t = pntInBrick[j].i;

			if (elemForPoint[t] == -1)
			{
				elemForPoint[t] = i;
				PointresForElem[i].push_back(t);
			}
		}

		pntInBrick.clear();
		pntInBrickTmp.clear();
	}

	for (i=0; i<n_pointres; i++)
	{
		if (elemForPoint[i]==-1)
		{
			char s[256];
			sprintf(s,"There is no element for receiver N %d",i);
			logfile<<s<<endl;
		}
	}
	
	return 0;
}
