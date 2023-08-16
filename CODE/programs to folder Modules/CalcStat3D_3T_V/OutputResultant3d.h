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
 *  This file contains the headers for outputting with smoothing in 3D tasks
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin                         
 *  Novosibirsk State Technical University,                 
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia       
 *  vdv_wk@mail.ru                                          
 *  Version 2.0 January 16, 2023                            
*/

#pragma once
#include "OutputArbitrary.h"
#include "Subdomain.h"
#include "AbstractFEM3D.h"
#include "PointVector.h"

struct _Plane_{
	pv::Vector N;
	double D;
	void set_D(pv::Point3D t){D=-(N.x()*t.x()+N.y()*t.y()+N.z()*t.z());}
};

// The class contains modules for smoothing the solution of a three-dimensional problem
class OutputResultant3d
{
private:
	vector<Subdomain> sub;					//   
	
	int n_pointres;
	double (*pointres)[3];

	AbstractFEM3D *TaskCalcMesh;			//     

	int levelNeighbors;						//      
	
	vector<int> ElemForPoint;				//      -,   
	vector< vector<int> > PointsForElem;	//   -   , -   
	vector<long_double> PointresXsorted; 
	vector<long_double> PointresYsorted;
	vector<long_double> PointresZsorted;

	int InitSubdomains();
	void PointresSort();
	int FindPointsForElems();

public:
	OutputResultant3d(AbstractFEM3D *TaskCalcMesh,const Res3DValueType& r_type);
	~OutputResultant3d();

	int Prepare();
	int Output(int itime);
	int OutputDirect(int itime);

	Res3DValueType ValueType;

	void SetLevelNeighbors(int val);
};
