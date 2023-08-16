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
 *  This file contains headers for time stamp functions
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 January 16, 2023
*/

#pragma once
#include<stdlib.h>
#include<time.h>
#include<vector>
#include<string>
using namespace std;


class TimeWatcher
{
	static TimeWatcher *TWInst;
	TimeWatcher();
	TimeWatcher(const TimeWatcher &twr);  
    TimeWatcher & operator = (TimeWatcher &twr);
public:
	static TimeWatcher & getInstance();
	void AddTimeSpot(const clock_t &ts,const char *descr);
	void AddMeshInfo(int nelem,int nodes,int nodes_term,int edges,int edges_term);
};
