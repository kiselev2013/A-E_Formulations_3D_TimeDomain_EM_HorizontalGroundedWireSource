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
 *  This file contains code for time stamp functions
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 January 16, 2023
*/

#include"stdafx.h"
#include"TimeWatcher.h"


TimeWatcher * TimeWatcher :: TWInst=NULL;

TimeWatcher :: TimeWatcher()
{
}

TimeWatcher & TimeWatcher :: getInstance()
{
   static TimeWatcher inst;
   return inst;
}

void TimeWatcher :: AddTimeSpot(const clock_t &t,const char *descr)
{
	ofstream ofp;
	ofp.open("TimeWatcher.log",ios::app);
	ofp<<(int)t<<"\t"<<descr<<'\n';
	ofp.close();
	ofp.clear();
}



void TimeWatcher :: AddMeshInfo(int nelem,int nodes,int nodes_term,int edges,int edges_term)
{
	ofstream ofp;
	ofp.open("TimeWatcher.log",ios::app);
	ofp<<"nelem= "<<nelem<<'\n';
	ofp<<"nodes= "<<nodes<<'\n';
	ofp<<"nodes_term= "<<nodes_term<<'\n';
	ofp<<"edges= "<<edges<<'\n';
	ofp<<"edges_term= "<<edges_term<<'\n';
	ofp.close();
	ofp.clear();
}
