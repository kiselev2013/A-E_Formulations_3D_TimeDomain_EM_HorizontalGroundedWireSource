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
 *  This file contains template classes and functions for work with lists
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin                                
 *  Novosibirsk State Technical University,                        
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia              
 *  Corresponding author: vdv_wk@mail.ru                           
 *  Version 2.0 January 16, 2023                                   
*/

#pragma once

template <class Type> class Element {
  public:
  Type v;
  Element<Type> *next;
  Element() { next=NULL; };
};

template <class T> class ListOfValues {
  private:
    int n;
    Element<T> *list;
    Element<T> *last;
	mutable Element<T> *elcur;
  public:
    ListOfValues() 
	{ 
		n=0; 
		list=new Element<T>; 
		last=list; 
	}
    ~ListOfValues() 
	{ 
		Element<T> *cur=list;
		Element<T> *tmp;
		last=NULL;
		do
		{
			tmp=cur->next;
			delete cur;
            cur=tmp;
		}
		while (cur);
	}
    int GetListLength() const 
	{ 
		return n; 
	}
	bool AddToList(const T& v)
	{
		if ((last->next=new Element<T>)!=NULL) {
			last=last->next;
			last->v=v;
			n++;
			return true;
		}
		return false;
	}
	void LoadList(T *mhr) const
	{
		int i=0;
		Element<T> *cur=list->next;
		while (cur!=NULL) {
			mhr[i++]=cur->v;
			cur=cur->next;
		}
	}
	bool FindInList(const T& v) const
	{
		Element<T>* plist=list->next;
		while(plist)
		{
			if (plist->v==v)
				return true;
			plist=plist->next;
		}
		return false;
	}

	
	void ResetScan() const
	{
		elcur=list->next;
	}
	void MoveScan() const
	{
		if (elcur!=NULL)
			elcur=elcur->next;
	}
	const Element<T>* GetCurrentElPtr() const
	{
		return elcur;
	}
};

template <class T> class QueueOfValues {
  private:
    int n;
    Element<T> *list;
    Element<T> *last;
  public:
    QueueOfValues() { n=0; list=last=NULL; }
    ~QueueOfValues() { list=last=NULL; }
	bool EmptyQueue() { return (list==NULL); }
	int GetQueueLength() { return n; }
	bool AddToQueue(const T& v)
	{
		if (list==NULL) {
			if ((list=new Element<T>)!=NULL) {
				list->v=v;
				last=list;
				n++;
				return true;
			}
		}
		else {
			if ((last->next=new Element<T>)!=NULL) {
				last=last->next;
				last->v=v;
				n++;
				return true;
			}
		}
		return false;
	}
	T GetFromQueue()
	{
		T ElementT; ElementT=list->v; n--;
		if (list==last) 
			list=last=NULL;
		else 
			list=list->next;
		return ElementT;
	}
	void LoadQueue(T *mhr)
	{
		int i=0;
		Element<T> *cur=list->next;
		while (cur!=NULL) {
			mhr[i++]=cur->v;
			cur=cur->next;
		}
	}
};
