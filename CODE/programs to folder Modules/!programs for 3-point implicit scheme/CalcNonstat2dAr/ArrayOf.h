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
 *  This file contains the code for working with template dynamic arrays
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin                            
 *  Novosibirsk State Technical University,                    
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia          
 *  Corresponding author: vdv_wk@mail.ru                       
 *  Version 2.0 January 16, 2023                                                  
*/

#pragma once

template<class T> class ArrayOf {
private:
	bool mem; 
	void InitArray()
	{
		n=0; 
		val=NULL;
		mem=false;
	}
public:
	int n;
	T* val;
	ArrayOf() 
	{ 
		InitArray();
	}
	ArrayOf(int size) 
	{
		InitArray();
		mem=((val=new T[n=size])!=NULL);
	}
	ArrayOf(int size, T* mass)
	{
		InitArray();
		if (!mass) return;
		mem=((val=new T[n=size])!=NULL);
		if (!mem) return;
		for (int i=0; i<n; i++)
			val[i]=mass[i];
	}
	ArrayOf(int size, T value)
	{
		InitArray();
		mem=((val=new T[n=size])!=NULL);
		if (!mem) return;
		for (int i=0; i<n; i++)
			val[i]=value;
	}
	~ArrayOf()
	{
		if (mem) delete [] val;
	}
};

template<class T> class ResizableArrayOf 
{
private:
	int AllocatedSize, Size;
	bool mem;
public:
	T* val;
	ResizableArrayOf(int size)
	{
		val=NULL;
		mem=((val=new T[AllocatedSize=size])!=NULL);
		Size=0;
	}
	~ResizableArrayOf()
	{
		if (val) delete [] val;
	}
	int SetSize(const int& newsize)
	{
		if (newsize<0||newsize>AllocatedSize) 
			return -1;
		Size=newsize;
		return Size;
	}
	const int& GetAllocatedSize() const 
	{
		return AllocatedSize;
	}
	const int& GetSize() const 
	{
		return Size;
	}
	int Add(const T& t)
	{
		if (Size+1<AllocatedSize)
			val[Size++]=t;
		else return -1;
		return Size;
	}
	bool Find(const T& t)
	{
		for (int i=0; i<Size; i++)
			if (val[i]==t)
				return true;
		return false;
	}
};
/*! @} */