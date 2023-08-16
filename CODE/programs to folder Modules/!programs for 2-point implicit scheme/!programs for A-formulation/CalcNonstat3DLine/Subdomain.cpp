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
 *  This file contains the code of functions for extracting smoothing subdomains and obtaining a solution in the elements belonging to them
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin                         
 *  Novosibirsk State Technical University,                 
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia       
 *  vdv_wk@mail.ru                                          
 *  Version 2.0 January 16, 2023                            
*/

#include "stdafx.h"

using namespace std;

#include "Subdomain.h"

#include "Hex_Local_Matrix.h"
#include "gauss3.h"

#define EpsComp 1e-6


extern ofstream logfile;

int Subdomain::subind=0;
int Subdomain::npls=0;
int Subdomain::ntimes=0;
int Subdomain::ftout=0;
int Subdomain::fdfct=0;

extern void CloseProgramm(int code);

struct double_eps
{
	double value;
	void operator = (double &a){value=a;}
	double get_value(){return value;}
};

bool operator == (double_eps &a,double_eps &b){return fabs(a.value-b.value)<EpsComp;}
bool operator > (double_eps &a,double_eps &b){return a.value>b.value+EpsComp;}
bool operator < (double_eps &a,double_eps &b){return a.value<b.value-EpsComp;}
bool operator >= (double_eps &a,double_eps &b){return (a>b || a==b);}
bool operator <= (double_eps &a,double_eps &b){return (a<b || a==b);}

double_eps operator + (double_eps &a,double_eps &b)
{
	double_eps t;
	t.value=a.value+b.value;
	return t;
}

double_eps operator - (double_eps &a,double_eps &b)
{
	double_eps t;
	t.value=a.value-b.value;
	return t;
}

double operator * (double_eps &a,double &b)
{
	return a.value*b;
}

double operator * (double_eps &a,double_eps &b)
{
	return a.value*b.value;
}

void Memory_allocation_error(const char *var, const char *func)
{
	string str;
	str = "MEMORY ALLOCATION ERROR for variable ";
	str = str + '\"' + var + '\"' + " in function " + '\"' + func + "\"\n";
	logfile << str << flush;
	cerr    << str << flush;
	cout << str << flush;
	throw logic_error(str);
}
void Cannot_open_file(const char *fname, const char *func)
{
	string str;
	str = "CANNOT OPEN FILE ";
	str = str + '\"' + fname + '\"' + " in function " + '\"' + func + "\"\n";
	logfile << str << flush;
	cerr    << str << flush;
	cout << str << flush;
	throw logic_error(str);
}
void Cannot_open_file_but_continue(const char *fname, const char *func)
{
	string str;
	str = "Cannot open file ";
	str = str + '\"' + fname + '\"' + " in function " + '\"' + func + "\"\n";
	logfile << str << flush;
	cerr    << str << flush;
}

int Sn[6][4]={{0,2,4,6},{0,1,4,5},{0,1,2,3},{1,3,5,7},{2,3,6,7},{4,5,6,7}};
Subdomain::Subdomain()
{
	nver = NULL;
	xyz = NULL;
	xyz_r = NULL;
	p = NULL;
	di = NULL;
	gg = NULL;
	pr = NULL;
	x = NULL;
	d = NULL;
	sg = NULL;

	subind++;
	my_ind=subind;
	nthreads=1;
}
Subdomain::~Subdomain()
{
	if (p)    {delete p; p = NULL;}
	if (pr)   {delete [] pr; pr = NULL;}
	if (x)   {delete [] x; x = NULL;}
	if (di)   {delete [] di; di = NULL;}
	if (gg)   {delete [] gg; gg = NULL;}
	
	if (xyz)  {delete [] xyz;  xyz  = NULL;}
	if (xyz_r){delete [] xyz_r;xyz_r= NULL;}
	if (nver) {delete [] nver; nver = NULL;}

	if (d)    {delete [] d; d = NULL;}
	if (sg)   {delete [] sg; sg = NULL;}
}
int Subdomain::Init(int material, int levelNeighbors, vector< vector<int> > &PointresForElem, AbstractFEM3D *TaskCalcMesh)
{
	int i, j, k, m, l, rr;
	int level;
	bool flag;
	int sz;
	vector<int> newElems, newElems2;
	vector<bool> isElemInSubdomain;
	vector<bool> isNodeInSubdomain;

	this->material = material;

	int emat;
	int nnoe=TaskCalcMesh->GetElementNodesNumber();
	int kpar=TaskCalcMesh->GetNumberOfElements();
	int kuzlov=TaskCalcMesh->GetNumberOfNodes();

	isElemInSubdomain.resize(kpar, false);
	isNodeInSubdomain.resize(kuzlov, false);

	vector<int> renumNodeFromOldToNew;
	renumNodeFromOldToNew.resize(kuzlov, -1);

	for (i=0; i<kpar; i++)
	{
		if (PointresForElem[i].size() > 0)
		{
			emat=TaskCalcMesh->GetElementMaterial(i);

			if (this->material != -1)
			{
				if (this->material != emat)
					continue;
			}

			isElemInSubdomain[i] = true;

			for (j=0; j<nnoe; j++)
			{
				k = TaskCalcMesh->GetNodeNumberOnElement(i,j);
				if (k >= 0) 
					isNodeInSubdomain[k] = true;
			}
		}
	}

	for (level=0; level < levelNeighbors; level++)
	{
		newElems.clear();
		newElems2.clear();

		for (i=0; i<kpar; i++)
		{
			if (isElemInSubdomain[i])
				continue;

			if (this->material != -1)
			{
				emat=TaskCalcMesh->GetElementMaterial(i);

				if (this->material != emat)
					continue;
			}

			flag = false;

			for (j=0; j<nnoe; j++)
			{
				k = TaskCalcMesh->GetNodeNumberOnElement(i,j);

				if (k >= 0)
				{
					if (isNodeInSubdomain[k])
					{
						flag = true;
						break;
					}
				}
			}

			if (flag)
				newElems.push_back(i);
		}

		std::sort(newElems.begin(), newElems.end());
		std::unique_copy(newElems.begin(), newElems.end(), back_inserter(newElems2));

		sz = (int)newElems2.size();
		for (i=0; i<sz; i++)
		{
			isElemInSubdomain[newElems2[i]] = true;

			for (j=0; j<nnoe; j++)
			{
				k = TaskCalcMesh->GetNodeNumberOnElement(newElems2[i],j);

				if (k >= 0)
					isNodeInSubdomain[k] = true;
			}
		}
	}

	this->n_elem = 0;
	for (i=0; i<kpar; i++)
	{
		if (isElemInSubdomain[i])
			this->n_elem++;
	}

	this->ValueInCenter.resize(this->n_elem);
	this->renumElemFromNewToOld.resize(this->n_elem, -1);

	this->n_nodes = 0;
	for (i=0; i<kuzlov; i++)
	{
		if (isNodeInSubdomain[i])
			this->n_nodes++;
	}

	this->renumNodeFromNewToOld.resize(this->n_nodes, -1);

	k = 0;
	for (i=0; i<kpar; i++)
	{
		if (isElemInSubdomain[i])
		{
			this->renumElemFromNewToOld[k] = i;
			k++;
		}
	}

	renumElemFromOldToNew.clear();
	renumElemFromOldToNew.resize(kpar,-1);

	for (i=0; i<this->n_elem; i++)renumElemFromOldToNew[renumElemFromNewToOld[i]]=i;
	
	if (this->xyz)  {delete [] this->xyz;  this->xyz = NULL;}
	this->xyz = new double[/*this->n_nodes*/kuzlov][3];

	if (this->xyz_r)  {delete [] this->xyz_r;  this->xyz_r = NULL;}
	this->xyz_r = new double[/*this->n_nodes*/kuzlov][3];

	if (this->nver) {delete [] this->nver; this->nver = NULL;}
	this->nver = new int[/*this->n_elem*/kpar][14];

	for (i=0; i</*this->n_elem*/kpar; i++){
		for (j=0; j<nnoe; j++){
			pv::Point3D TempPoint;
			k = TaskCalcMesh->GetNodeNumberOnElement(i,j);
			nver[i][j]=k;
			TempPoint=TaskCalcMesh->GetNode(k);
			this->xyz[k][0]=TempPoint.getx();
			this->xyz[k][1]=TempPoint.gety();
			this->xyz[k][2]=TempPoint.getz();
		}
	}


	this->n_nodes_c = this->n_nodes;

	k = 0;	
	j = this->n_nodes_c;
	for (i=0; i<kuzlov; i++)
	{
		if (isNodeInSubdomain[i])
		{
			this->renumNodeFromNewToOld[k] = i;
			renumNodeFromOldToNew[i] = k;
			k++;
		}
	}

	if (this->xyz_r)  {delete [] this->xyz_r;  this->xyz_r = NULL;}

	if (this->xyz)  {delete [] this->xyz;  this->xyz = NULL;}
	this->xyz = new double[this->n_nodes][3];

	for (i=0; i<this->n_nodes; i++)
	{
		pv::Point3D TempPoint;
		TempPoint=TaskCalcMesh->GetNode(this->renumNodeFromNewToOld[i]);
		this->xyz[i][0]=TempPoint.getx();
		this->xyz[i][1]=TempPoint.gety();
		this->xyz[i][2]=TempPoint.getz();
	}

	if (this->nver) {delete [] this->nver; this->nver = NULL;}
	this->nver = new int[this->n_elem][14];

	for (i=0; i<this->n_elem; i++){

		int elem = this->renumElemFromNewToOld[i];

		this->nver[i][13] = TaskCalcMesh->GetTypeOfElement(elem);

		for (j=0; j<8; j++){
			this->nver[i][j] = renumNodeFromOldToNew[TaskCalcMesh->GetNodeNumberOnElement(elem,j)];
		}
	}

	logfile<<"Resultant mesh: "<<this->n_nodes<<" nodes; "<<this->n_elem<<" elements; "<<this->n_nodes-this->n_nodes_c<<" terminal nodes;"<<'\n';
	cout<<"Resultant mesh: "<<this->n_nodes<<" nodes; "<<this->n_elem<<" elements; "<<this->n_nodes-this->n_nodes_c<<" terminal nodes;"<<'\n';

	this->p = new Portret(this->nver, this->n_elem, this->n_nodes, this->n_nodes_c);
	this->p->Gen_T_Portrait2();

	vLV.resize(n_elem);

	for (i=0; i<this->n_nodes; i++)
	{
		pv::Point3D TempPoint;
		TempPoint=TaskCalcMesh->GetNodeTrue(this->renumNodeFromNewToOld[i]);
		this->xyz[i][0]=TempPoint.getx();
		this->xyz[i][1]=TempPoint.gety();
		this->xyz[i][2]=TempPoint.getz();
	}

	AsmGlobalMatrix();

	return 0;
}

void Subdomain::AsmGlobalMatrix()
{
	int i, j, k, m, it, jt, i_mu, j_nu;
	int ii, jj; //  
	int ig_n_1 = p->ig[n_nodes_c];

	if ((pr = new double[n_nodes_c*npls*ntimes*ncmp]) == 0) Memory_allocation_error("pr", "AsmGlobalMatrix");
	if ((x = new double[n_nodes_c*npls*ntimes*ncmp]) == 0) Memory_allocation_error("x", "AsmGlobalMatrix");
	if ((di = new double[n_nodes_c]) == 0) Memory_allocation_error("di", "AsmGlobalMatrix");
	if ((gg = new double[ig_n_1]) == 0) Memory_allocation_error("gg", "AsmGlobalMatrix");

	for (i=0; i<n_nodes_c; i++)
	{
		di[i] = 0;
	}

	for (i=0; i<n_nodes_c*npls*ntimes*ncmp; i++)
	{
		pr[i] = 0;
		x[i] = 0;
	}

	for (i=0; i<ig_n_1; i++)
		gg[i] = 0;

	for(i=0; i<n_elem; i++)
	{
		Hex_Local_Matrix L(i, nver, xyz);
		L.CalcMassMatrix();

		for(j=0; j<8; j++)
		{
			vLV[i].g[j]=0.0;
			for(k=0; k<8; k++){vLV[i].g[j]+=L.b[j][k];}

			ii = nver[i][j];
			di[ii] += L.b[j][j];
			for(k=0; k<8; k++)
			{
				jj = nver[i][k];
				if(jj < ii) 
				{
					for(m=p->ig[ii]; m<=p->ig[ii+1]-1; m++)
					{
						if(p->jg[m] == jj)
							gg[m] += L.b[j][k];
					}
				}
			}// k
		}// j
	}// i

	prds.factorize(n_nodes_c,p->ig,p->jg,gg,di,nthreads);
}
void Subdomain::AsmGlobalVector(int ipls,int itime,int icmp)
{
	int i, j, it, i_mu, offset;
	int ii; //  
	int ig_n_1 = p->ig[n_nodes_c];

	offset=((icmp*ntimes+itime)*npls+ipls)*n_nodes_c;

	for (i=0; i<n_nodes_c; i++)
		pr[offset+i] = 0;

	for(i=0; i<n_elem; i++)
	{


		Local_Vector &lv=vLV[i];

		for(j=0; j<8; j++)
		{
			ii = nver[i][j];
			pr[offset+ii] += /*L.g8[j]*/lv.g[j]*ValueInCenter[i];
		}// j
	}// i
}
void Subdomain::CalcRightPartVect(T_Brick &L)
{
	int i,j;
	double BzOnElem;

	BzOnElem=ValueInCenter[L.num];

	Hex_Local_Matrix LM(L.num,nver,xyz);

	LM.Calc_local_matrix_b_for_parallelepiped();

	for(i=0;i<8;i++){
		L.g8[i]=0;
		for(j=0;j<8;j++){
			L.g8[i]+=LM.b[i][j];
		}
		L.g8[i]*=BzOnElem;
	}
}
