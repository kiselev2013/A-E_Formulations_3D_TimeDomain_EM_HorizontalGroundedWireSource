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
 *  This file contains the headers for application the main homogeneous boundary conditions in 3D VFEM
 *  
 *  
 *  Written by Ph.D. Denis V. Vagin                         
 *  Novosibirsk State Technical University,                 
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia       
 *  vdv_wk@mail.ru                                          
 *  Version 2.0 January 16, 2023                            
*/

// The class contains procedures for taking into account the main boundary conditions for a 3D VFEM problem
class T_Bound_Cond
{
public:
	double (*xyz)[3]; //  
	int (*edges)[2]; // ,  2- 
	int (*ed)[25];   // ,   
	int (*nver)[14]; // ,   
	int *bound_nodes; //  ,    
	int n_edges;  //   
	int n_nodes;  //   ()
	int n_bound_nodes; //  ,    
	int n_elem; //  
	double time; //   

	int n_bound_edges; //  ,    
	int *is_edge_bound; //      
	int *is_node_bound; //      
	double *dirichlet; //   1-   

	T_Bound_Cond(int n_nodes, int n_edges, int n_bound_nodes,
		int *bound_nodes, double (*xyz)[3], int (*edges)[2],
		int n_elem, int (*ed)[25], int (*nver)[14]);

	void Make_List_Of_Bound_Edges();
		
	void Set_Dirichlet_Cond(int *ig, int *jg, double *di, double *ggl, double *pr, 
		double time, int ntest);

	void Set_Dirichlet_Cond_Pr(double *pr, double time, int ntest, int ipls);

	void Make_Dirichlet();

	~T_Bound_Cond();
};
