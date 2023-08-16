/**
 *  This code is freely available under the following conditions:
 *
 *  1) The code is to be used only for non-commercial purposes.
 *  2) No changes and modifications to the code without prior permission of the developer.
 *  3) No forwarding the code to a third party without prior permission of the developer.
 *
 *
 * The file contains basic utility functions and simple operations
 *
 *  Written by Ph.D. Dmitry S. Kiselev
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  harlequin_00@mail.ru
  *  Version 2.0 16 January, 2023
*/

#pragma once

#include <algorithm>
#include <map>
#include <set>
#include <vector>

#include "Mesh1D.h"
#include "Logging.h"


namespace Mesh3D
{
	using namespace std;

	// Common settings class
	struct Settings
	{
	public:
		double steps[3]; // Mesh initial steps
		double sparse[3]; // Sparce coefficient
		double farBound[3]; // Distance from receivers bounding box+gap to the boundary of the calculation domain
		double gapFromReceivers[3]; // Gap between receivers bounding box and sparced mesh area
		double eps = 0.1; // Gate for determining if coordinate lines are equal
		int receiversRefinementsCount; // Mesh refinements around receivers
	};


	// Get cell (with ix, iy, iz 1D mesh indeces) neighbours
	void GetNeighbours(int nx, int ny, int nz, int ix, int iy, int iz, int cell, int *neighbours);

	// Build base 3d mesh using model bounds
	int BuildBaseMesh(Settings &meshSettings, double **modelBounds, vector<double> *mesh);

	// Build template base 3d mesh using model bounds
	int BuildBaseMesh3DTemplate(double **modelBounds, vector<double> &layersBounds, vector<double> *mesh, vector<double> *meshTemplate);

	// Align 3D mesh begin and end with bounds 
	int AlignMeshBounds(double **modelBounds, vector<double> *mesh);

	// Find edge by nodes
	int FindEdge(int node1, int node2, map<pair<int, int>, int> &edgesToNumbers);

	// Colelct mesh edges
	int BuildEdges(vector<double> *mesh, set<pair<int, int>> &edges);

	// Colelct mesh elements described with edges
	int BuildElementsByEdges(vector<double> *mesh, set<pair<int, int>> &edges, vector<vector<int>> &elementsByEdges);

	// Get total edges
	int CalculateEdgesCount(vector<double> *mesh);
}