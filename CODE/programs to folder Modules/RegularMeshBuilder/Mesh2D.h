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

#include <vector>

namespace Mesh2D
{
	using namespace std;

	struct Settings
	{
	public:
		double steps[2]; // Mesh initial steps
		double sparse[2]; // Sparce coefficient
		double farBound[2]; // Distance from receivers bounding box+gap to the boundary of the calculation domain
		double gapFromReceivers[2]; // Gap between receivers bounding box and sparced mesh area
		double eps = 1e-3; // Gate for determining if coordinate lines are equal
	};

	int BuildBaseMesh(Settings &meshSettings, double **modelBounds, vector<double> *mesh);

	// Align 2D mesh begin and end with bounds 
	int AlignMeshBounds(double **modelBounds, vector<double> *mesh);

	vector<int> GetFirstBoundaryConditionNodesNumbers(int nx, int ny, int index);
}