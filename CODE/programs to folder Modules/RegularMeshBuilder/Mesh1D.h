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

#include "Logging.h"

namespace Mesh1D
{
	using namespace std;

	// Build 1D mesh within bounds with constant step
	int Build1DUniformMesh(double step, double *bounds, vector<double> &mesh);

	// Build 1D mesh within bounds with increasing by sparce step. Direction > 0 is from left to right, otherwise - from right to left
	int Build1DSparsedMesh(double step, double sparse, int direction, double *bounds, vector<double> &mesh);

	// Build base 1d mesh using model bounds
	int BuildBaseMesh1D(double step, double sparse, double *modelBounds, vector<double> &mesh);

	// Align 1D mesh begin and end with bounds 
	int AlignMeshBounds(double *bounds, vector<double> &mesh);

	// Build mesh for 1D problem
	int Build1DMesh(vector<double> &mesh1D);

	// Build mesh for 1D problem
	int Build1DMesh(vector<double> &mesh1D);

	// Refine mesh around points
	void RefineMesh(vector<double> &mesh1D, vector<double> &points, int receiversRefinementsCount, double precision);

	// Get interval between model bounds where coordinate is placed
	size_t GetLocalCoordinateInBoundsInterval(const double coordinate, const vector<double> &modelBounds);

	// Get local coordinate between model bounds where coordinate is placed
	double GetLocalCoordinateInBounds(const double coordinate, const vector<double> &modelBounds);
}