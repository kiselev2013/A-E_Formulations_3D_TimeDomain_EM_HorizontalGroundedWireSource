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

#include "Mesh1D.h"
#include "BasicOperations.h"
#include <stdexcept>
#include <algorithm>

namespace Mesh1D
{
	// Build 1D mesh within bounds with constant step
	int Build1DUniformMesh(double step, double *bounds, vector<double> &mesh)
	{
		try
		{
			if (bounds[1] < bounds[0]) { write_to_log("Error : Build1DUniformMesh : Mesh building bounds are inadequate\n"); return 1; }
			if (step < 0) { write_to_log("Error : Build1DUniformMesh : Step is less tha zero\n"); return 1; }
			int stepsCount = (bounds[1] - bounds[0]) / step + 0.5;

			if (fabs(step * stepsCount - (bounds[1] - bounds[0])) > 1e-5)
				stepsCount++;

			step = (bounds[1] - bounds[0]) / stepsCount;

			mesh.resize(stepsCount + 1);
			for (size_t i = 0; i < stepsCount + 1; i++)
				mesh[i] = bounds[0] + i * step;

			return 0;
		}
		catch (exception ex)
		{
			write_to_log("Error : Build1DUniformMesh : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Build 1D mesh within bounds with increasing by sparce step. Direction > 0 is from left to right, otherwise - from right to left
	int Build1DSparsedMesh(double step, double sparse, int direction, double *bounds, vector<double> &mesh)
	{
		try
		{
			if (bounds[1] < bounds[0]) { write_to_log("Error : Build1DSparsedMesh : Mesh building bounds are inadequate\n"); return 1; }
			if (step < 0) { write_to_log("Error : Build1DSparsedMesh : Step is less tha zero\n"); return 1; }
			if (sparse < 1.0) { write_to_log("Error : Build1DSparsedMesh : Sparse coefficient is less tha 1.0\n"); return 1; }

			double coord;
			if (direction > 0)
			{
				coord = bounds[0];
				mesh.push_back(coord);
				do {
					coord += step;
					step *= sparse;
					mesh.push_back(coord);
				} while (coord < bounds[1]);
			}
			else
			{
				vector<double> meshReverse;
				coord = bounds[1];
				meshReverse.push_back(coord);
				do {
					coord -= step;
					step *= sparse;
					meshReverse.push_back(coord);
				} while (coord > bounds[0]);

				mesh.resize(meshReverse.size());
				int meshIndex = 0;
				for (int meshReverseIndex = meshReverse.size() - 1; meshReverseIndex >= 0; meshReverseIndex--, meshIndex++)
					mesh[meshIndex] = meshReverse[meshReverseIndex];
			}

			return 0;
		}
		catch (exception ex)
		{
			write_to_log("Error : Build1DSparsedMesh : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Build base 1d mesh using model bounds
	int BuildBaseMesh1D(double step, double sparse, double *modelBounds, vector<double> &mesh)
	{
		try
		{
			vector<double> meshL, meshM, meshR;

			if (Build1DUniformMesh(step, modelBounds + 1, meshM) != 0) { write_to_log("Error : BuildBaseMesh1D : Could not build middle uniform mesh\n"); return 1; }
			if (Build1DSparsedMesh(step, sparse, -1, modelBounds + 0, meshL) != 0) { write_to_log("Error : BuildBaseMesh1D : Could not build left mesh\n"); return 1; }
			if (Build1DSparsedMesh(step, sparse, 1, modelBounds + 2, meshR) != 0) { write_to_log("Error : BuildBaseMesh1D : Could not build roght mesh\n"); return 1; }

			mesh.resize(meshL.size() + meshM.size() + meshR.size() - 2);

			int meshIndex = 0;
			for (int tmpIndex = 0; tmpIndex < meshL.size() - 1; tmpIndex++, meshIndex++)
				mesh[meshIndex] = meshL[tmpIndex];

			for (int tmpIndex = 0; tmpIndex < meshM.size() - 1; tmpIndex++, meshIndex++)
				mesh[meshIndex] = meshM[tmpIndex];

			for (int tmpIndex = 0; tmpIndex < meshR.size(); tmpIndex++, meshIndex++)
				mesh[meshIndex] = meshR[tmpIndex];

			return 0;
		}
		catch (exception ex)
		{
			write_to_log("Error : BuildBaseMesh1D : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Align 1D mesh begin and end with bounds 
	int AlignMeshBounds(double *bounds, vector<double> &mesh)
	{
		try
		{
			mesh[0] = bounds[0];
			if ((mesh[1] - mesh[0]) < (mesh[2] - mesh[1]) * 0.3)
				mesh.erase(mesh.begin() + 1);

			mesh[mesh.size() - 1] = bounds[3];
			if ((mesh[mesh.size() - 1] - mesh[mesh.size() - 2]) < (mesh[mesh.size() - 2] - mesh[mesh.size() - 3]) * 0.3)
				mesh.erase(mesh.begin() + mesh.size() - 2);

			return 0;
		}
		catch (exception ex)
		{
			write_to_log("Error : AlignMeshBounds : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Build mesh for 1D problem
	int Build1DMesh(vector<double> &mesh1D)
	{
		try
		{
			mesh1D.reserve(1000);

			double z = 0;
			double step = 1e-2;
			double sparse = 1.02;
			while (z > -1e+7)
			{
				mesh1D.push_back(z);
				z -= step;
				step *= sparse;
			}
			mesh1D.push_back(z);
			return 0;
		}
		catch (exception ex)
		{
			write_to_log("Error : Build1DMesh : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}


	pair<int, int> GetRefinementInterval(vector<double> &mesh1D, double point, int refinmentElementsCount)
	{
		int intervalIndex = std::upper_bound(mesh1D.begin(), mesh1D.end(), point) - mesh1D.begin() - 1;
		return {
			max(0, intervalIndex - refinmentElementsCount),
			min((int)mesh1D.size() - 1, intervalIndex + refinmentElementsCount + 1)
		};
	}

	void AddRefinementCoordinates(vector<double> &mesh1D, pair<int, int> &interval, vector<double> &coordinates)
	{
		for (int i = interval.first; i < interval.second; i++)
		{
			coordinates.push_back((mesh1D[i + 1] + mesh1D[i]) * 0.5);
		}
	}

	// Refine mesh around points
	void RefineMesh(vector<double> &mesh1D, vector<double> &points, int refinmentsCount, double precision)
	{
		vector<double> coordinates;

		int refinementIndex = 0;
		for (int refinementIndex = 0; refinementIndex < refinmentsCount; ++refinementIndex)
		{
			coordinates.clear();
			for (auto point : points)
				AddRefinementCoordinates(mesh1D, GetRefinementInterval(mesh1D, point, 1), coordinates);

			for (auto coordinate : coordinates)
				BasicOperations::InsertCoordinateIntoSortedArray(mesh1D, coordinate, precision);
		}
	}


	size_t GetLocalCoordinateInBoundsInterval(const double coordinate, const vector<double> &modelBounds)
	{
		const double eps = 1e-5;
		auto found = upper_bound(modelBounds.begin(), modelBounds.end(), coordinate);

		if (found == modelBounds.begin())
			return 0;

		if (found == modelBounds.end())
			return modelBounds.size() - 2;

		return min(static_cast<size_t>(distance(modelBounds.begin(), found)) - 1, modelBounds.size() - 2);
	}

	double GetLocalCoordinateInBounds(const double coordinate, const vector<double> &modelBounds)
	{
		size_t interval = GetLocalCoordinateInBoundsInterval(coordinate, modelBounds);
		const double bottom = modelBounds[interval];
		const double top = modelBounds[interval + 1];
		double zeta = (coordinate - bottom) / (top - bottom);
		return max(0.0, min(zeta + interval, static_cast<double>(modelBounds.size() - 1)));
	}
}