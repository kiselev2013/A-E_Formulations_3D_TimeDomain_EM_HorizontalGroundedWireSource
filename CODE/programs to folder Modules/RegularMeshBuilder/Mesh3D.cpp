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

#include "Mesh3D.h"

#include "Mesh1D.h"

namespace Mesh3D
{
	// Get cell (with ix, iy, iz 1D mesh indeces) neighbours
	void GetNeighbours(int nx, int ny, int nz, int ix, int iy, int iz, int cell, int *neighbours)
	{
		int nxy = (nx - 1) * (ny - 1);
		neighbours[0] = ix == 0 ? -1 : 1 + cell - 1;
		neighbours[3] = ix == nx - 2 ? -1 : 1 + cell + 1;
		neighbours[1] = iy == 0 ? -1 : 1 + cell - (nx - 1);
		neighbours[4] = iy == ny - 2 ? -1 : 1 + cell + (nx - 1);
		neighbours[2] = iz == 0 ? -1 : 1 + cell - nxy;
		neighbours[5] = iz == nz - 2 ? -1 : 1 + cell + nxy;
	}

	// Build base 3d mesh using model bounds
	int BuildBaseMesh(Settings &meshSettings, double **modelBounds, vector<double> *mesh)
	{
		try
		{
			int direction;
			direction = 0; if (Mesh1D::BuildBaseMesh1D(meshSettings.steps[direction], meshSettings.sparse[direction], modelBounds[direction], mesh[direction]) != 0) { write_to_log("Error : BuildBaseMesh3D : Could not build X base mesh\n"); return 1; }
			direction = 1; if (Mesh1D::BuildBaseMesh1D(meshSettings.steps[direction], meshSettings.sparse[direction], modelBounds[direction], mesh[direction]) != 0) { write_to_log("Error : BuildBaseMesh3D : Could not build Y base mesh\n"); return 1; }
			direction = 2; if (Mesh1D::BuildBaseMesh1D(meshSettings.steps[direction], meshSettings.sparse[direction], modelBounds[direction], mesh[direction]) != 0) { write_to_log("Error : BuildBaseMesh3D : Could not build Z base mesh\n"); return 1; }
			return 0;
		}
		catch (exception ex)
		{
			write_to_log("Error : BuildBaseMesh3D : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Build template base 3d mesh using model bounds
	int BuildBaseMesh3DTemplate(double **modelBounds, vector<double> &layersBounds, vector<double> *mesh, vector<double> *meshTemplate)
	{
		try
		{
			vector<double> bounds[3];
			for (size_t i = 0; i < 3; i++)
				for (size_t j = 0; j < 4; j++)
					bounds[i].push_back(modelBounds[i][j]);
			bounds[2] = layersBounds;

			for (int direction = 0; direction < 3; direction++)
			{
				auto &mesh1D = mesh[direction];
				auto &mesh1DTemplate = meshTemplate[direction];
				mesh1DTemplate.resize(mesh1D.size());
				for (int meshIndex = 0; meshIndex < mesh1D.size(); meshIndex++)
				{
					double coordinate = mesh1D[meshIndex];
					mesh1DTemplate[meshIndex] = Mesh1D::GetLocalCoordinateInBounds(coordinate, bounds[direction]);
				}
			}

			return 0;
		}
		catch (exception ex)
		{
			write_to_log("Error : BuildBaseMesh3DTemplate : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Align 3D mesh begin and end with bounds 
	int AlignMeshBounds(double **modelBounds, vector<double> *mesh)
	{
		try
		{
			int direction;
			direction = 0; if (Mesh1D::AlignMeshBounds(modelBounds[direction], mesh[direction]) != 0) { write_to_log("Error : AlignMeshBounds : Could not align X mesh\n"); return 1; }
			direction = 1; if (Mesh1D::AlignMeshBounds(modelBounds[direction], mesh[direction]) != 0) { write_to_log("Error : AlignMeshBounds : Could not align Y mesh\n"); return 1; }
			direction = 2; if (Mesh1D::AlignMeshBounds(modelBounds[direction], mesh[direction]) != 0) { write_to_log("Error : AlignMeshBounds : Could not align Z mesh\n"); return 1; }
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

	// Find edge by nodes
	int FindEdge(int node1, int node2, map<pair<int, int>, int> &edgesToNumbers)
	{
		static pair<int, int> edge;
		edge.first = node1;
		edge.second = node2;

		auto &found = edgesToNumbers.find(edge);
		if (found == edgesToNumbers.end())
		{
			char buf[2048];
			sprintf(buf, "Error : Could not find edge for nodes %d and %d\n", node1, node2);
			write_to_log(buf);
			exit(1);
		}
		else
			return found->second;
	}

	// Colelct mesh edges
	int BuildEdges(vector<double> *mesh, set<pair<int, int>> &edges)
	{
		try
		{
			pair<int, int> edge;
			int nx, ny, nz;
			nx = mesh[0].size();
			ny = mesh[1].size();
			nz = mesh[2].size();
			int nxy = nx * ny;
			nx = mesh[0].size();
			ny = mesh[1].size();
			nz = mesh[2].size();

			int cell = 0;
			int node;
			for (int iz = 0; iz < nz - 1; iz++)
			{
				for (int iy = 0; iy < ny - 1; iy++)
				{
					node = iz * nxy + iy * nx + 1;
					for (int ix = 0; ix < nx - 1; ix++, node++)
					{
						edges.insert(pair<int, int>(node, node + 1));
						edges.insert(pair<int, int>(node + nx, node + nx + 1));
						edges.insert(pair<int, int>(node + nxy, node + nxy + 1));
						edges.insert(pair<int, int>(node + nxy + nx, node + nxy + nx + 1));

						edges.insert(pair<int, int>(node, node + nx));
						edges.insert(pair<int, int>(node + 1, node + nx + 1));
						edges.insert(pair<int, int>(node + nxy, node + nxy + nx));
						edges.insert(pair<int, int>(node + nxy + 1, node + nxy + nx + 1));

						edges.insert(pair<int, int>(node, node + nxy));
						edges.insert(pair<int, int>(node + 1, node + 1 + nxy));
						edges.insert(pair<int, int>(node + nx, node + nx + nxy));
						edges.insert(pair<int, int>(node + nx + 1, node + nx + 1 + nxy));
					}
				}
			}
			return 0;
		}
		catch (exception ex)
		{
			write_to_log("Error : SetObjectsMaterials : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Colelct mesh elements described with edges
	int BuildElementsByEdges(vector<double> *mesh, set<pair<int, int>> &edges, vector<vector<int>> &elementsByEdges)
	{
		try
		{
			pair<int, int> edge;
			int nx, ny, nz;
			nx = mesh[0].size();
			ny = mesh[1].size();
			nz = mesh[2].size();
			int nxy = nx * ny;
			nx = mesh[0].size();
			ny = mesh[1].size();
			nz = mesh[2].size();
			elementsByEdges.resize((nx - 1) * (ny - 1) * (nz - 1));
			map<pair<int, int>, int> edgesToNumbers;
			int number = 1;
			for (auto edge = edges.begin(); edge != edges.end(); ++edge, number++)
				edgesToNumbers[*edge] = number;

			int cell = 0;
			int node;
			for (int iz = 0; iz < nz - 1; iz++)
			{
				for (int iy = 0; iy < ny - 1; iy++)
				{
					node = iz * nxy + iy * nx + 1;
					for (int ix = 0; ix < nx - 1; ix++, node++, cell++)
					{
						elementsByEdges[cell].resize(12);

						elementsByEdges[cell][0] = FindEdge(node, node + 1, edgesToNumbers);
						elementsByEdges[cell][1] = FindEdge(node + nx, node + nx + 1, edgesToNumbers);
						elementsByEdges[cell][2] = FindEdge(node + nxy, node + nxy + 1, edgesToNumbers);
						elementsByEdges[cell][3] = FindEdge(node + nxy + nx, node + nxy + nx + 1, edgesToNumbers);

						elementsByEdges[cell][4] = FindEdge(node, node + nx, edgesToNumbers);
						elementsByEdges[cell][5] = FindEdge(node + nxy, node + nxy + nx, edgesToNumbers);
						elementsByEdges[cell][6] = FindEdge(node + 1, node + nx + 1, edgesToNumbers);
						elementsByEdges[cell][7] = FindEdge(node + nxy + 1, node + nxy + nx + 1, edgesToNumbers);

						elementsByEdges[cell][8] = FindEdge(node, node + nxy, edgesToNumbers);
						elementsByEdges[cell][9] = FindEdge(node + 1, node + 1 + nxy, edgesToNumbers);
						elementsByEdges[cell][10] = FindEdge(node + nx, node + nx + nxy, edgesToNumbers);
						elementsByEdges[cell][11] = FindEdge(node + nx + 1, node + nx + 1 + nxy, edgesToNumbers);
					}
				}
			}
			return 0;
		}
		catch (exception ex)
		{
			write_to_log("Error : SetObjectsMaterials : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Get total edges
	int CalculateEdgesCount(vector<double> *mesh)
	{
		try
		{
			int nx, ny, nz;
			int edgesCount = 0;
			nx = mesh[0].size();
			ny = mesh[1].size();
			nz = mesh[2].size();
			edgesCount += (nx - 1) * ny * nz;
			edgesCount += (ny - 1) * nx * nz;
			edgesCount += (nz - 1) * nx * ny;
			return edgesCount;
		}
		catch (exception ex)
		{
			write_to_log("Error : WriteTSize3DNode : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 0;
		}
	}
}