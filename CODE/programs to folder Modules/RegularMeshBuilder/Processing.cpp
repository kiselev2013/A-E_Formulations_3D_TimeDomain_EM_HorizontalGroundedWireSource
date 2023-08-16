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

#include "Processing.h"

#include <chrono>
#include <filesystem>

#include "BasicOperations.h"
#include "Model.h"
#include "Mesh2D.h"
#include "Mesh3D.h"
#include "MeshTime.h"
#include "ReadWrite.h"

namespace Processing
{

	// Create additional mesh lines for objects
	int InsertObjectsIntoMesh(vector<Model::GeoBox<double>> &objects, double eps, vector<double> *mesh)
	{
		try
		{
			char buf[2048];
			int objectNumber = 1;
			for (auto obj = objects.begin(); obj != objects.end(); ++obj, objectNumber++)
			{
				for (int bound = 0; bound < 6; bound++)
				{
					if (BasicOperations::InsertCoordinateIntoSortedArray(mesh[bound / 2], obj->coordinates[bound], eps) != 0)
					{
						sprintf(buf, "Error : InsertObjectsIntoMesh : Could not insert object %d bound %d into mesh\n", objectNumber, bound + 1); write_to_log(buf);
						write_to_log(buf);
						return 1;
					}
				}
			}
			return 0;
		}
		catch (exception ex)
		{
			write_to_log("Error : InsertObjectsIntoMesh : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Create additional mesh lines for layers
	int InsertLayersIntoMesh(vector<Model::GeoLayer<double>> &layers, double eps, vector<double> &meshZ)
	{
		try
		{
			char buf[2048];
			int layerNumber = 1;
			for (auto layer = layers.begin(); layer != layers.end(); ++layer, layerNumber++)
			{
				if (layerNumber == layers.size())
					break;

				if (BasicOperations::InsertCoordinateIntoSortedArray(meshZ, layer->bottom, eps) != 0)
				{
					sprintf(buf, "Error : InsertLayersIntoMesh : Could not insert object %d bound %d into mesh\n", layerNumber, layer->bottom); write_to_log(buf);
					write_to_log(buf);
					return 1;
				}
			}
			return 0;
		}
		catch (exception ex)
		{
			write_to_log("Error : InsertLayersIntoMesh : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Calculate receivers bounding box with rounding
	int CalculateModelInternalBounds(vector<Model::Receiver> &receivers, vector<Model::Generator> &generators, double *gaps, double *steps, Model::GeoBox<double> &bounds)
	{
		try
		{
			if (receivers.size() == 0 && generators.size() == 0) { write_to_log("Error : CalculateModelInternalBounds : No receivers and generators\n"); return 1; }

			bounds.coordinates[0] = bounds.coordinates[2] = bounds.coordinates[4] = DBL_MAX;
			bounds.coordinates[1] = bounds.coordinates[3] = bounds.coordinates[5] = -DBL_MAX;

			for (auto receiver = receivers.begin(); receiver != receivers.end(); ++receiver)
			{
				for (int direction = 0; direction < 3; direction++)
				{
					if (bounds.coordinates[direction * 2 + 0] > receiver->coordinates[direction]) bounds.coordinates[direction * 2 + 0] = receiver->coordinates[direction];
					if (bounds.coordinates[direction * 2 + 1] < receiver->coordinates[direction]) bounds.coordinates[direction * 2 + 1] = receiver->coordinates[direction];
				}
			}

			for (auto generator = generators.begin(); generator != generators.end(); ++generator)
			{
				for (int direction = 0; direction < 3; direction++)
				{
					if (bounds.coordinates[direction * 2 + 0] > generator->coordinates[direction]) bounds.coordinates[direction * 2 + 0] = generator->coordinates[direction];
					if (bounds.coordinates[direction * 2 + 1] < generator->coordinates[direction]) bounds.coordinates[direction * 2 + 1] = generator->coordinates[direction];
				}
			}

			for (int direction = 0; direction < 3; direction++)
			{
				bounds.coordinates[direction * 2 + 0] -= gaps[direction];
				bounds.coordinates[direction * 2 + 1] += gaps[direction];

				double d = bounds.coordinates[direction * 2 + 1] - bounds.coordinates[direction * 2 + 0];
				int stepsCount = (int)(d / steps[direction] + 0.5);
				if (fabs(stepsCount * steps[direction] - d) > 1e-7)
				{
					stepsCount++;
					bounds.coordinates[direction * 2 + 0] = BasicOperations::RoundDouble_(bounds.coordinates[direction * 2 + 0], steps[direction]);
					bounds.coordinates[direction * 2 + 1] = bounds.coordinates[direction * 2 + 0] + stepsCount * steps[direction];
				}
			}

			return 0;
		}
		catch (exception ex)
		{
			write_to_log("Error : CalculateModelInternalBounds : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}

	}

	// Calculate receivers bounding box
	int CalculateModelInternalBounds(vector<Model::Receiver> &receivers, vector<Model::Generator> &generators, Model::GeoBox<double> &bounds)
	{
		try
		{
			if (receivers.size() == 0 && generators.size() == 0) { write_to_log("Error : CalculateModelInternalBounds : No receivers and generators\n"); return 1; }

			bounds.coordinates[0] = bounds.coordinates[2] = bounds.coordinates[4] = DBL_MAX;
			bounds.coordinates[1] = bounds.coordinates[3] = bounds.coordinates[5] = -DBL_MAX;

			for (auto receiver = receivers.begin(); receiver != receivers.end(); ++receiver)
			{
				for (int direction = 0; direction < 3; direction++)
				{
					if (bounds.coordinates[direction * 2 + 0] > receiver->coordinates[direction]) bounds.coordinates[direction * 2 + 0] = receiver->coordinates[direction];
					if (bounds.coordinates[direction * 2 + 1] < receiver->coordinates[direction]) bounds.coordinates[direction * 2 + 1] = receiver->coordinates[direction];
				}
			}

			for (auto generator = generators.begin(); generator != generators.end(); ++generator)
			{
				for (int direction = 0; direction < 3; direction++)
				{
					if (bounds.coordinates[direction * 2 + 0] > generator->coordinates[direction]) bounds.coordinates[direction * 2 + 0] = generator->coordinates[direction];
					if (bounds.coordinates[direction * 2 + 1] < generator->coordinates[direction]) bounds.coordinates[direction * 2 + 1] = generator->coordinates[direction];
				}
			}

			return 0;
		}
		catch (exception ex)
		{
			write_to_log("Error : CalculateModelInternalBounds : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}

	}

	// Calculatemodel bounding box
	int CalculateModelBounds3D(vector<Model::Receiver> &receivers, vector<Model::Generator> &generators, vector<Model::GeoLayer<double>> &layers, Mesh3D::Settings &meshSettings, double **bounds)
	{
		try
		{
			Model::GeoBox<double> internalBoundsNotRounded;
			Model::GeoBox<double> internalBounds;
			if (CalculateModelInternalBounds(receivers, generators, internalBoundsNotRounded) != 0) { write_to_log("Error : CalculateModelBounds3D : Could not calculate receiver bounds\n"); return 1; }
			if (CalculateModelInternalBounds(receivers, generators, meshSettings.gapFromReceivers, meshSettings.steps, internalBounds) != 0) { write_to_log("Error : CalculateModelBounds3D : Could not calculate receiver bounds\n"); return 1; }

			for (int direction = 0; direction < 3; direction++)
			{
				bounds[direction][0] = internalBoundsNotRounded.coordinates[direction * 2] - meshSettings.farBound[direction];
				bounds[direction][1] = internalBounds.coordinates[direction * 2];
				bounds[direction][2] = internalBounds.coordinates[direction * 2 + 1];
				bounds[direction][3] = internalBoundsNotRounded.coordinates[direction * 2 + 1] + meshSettings.farBound[direction];
			}

			if (layers[0].top < bounds[2][0])
				bounds[2][0] = layers[0].top - (layers[1].top - layers[1].bottom);


			return 0;
		}
		catch (exception ex)
		{
			write_to_log("Error : CalculateModelBounds3D : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Calculatemodel bounding box
	int CalculateModelBounds2D(vector<Model::Receiver> &receivers, vector<Model::Generator> &generators, Mesh2D::Settings &meshSettings, double **bounds)
	{
		try
		{
			double steps[3]{ meshSettings.steps[0], meshSettings.steps[0], meshSettings.steps[1] };
			double gaps[3]{ meshSettings.gapFromReceivers[0], meshSettings.gapFromReceivers[0], meshSettings.gapFromReceivers[1] };
			Model::GeoBox<double> internalBoundsNotRounded;
			Model::GeoBox<double> internalBounds;
			if (CalculateModelInternalBounds(vector<Model::Receiver>(), generators, internalBoundsNotRounded) != 0) { write_to_log("Error : CalculateModelBounds2D : Could not calculate receiver bounds\n"); return 1; }
			if (CalculateModelInternalBounds(vector<Model::Receiver>(), generators, gaps, steps, internalBounds) != 0) { write_to_log("Error : CalculateModelBounds2D : Could not calculate receiver bounds\n"); return 1; }

			for (int direction = 0; direction < 3; direction += 2)
			{
				bounds[direction == 0 ? 0 : 1][0] = internalBoundsNotRounded.coordinates[direction * 2] - meshSettings.farBound[direction == 0 ? 0 : 1];
				bounds[direction == 0 ? 0 : 1][1] = internalBounds.coordinates[direction * 2];
				bounds[direction == 0 ? 0 : 1][2] = internalBounds.coordinates[direction * 2 + 1];
				bounds[direction == 0 ? 0 : 1][3] = internalBoundsNotRounded.coordinates[direction * 2 + 1] + meshSettings.farBound[direction == 0 ? 0 : 1];
			}


			return 0;
		}
		catch (exception ex)
		{
			write_to_log("Error : CalculateModelBounds2D : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Apply layers materials to mesh
	int SetLayersMaterials3D(vector<Model::GeoLayer<double>> &layers, vector<Model::Material> &layerMaterials, vector<double> *mesh)
	{
		try
		{
			int topIndex, bottomIndex, cell;
			int nx = mesh[0].size() - 1;
			int nxy = nx * (mesh[1].size() - 1);
			layerMaterials.resize(nxy * (mesh[2].size() - 1));

			for (auto layer = layers.begin(); layer != layers.end(); ++layer)
			{
				topIndex = BasicOperations::BinarySearchInVectorSorted(mesh[2], layer->top);
				bottomIndex = BasicOperations::BinarySearchInVectorSorted(mesh[2], layer->bottom);

				for (int iz = bottomIndex; iz < topIndex; iz++)
					for (int iy = 0; iy < mesh[1].size() - 1; iy++)
					{
						cell = iz * nxy + iy * nx;
						for (int ix = 0; ix < nx; ix++, cell++)
							layerMaterials[cell] = layer->materialValue;
					}
			}
			return 0;
		}
		catch (exception ex)
		{
			write_to_log("Error : SetLayersMaterials : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Apply layers materials to mesh
	int SetLayersMaterials2D(vector<Model::GeoLayer<double>> &layers, vector<Model::Material> &layerMaterials, vector<double> *mesh)
	{
		try
		{
			int topIndex, bottomIndex, cell;
			int nx = mesh[0].size() - 1;
			int nxy = nx * (mesh[1].size() - 1);
			layerMaterials.resize(nxy);

			for (auto layer = layers.begin(); layer != layers.end(); ++layer)
			{
				topIndex = BasicOperations::BinarySearchInVectorSorted(mesh[1], layer->top);
				bottomIndex = BasicOperations::BinarySearchInVectorSorted(mesh[1], layer->bottom);

				for (int iz = bottomIndex; iz < topIndex; iz++)
				{
					cell = iz * nx;
					for (int ix = 0; ix < nx; ix++, cell++)
						layerMaterials[cell] = layer->materialValue;
				}
			}
			return 0;
		}
		catch (exception ex)
		{
			write_to_log("Error : SetLayersMaterials : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Apply objects materials to mesh
	int SetObjectsMaterials(vector<Model::GeoBox<double>> &objects, vector<Model::Material> &objectMaterials, vector<double> *mesh)
	{
		try
		{
			int ind[6], cell;
			int nx = mesh[0].size() - 1;
			int nxy = nx * (mesh[1].size() - 1);

			for (auto obj = objects.begin(); obj != objects.end(); ++obj)
			{
				for (int i = 0; i < 6; i++)
					ind[i] = BasicOperations::BinarySearchInVectorSorted(mesh[i / 2], obj->coordinates[i]);

				for (int iz = ind[4]; iz < ind[5]; iz++)
					for (int iy = ind[2]; iy < ind[3]; iy++)
						for (int ix = ind[0]; ix < ind[1]; ix++)
							objectMaterials[iz * nxy + iy * nx + ix] = obj->materialValue;
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

	// Set primary field calculation area metarial parameters
	int SetNormalMaterials(vector<Model::Material> &layerMaterials, vector<Model::Material> &objectMaterials)
	{
		try
		{
			for (int materialIndex = 0; materialIndex < objectMaterials.size(); materialIndex++)
				if (!objectMaterials[materialIndex].SigmaH.Equal(layerMaterials[materialIndex].SigmaH))
					objectMaterials[materialIndex].SigmaN.constantValue = layerMaterials[materialIndex].SigmaH.constantValue;
			return 0;
		}
		catch (exception ex)
		{
			write_to_log("Error : SetNormalMaterials : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Set material numbers to cells
	int SetMaterialNumbers(vector<Model::Material> &objectMaterials, vector<Model::Material> &layerMaterials, vector<int> &materialNumbers, vector<Model::Material> &uniqueMaterials, vector<double> &meshZ, vector<Model::GeoLayer<double>> &layers, int nxy)
	{
		vector<bool> done;
		done.assign(objectMaterials.size(), false);
		materialNumbers.resize(objectMaterials.size());
		try
		{
			int materialNumber = 1;

			// air
			int index1 = objectMaterials.size() - 1;
			materialNumbers[index1] = materialNumber;

			auto *material1 = &objectMaterials[index1];
			material1->numberInNormalTask = 1;
			uniqueMaterials.push_back(*material1);

			for (int index2 = index1 - 1; index2 >= 0; index2--)
			{
				if (done[index2])
					continue;

				if (!material1->Equal(objectMaterials[index2]))
					continue;

				materialNumbers[index2] = materialNumber;
				done[index2] = true;
			}
			done[index1] = true;
			materialNumber++;

			// layers
			int zIndex = 0;
			int layerIndex = layers.size() - 1;
			material1 = &layers[layerIndex].materialValue;
			material1->numberInNormalTask = layerIndex + 1;
			uniqueMaterials.push_back(*material1);
			
			for (index1 = 0; index1 < layerMaterials.size(); index1 += nxy, zIndex++)
			{
				if (done[index1])
					continue;

				
				double z = (meshZ[zIndex] + (zIndex == meshZ.size() - 1 ? meshZ[zIndex] : meshZ[zIndex + 1])) * 0.5;
				while (z > layers[layerIndex].top)
				{
					materialNumber++;
					layerIndex--;
					material1 = &layers[layerIndex].materialValue;
					material1->numberInNormalTask = layerIndex + 1;
					uniqueMaterials.push_back(*material1);
				}

				for (int index2 = index1 + 1; index2 < index1 + nxy; index2++)
				{
					if (done[index2]) continue;
					if (!material1->Equal(objectMaterials[index2])) continue;

					materialNumbers[index2] = materialNumber;
					done[index2] = true;
				}

				if (material1->Equal(objectMaterials[index1]))
				{
					materialNumbers[index1] = materialNumber;
					done[index1] = true;
				}
			}

			materialNumber++;
			// other

			layerIndex = layers.size() - 1;
			for (index1 = 0; index1 < objectMaterials.size(); index1++)
			{
				if (done[index1])
					continue;


				zIndex = index1 / nxy;
				double z = (meshZ[zIndex] + (zIndex == meshZ.size() - 1 ? meshZ[zIndex] : meshZ[zIndex + 1])) * 0.5;
				while (z > layers[layerIndex].top)
					layerIndex--;

				materialNumbers[index1] = materialNumber;

				auto &material1 = objectMaterials[index1];
				material1.numberInNormalTask = layerIndex + 1;
				uniqueMaterials.push_back(material1);

				for (int index2 = index1 + 1; index2 < objectMaterials.size(); index2++)
				{
					if (done[index2])
						continue;

					if (!material1.Equal(objectMaterials[index2]))
						continue;

					materialNumbers[index2] = materialNumber;
					done[index2] = true;
				}
				done[index1] = true;
				materialNumber++;
			}
			return 0;
		}
		catch (exception ex)
		{
			write_to_log("Error : SetMaterialNumbers : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Set material numbers to cells
	int SetMaterialNumbers2D(vector<Model::Material> &layerMaterials, vector<int> &materialNumbers, vector<Model::Material> &uniqueMaterials, vector<double> &meshZ, vector<Model::GeoLayer<double>> &layers, int nx, bool air)
	{
		materialNumbers.resize(layerMaterials.size());
		try
		{
			int materialNumber = 1;
			
			int indexZTop = meshZ.size() - 2;
			for (auto layer = layers.begin(); layer != layers.end(); ++layer)
			{
				uniqueMaterials.push_back(layer->materialValue);
				int indexZBottom = indexZTop;
				double z = (meshZ[indexZBottom] + meshZ[indexZBottom + 1]) * 0.5;
				while (z > layer->bottom && indexZBottom >= 0)
				{
					--indexZBottom;
					if (indexZBottom >= 0)
						z = (meshZ[indexZBottom] + meshZ[indexZBottom + 1]) * 0.5;
				}

				for (int index2 = (indexZBottom + 1) * nx; index2 < (indexZTop + 1) * nx; ++index2)
					materialNumbers[index2] = materialNumber;
				++materialNumber;

				indexZTop = indexZBottom;
			}

			return 0;
		}
		catch (exception ex)
		{
			write_to_log("Error : SetMaterialNumbers : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Set materials to mesh
	int SetMaterials3D(vector<Model::GeoLayer<double>> &layers, vector<Model::GeoBox<double>> &objects, vector<int> &materialNumbers, vector<double> *mesh, vector<Model::Material> &uniqueMaterials)
	{
		try
		{
			vector<Model::Material> layerMaterials;
			vector<Model::Material> objectMaterials;
			if (SetLayersMaterials3D(layers, layerMaterials, mesh) != 0) { write_to_log("Error : SetMaterials : Could not set layers materials to mesh\n"); return 1; }
			objectMaterials = layerMaterials;
			if (SetObjectsMaterials(objects, objectMaterials, mesh) != 0) { write_to_log("Error : SetMaterials : Could not set objects materials to mesh\n"); return 1; }
			if (SetNormalMaterials(layerMaterials, objectMaterials) != 0) { write_to_log("Error : SetMaterials : Could not set normal parameters to materials\n"); return 1; }
			if (SetMaterialNumbers(objectMaterials, layerMaterials, materialNumbers, uniqueMaterials, mesh[2], layers, (mesh[0].size() - 1) * (mesh[1].size() - 1)) != 0) { write_to_log("Error : SetMaterials : Could not set material numbers\n"); return 1; }

			return 0;
		}
		catch (exception ex)
		{
			write_to_log("Error : SetMaterials : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Set materials to mesh
	int SetMaterials2D(vector<Model::GeoLayer<double>> &layers, vector<int> &materialNumbers, vector<double> *mesh, vector<Model::Material> &uniqueMaterials, bool air)
	{
		try
		{
			vector<Model::Material> layerMaterials;
			if (SetLayersMaterials2D(layers, layerMaterials, mesh) != 0) { write_to_log("Error : SetMaterials2D : Could not set layers materials to mesh\n"); return 1; }
			if (SetMaterialNumbers2D(layerMaterials, materialNumbers, uniqueMaterials, mesh[1], layers, (mesh[0].size() - 1), air) != 0) { write_to_log("Error : SetMaterials2D : Could not set material numbers\n"); return 1; }

			return 0;
		}
		catch (exception ex)
		{
			write_to_log("Error : SetMaterials : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	void RefineMesh(vector<double> *meshes1D, vector<Model::Receiver> &receivers, Mesh3D::Settings &meshSettings)
	{
		vector<double> coordinates(receivers.size());

		for (size_t i = 0; i < 2; i++)
		{
			transform(receivers.begin(), receivers.end(), coordinates.begin(), [i](const Model::Receiver& rec) {return rec.coordinates[i]; });
			Mesh1D::RefineMesh(meshes1D[i], coordinates, meshSettings.receiversRefinementsCount, meshSettings.eps);
		}
	}

	// Mesh building procedure
	int BuildMesh3D(vector<Model::Receiver> &receivers, vector<Model::Generator> &generators, vector<Model::GeoLayer<double>> &layers, vector<Model::GeoBox<double>> &objects, Mesh3D::Settings &meshSettings, vector<double> *meshes1D, vector<double> *meshesTemplate1D, vector<int> &materialNumbers, set<pair<int, int>> &edges, vector<vector<int>> &elementsByEdges, vector<Model::Material> &uniqueMaterials)
	{
		try
		{
			vector<double> zBounds(layers.size() + 1);
			for (size_t layerIndex = 0; layerIndex < layers.size(); layerIndex++)
				zBounds[layers.size() - layerIndex] = layers[layerIndex].top;

			double *bounds[3] = { new double[4], new double[4], new double[4] };
			if (CalculateModelBounds3D(receivers, generators, layers, meshSettings, bounds) != 0) { write_to_log("Error : BuildMesh3D : Could not calculate model bounds\n"); return 1; }
			zBounds.front() = bounds[2][0];
			zBounds.back() = bounds[2][3];
			if (Mesh3D::BuildBaseMesh(meshSettings, bounds, meshes1D) != 0) { write_to_log("Error : BuildMesh3D : Could not build base mesh\n"); return 1; }
			if (Mesh3D::AlignMeshBounds(bounds, meshes1D) != 0) { write_to_log("Error : BuildMesh3D : Could not align base mesh\n"); return 1; }
			if (InsertObjectsIntoMesh(objects, meshSettings.eps, meshes1D) != 0) { write_to_log("Error : BuildMesh3D : Could not insert objects into mesh\n"); return 1; }
			if (InsertLayersIntoMesh(layers, meshSettings.eps, meshes1D[2]) != 0) { write_to_log("Error : BuildMesh3D : Could not insert layers into mesh\n"); return 1; }
			RefineMesh(meshes1D, receivers, meshSettings);
			if (Mesh3D::BuildBaseMesh3DTemplate(bounds, zBounds, meshes1D, meshesTemplate1D) != 0) { write_to_log("Error : BuildMesh3D : Could not build template mesh\n"); return 1; }
			if (SetMaterials3D(layers, objects, materialNumbers, meshes1D, uniqueMaterials) != 0) { write_to_log("Error : BuildMesh3D : Could not set material numbers to mesh\n"); return 1; }			
			if (Mesh3D::BuildEdges(meshes1D, edges) != 0) { write_to_log("Error : BuildMesh3D : Could not build edges\n"); return 1; }
			if (Mesh3D::BuildElementsByEdges(meshes1D, edges, elementsByEdges) != 0) { write_to_log("Error : BuildMesh3D : Could not build elements by edges\n"); return 1; }

			

			delete[] bounds[0];
			delete[] bounds[1];
			delete[] bounds[2];
			return 0;
		}
		catch (exception ex)
		{
			write_to_log("Error : BuildMesh : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Mesh building procedure
	int BuildMesh2D(vector<Model::Receiver> &receivers, vector<Model::Generator> &generators, vector<Model::GeoLayer<double>> &layers, Mesh2D::Settings &meshSettings, vector<double> *meshes1D, vector<int> &materialNumbers, vector<Model::Material> &uniqueMaterials, bool air)
	{
		try
		{
			double *bounds[2] = { new double[4], new double[4] };
			if (CalculateModelBounds2D(receivers, generators, meshSettings, bounds) != 0) { write_to_log("Error : BuildMesh2D : Could not calculate model bounds\n"); return 1; }

			bounds[0][0] = 0.5;
			bounds[0][1] = 0.5;
			bounds[0][2] = 0.5;
			bounds[0][3] = meshSettings.farBound[0];

			if (!air)
			{
				bounds[1][3] = 0.0;
				if (bounds[1][2] >= 0.0)
					bounds[1][2] = 0.0;
			}

			if (Mesh2D::BuildBaseMesh(meshSettings, bounds, meshes1D) != 0) { write_to_log("Error : BuildMesh2D : Could not build base mesh\n"); return 1; }
			if (Mesh2D::AlignMeshBounds(bounds, meshes1D) != 0) { write_to_log("Error : BuildMesh2D : Could not align base mesh\n"); return 1; }
			if (InsertLayersIntoMesh(layers, meshSettings.eps, meshes1D[1]) != 0) { write_to_log("Error : BuildMesh2D : Could not insert layers into mesh\n"); return 1; }
			if (SetMaterials2D(layers, materialNumbers, meshes1D, uniqueMaterials, air) != 0) { write_to_log("Error : BuildMesh2D : Could not set material numbers to mesh\n"); return 1; }
			
			delete[] bounds[0];
			delete[] bounds[1];
			return 0;
		}
		catch (exception ex)
		{
			write_to_log("Error : BuildMesh : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	int PrepareMesh3D(Mesh3D::Settings &meshSettings , const std::chrono::high_resolution_clock::time_point &tProgramStart, vector<Model::Receiver> &receivers, vector<Model::Generator> &generators
		, vector<Model::GeoLayer<double>> &layers, vector<Model::GeoBox<double>> &objects)
	{
		
		vector<Model::Material> materials;
		vector<double> meshes1D[3];
		vector<double> meshesTemplate1D[3];
		vector<int> materialNumbers;
		set<pair<int, int>> edges;
		vector<vector<int>> elementsByEdges;

		fprintf(log_file, "%lf :", std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> >>(std::chrono::high_resolution_clock::now() - tProgramStart).count());
		write_to_log("Building mesh\n");
		if (BuildMesh3D(receivers, generators, layers, objects, meshSettings, meshes1D, meshesTemplate1D, materialNumbers, edges, elementsByEdges, materials) != 0)
		{
			write_to_log("Error : MainProcedure : Could not build mesh\n");
			return 1;
		}

		fprintf(log_file, "%lf :", std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> >>(std::chrono::high_resolution_clock::now() - tProgramStart).count());
		write_to_log("Writing materials\n");
		if (ReadWrite::WriteMaterials3D(materials) != 0)
		{
			write_to_log("Error : MainProcedure : Could not write materials\n");
			return 1;
		}

		fprintf(log_file, "%lf :", std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> >>(std::chrono::high_resolution_clock::now() - tProgramStart).count());
		write_to_log("Writing mesh\n");
		if (ReadWrite::WriteMesh3D(meshes1D, meshesTemplate1D, materialNumbers, edges, elementsByEdges) != 0)
		{
			write_to_log("Error : MainProcedure : Could not write mesh\n");
			return 1;
		}

		return 0;
	}

	int PrepareMesh2D(const std::chrono::high_resolution_clock::time_point &tProgramStart, Mesh2D::Settings &meshSettings, vector<Model::Receiver> &receivers, vector<Model::Generator> &generators, vector<Model::GeoLayer<double>> &layers, int index)
	{
		vector<Model::Material> materials;
		vector<double> meshes1D[2];
		vector<double> meshesTemplate1D[2];
		vector<int> materialNumbers;

		fprintf(log_file, "%lf :", std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> >>(std::chrono::high_resolution_clock::now() - tProgramStart).count());
		write_to_log("Building mesh\n");
		if (BuildMesh2D(receivers, generators, layers, meshSettings, meshes1D, materialNumbers, materials, index == 2 || index == 3) != 0)
		{
			write_to_log("Error : MainProcedure : Could not build mesh\n");
			return 1;
		}

		fprintf(log_file, "%lf :", std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> >>(std::chrono::high_resolution_clock::now() - tProgramStart).count());
		write_to_log("Writing materials\n");
		if (ReadWrite::WriteMaterials2D(materials) != 0)
		{
			write_to_log("Error : MainProcedure : Could not write materials\n");
			return 1;
		}

		fprintf(log_file, "%lf :", std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> >>(std::chrono::high_resolution_clock::now() - tProgramStart).count());
		write_to_log("Writing mesh\n");
		if (ReadWrite::WriteMesh2D(meshes1D, materialNumbers, index) != 0)
		{
			write_to_log("Error : MainProcedure : Could not write mesh\n");
			return 1;
		}

		return 0;
	}

	int PrepareMeshes2D(Mesh2D::Settings &mesh2DSettings, const std::chrono::high_resolution_clock::time_point &tProgramStart, vector<Model::Receiver> &receivers, vector<Model::Generator> &generators, vector<Model::GeoLayer<double>> &layers)
	{
		/*
		Mesh2D::Settings meshSettings;
		meshSettings.farBound[0] = meshSettings.farBound[1] = 2e+5;
		meshSettings.gapFromReceivers[0] = meshSettings.gapFromReceivers[1] = 0;
		meshSettings.sparse[0] = 1.05;
		meshSettings.sparse[1] = 1.05;
		meshSettings.steps[0] = meshSettings.steps[1] = 0.2;
		*/

		layers.front().materialValue.SigmaH.constantValue
			= layers.front().materialValue.SigmaV.constantValue
			= layers.front().materialValue.SigmaN.constantValue
			= 0.0;


		char buf[2048];
		
		for (int i = 1; i < 5; i++)
		{
			fprintf(log_file, "%lf :", std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> >>(std::chrono::high_resolution_clock::now() - tProgramStart).count());
			sprintf(buf, "Building 2D mesh %d\n", i);
			write_to_log(buf);

			if (PrepareMesh2D(tProgramStart, mesh2DSettings, receivers, generators, layers, i) != 0)
				return 1;
		}

		return 0;
	}

	int PrepareMeshTime(const std::chrono::high_resolution_clock::time_point &tProgramStart, MeshTime::Settings &settings)
	{
		fprintf(log_file, "%lf :", std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> >>(std::chrono::high_resolution_clock::now() - tProgramStart).count());
		write_to_log("Building time mesh\n");

		vector<double> mesh;
		double bounds[2] {0, settings.timeLast };
		Mesh1D::Build1DSparsedMesh(settings.timeStep, settings.timeSparse, 1, bounds, mesh);

		if (!std::filesystem::exists("currentfunction") && ReadWrite::WriteCurrentFunction("currentfunction", mesh) != 0)
		{
			write_to_log("Error : MainProcedure : Could not write currentfunction\n");
			return 1;
		}

		if (!std::filesystem::exists("deltafunction") && ReadWrite::WriteDeltaFunction("deltafunction", mesh) != 0)
		{
			write_to_log("Error : MainProcedure : Could not write deltafunction\n");
			return 1;
		}

		if (!std::filesystem::exists("infite.0") && ReadWrite::WriteInfite0("infite.0", mesh) != 0)
		{
			write_to_log("Error : MainProcedure : Could not write infite.0\n");
			return 1;
		}

		return 0;
	}
	
	// Main mesh building and files writing procedure
	int MainProcedure()
	{
		char buf[2048];
		vector<Model::Receiver> receivers;
		vector<Model::Generator> generators;
		vector<Model::GeoLayer<double>> layers;
		vector<Model::GeoBox<double>> objects;
		

		// Program start time
		std::chrono::high_resolution_clock::time_point tProgramStart = std::chrono::high_resolution_clock::now();;


		fprintf(log_file, "%lf :", std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> >>(std::chrono::high_resolution_clock::now() - tProgramStart).count());
		write_to_log("Reading layers\n");

		if (ReadWrite::ReadZSig2D("z_sig_2d", layers) != 0)
		{
			if (ReadWrite::ReadZSig2D("../z_sig_2d", layers) != 0)
			{
				write_to_log("Error : MainProcedure : Could not read z_sig_2d\n");
				return 1;
			}
		}

		fprintf(log_file, "%lf :", std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> >>(std::chrono::high_resolution_clock::now() - tProgramStart).count());
		write_to_log("Reading objects\n");
		if (ReadWrite::ReadObjects("objects", objects) != 0)
		{
			if (ReadWrite::ReadObjects("../objects", objects) != 0)
			{
				write_to_log("Error : MainProcedure : Could not read objects\n");
				return 1;
			}
		}

		fprintf(log_file, "%lf :", std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> >>(std::chrono::high_resolution_clock::now() - tProgramStart).count());
		write_to_log("Reading receivers\n");
		if (ReadWrite::ReadReceivers("xyzmn", receivers) != 0)
		{
			if (ReadWrite::ReadReceivers("../xyzmn", receivers) != 0)
			{
				write_to_log("Error : MainProcedure : Could not read xyzmn\n");
				return 1;
			}
		}

		if (ReadWrite::ReadGenerators("sours", generators) != 0)
		{
			if (ReadWrite::ReadGenerators("../sours", generators) != 0)
			{
				write_to_log("Error : MainProcedure : Could not read sours\n");
				return 1;
			}
		}

		int refineStepsCount;
		if (ReadWrite::ReadInt("numberofmnpoints", refineStepsCount) != 0)
		{
			write_to_log("Error : MainProcedure : Could not read numberofmnpoints\n");
			return 1;
		}
		
		if (ReadWrite::WritexyzVectorE("xyzVectorE", receivers, refineStepsCount) != 0)
		{
			write_to_log("Error : MainProcedure : Could not write xyzVectorE\n");
			return 1;
		}

		if (ReadWrite::WriteRecvsE("recvse", receivers, refineStepsCount, generators.size() / 2) != 0)
		{
			write_to_log("Error : MainProcedure : Could not write recvse\n");
			return 1;
		}

		if (ReadWrite::WriteGroup("group", generators) != 0)
		{
			write_to_log("Error : MainProcedure : Could not write group\n");
			return 1;
		}

		if (ReadWrite::WriteGens("gens", generators) != 0)
		{
			write_to_log("Error : MainProcedure : Could not write gens\n");
			return 1;
		}
		

		Mesh3D::Settings mesh3DSettings;
		Mesh2D::Settings mesh2DSettings;
		MeshTime::Settings timeSettings;
		
		fprintf(log_file, "%lf :", std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> >>(std::chrono::high_resolution_clock::now() - tProgramStart).count());
		write_to_log("Reading settings\n");
		if (ReadWrite::ReadSettings("RegularMeshBuilderSettings.cfg", mesh3DSettings, mesh2DSettings, timeSettings) != 0)
		{
			if (ReadWrite::ReadSettings("../RegularMeshBuilderSettings.cfg", mesh3DSettings, mesh2DSettings, timeSettings) != 0)
			{
				write_to_log("Error : MainProcedure : Could not read RegularMeshBuilderSettings\n");
				return 1;
			}
		}
		/*
		fprintf(log_file, "%lf :", std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> >>(std::chrono::high_resolution_clock::now() - tProgramStart).count());
		write_to_log("Preparing 3D mesh\n");
		if (PrepareMesh3D(mesh3DSettings, tProgramStart, receivers, generators, layers, objects) != 0)
		{
			write_to_log("Error : MainProcedure : Could not prepare 3D mesh\n");
			return 1;
		}*/
		
		fprintf(log_file, "%lf :", std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> >>(std::chrono::high_resolution_clock::now() - tProgramStart).count());
		write_to_log("Preparing 2D meshes\n");
		if (PrepareMeshes2D(mesh2DSettings, tProgramStart, receivers, generators, layers) != 0)
		{
			write_to_log("Error : MainProcedure : Could not prepare 2D meshes\n");
			return 1;
		}

		if (PrepareMeshTime(tProgramStart, timeSettings) != 0)
		{
			write_to_log("Error : MainProcedure : Could not prepare time meshe\n");
			return 1;
		}

		return 0;
	}
}