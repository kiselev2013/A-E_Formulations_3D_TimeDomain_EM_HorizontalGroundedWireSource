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

#include "ReadWrite.h"
#include "Logging.h"


namespace ReadWrite
{
	int ReadInt(char *file_name, int &val)
	{
		FILE *file_in = NULL;
		try
		{
			if (!(file_in = fopen(file_name, "r"))) return 1;

			fscanf(file_in, "%d", &val);

			fclose(file_in);
			return 0;
		}
		catch (exception ex)
		{
			if (file_in != NULL)
				fclose(file_in);
			write_to_log("Error : ReadInt : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Read objects from file
	int ReadObjects(char *file_name, vector<Model::GeoBox<double>> &objects)
	{
		FILE *file_in = NULL;
		try
		{
			if (!(file_in = fopen(file_name, "r"))) return 1;

			Model::Material newMaterial;
			Model::GeoBox<double> newObject;

			int objectsCount;
			fscanf(file_in, "%d", &objectsCount);
			objects.resize(objectsCount);

			for (int objectIndex = 0; objectIndex < objectsCount; objectIndex++)
			{
				fscanf(file_in, "%lf%lf%lf%lf%lf%lf%lf"
					, newObject.coordinates
					, newObject.coordinates + 1
					, newObject.coordinates + 2
					, newObject.coordinates + 3
					, newObject.coordinates + 4
					, newObject.coordinates + 5
					, &newMaterial.SigmaH.constantValue);
				newMaterial.SigmaV.constantValue = newMaterial.SigmaH.constantValue;
				newMaterial.SigmaN.constantValue = newMaterial.SigmaH.constantValue;
				newObject.materialValue = newMaterial;
				objects[objectIndex] = newObject;
			}

			fclose(file_in);
			return 0;
		}
		catch (exception ex)
		{
			if (file_in != NULL)
				fclose(file_in);
			write_to_log("Error : ReadObjects : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Read layers from file
	int ReadZSig2D(char *file_name, vector<Model::GeoLayer<double>> &layers)
	{
		FILE *file_in = NULL;
		try
		{
			Model::Material newMaterial;
			Model::GeoLayer<double> newLayer;
			double z, sigma;
			int layersCount;
			if (!(file_in = fopen(file_name, "r"))) return 1;

			fscanf(file_in, "%d%", &layersCount);

			layers.resize(layersCount + 1);

			CreateAirMaterial(newMaterial);
			layers[0].materialValue = newMaterial;
			layers[0].bottom = 0;
			layers[0].top = 1e+30;
			layers[0].material = 0;
			for (int layerIndex = 0; layerIndex < layersCount; layerIndex++)
			{
				fscanf(file_in, "%lf%lf", &z, &sigma);
				newLayer.top = z;
				newLayer.bottom = layerIndex == 0 ? -1e+30 : layers[layers.size() - layerIndex].top;
				newMaterial.SigmaH.constantValue = sigma;
				newMaterial.SigmaV.constantValue = sigma;
				newMaterial.SigmaN.constantValue = sigma;
				newLayer.materialValue = newMaterial;
				layers[layers.size() - layerIndex - 1] = newLayer;
			}

			fclose(file_in);
			return 0;
		}
		catch (exception ex)
		{
			if (file_in != NULL)
				fclose(file_in);
			write_to_log("Error : ReadZSig2D : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Read receivers from file
	int ReadReceivers(char *file_name, vector<Model::Receiver> &receivers)
	{
		FILE *file_in = NULL;
		try
		{
			if (!(file_in = fopen(file_name, "r"))) return 1;

			Model::Receiver newReceiver;

			int receiversCount;
			fscanf(file_in, "%d", &receiversCount);
			receiversCount *= 2;
			receivers.resize(receiversCount);

			for (int receiverIndex = 0; receiverIndex < receiversCount; receiverIndex++)
			{
				fscanf(file_in, "%lf%lf%lf", newReceiver.coordinates, newReceiver.coordinates + 1, newReceiver.coordinates + 2);
				receivers[receiverIndex] = newReceiver;
			}

			fclose(file_in);
			return 0;
		}
		catch (exception ex)
		{
			if (file_in != NULL)
				fclose(file_in);
			write_to_log("Error : ReadReceivers : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Read generators from file
	int ReadGenerators(char *file_name, vector<Model::Generator> &generators)
	{
		FILE *file_in = NULL;
		try
		{
			if (!(file_in = fopen(file_name, "r"))) return 1;

			Model::Generator newGenerator;
			while (fscanf(file_in, "%lf%lf%lf", newGenerator.coordinates, newGenerator.coordinates + 1, newGenerator.coordinates + 2) == 3)
			{
				generators.push_back(newGenerator);
			}

			fclose(file_in);
			return 0;
		}
		catch (exception ex)
		{
			if (file_in != NULL)
				fclose(file_in);
			write_to_log("Error : ReadReceivers : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Read settings from file
	int ReadSettings(char *file_name, Mesh3D::Settings &mesh3DSettings, Mesh2D::Settings &mesh2DSettings, MeshTime::Settings &timeSettings)
	{
		FILE *file_in = NULL;
		try
		{
			if (!(file_in = fopen(file_name, "r"))) return 1;

			fscanf(file_in, "%lf%lf", mesh3DSettings.steps, mesh3DSettings.steps + 2); mesh3DSettings.steps[1] = mesh3DSettings.steps[0];
			fscanf(file_in, "%lf%lf", mesh3DSettings.sparse, mesh3DSettings.sparse + 2); mesh3DSettings.sparse[1] = mesh3DSettings.sparse[0];
			fscanf(file_in, "%lf%lf", mesh3DSettings.gapFromReceivers, mesh3DSettings.gapFromReceivers + 2); mesh3DSettings.gapFromReceivers[1] = mesh3DSettings.gapFromReceivers[0];
			fscanf(file_in, "%lf%lf", mesh3DSettings.farBound, mesh3DSettings.farBound + 2); mesh3DSettings.farBound[1] = mesh3DSettings.farBound[0];
			fscanf(file_in, "%d", &mesh3DSettings.receiversRefinementsCount);

			fscanf(file_in, "%lf%lf", mesh2DSettings.steps, mesh2DSettings.steps + 1);
			fscanf(file_in, "%lf%lf", mesh2DSettings.sparse, mesh2DSettings.sparse + 1);
			fscanf(file_in, "%lf%lf", mesh2DSettings.gapFromReceivers, mesh2DSettings.gapFromReceivers + 1);
			fscanf(file_in, "%lf%lf", mesh2DSettings.farBound, mesh2DSettings.farBound + 1);

			
			fscanf(file_in, "%lf%lf%lf", &timeSettings.timeStep, &timeSettings.timeSparse, &timeSettings.timeLast);

			
			
			fclose(file_in);
			return 0;
		}
		catch (exception ex)
		{
			if (file_in != NULL)
				fclose(file_in);
			write_to_log("Error : ReadSettings : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Write mesh file with nodes and cells count
	int WriteInftry(char *file_name, int nx, int ny, int nz)
	{
		FILE *file_out = NULL;
		try
		{
			if (!(file_out = fopen(file_name, "w"))) return 1;

			int nodesCount = nx * ny * nz;
			int cellsCount = (nx - 1) * (ny - 1) * (nz - 1);
			int bc1NodesCount = (nx * ny + nx * (nz - 2) + (ny - 2) * (nz - 2)) * 2;

			fprintf(file_out, " ISLAU= %d", 0);
			fprintf(file_out, " INDKU1= %d", 0);
			fprintf(file_out, " INDFPO= %d", 0);
			fprintf(file_out, "\nKUZLOV= %d", nodesCount);
			fprintf(file_out, "   KPAR= %d", cellsCount);
			fprintf(file_out, "    KT1= %d", bc1NodesCount);
			fprintf(file_out, "   KTR2= %d", 0);
			fprintf(file_out, "   KTR3= %d\n", 0);
			fprintf(file_out, "   KT8= %d", 0);
			fprintf(file_out, "   KT9= %d\n", 0);
			fprintf(file_out, "KISRS1= %d", 0);
			fprintf(file_out, " KISRS2= %d", 0);
			fprintf(file_out, " KISRS3= %d", 0);
			fprintf(file_out, "   KBRS= %d\n", 0);
			fprintf(file_out, "   KT7= %d", 0);
			fprintf(file_out, "   KT10= %d", 0);
			fprintf(file_out, "  KTR4= %d", 0);
			fprintf(file_out, "  KTSIM= %d\n", 0);
			fprintf(file_out, "   KT6= %d", 0);
			fclose(file_out);

			return 0;
		}
		catch (exception ex)
		{
			if (file_out != NULL)
				fclose(file_out);
			write_to_log("Error : WriteInftry : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Write mesh nodes coordinates
	int WriteXyz(char *file_name, vector<double> *mesh)
	{
		FILE *file_out = NULL;
		try
		{
			double point[3];

			if (!(file_out = fopen(file_name, "wb"))) return 1;

			auto &meshX = mesh[0];
			auto &meshY = mesh[1];
			auto &meshZ = mesh[2];

			for (int iz = 0; iz < meshZ.size(); iz++)
			{
				point[2] = meshZ[iz];
				for (int iy = 0; iy < meshY.size(); iy++)
				{
					point[1] = meshY[iy];
					for (int ix = 0; ix < meshX.size(); ix++)
					{
						point[0] = meshX[ix];
						fwrite(point, sizeof(double), 3, file_out);
					}
				}
			}

			fclose(file_out);
			return 0;
		}
		catch (exception ex)
		{
			if (file_out != NULL)
				fclose(file_out);
			write_to_log("Error : WriteXyz : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Write mesh cells described with node numbers
	int WriteNver(char *file_name, int nx, int ny, int nz)
	{
		FILE *file_out = NULL;
		try
		{
			int element[14];
			element[8] = 0;
			element[9] = 0;
			element[10] = 0;
			element[11] = 0;
			element[12] = 0;
			element[13] = 0;
			int nxy = nx * ny;

			if (!(file_out = fopen(file_name, "wb"))) return 1;

			int node;
			for (int iz = 0; iz < nz - 1; iz++)
			{
				for (int iy = 0; iy < ny - 1; iy++)
				{
					node = nxy * iz + nx * iy + 1;
					element[0] = node;
					element[1] = node + 1;
					element[2] = node + nx;
					element[3] = node + nx + 1;
					element[4] = node + nxy;
					element[5] = node + nxy + 1;
					element[6] = node + nxy + nx;
					element[7] = node + nxy + nx + 1;
					for (int ix = 0; ix < nx - 1; ix++)
					{
						fwrite(element, sizeof(int), 14, file_out);

						element[0]++;
						element[1]++;
						element[2]++;
						element[3]++;
						element[4]++;
						element[5]++;
						element[6]++;
						element[7]++;
					}
				}
			}
			fclose(file_out);

			return 0;
		}
		catch (exception ex)
		{
			if (file_out != NULL)
				fclose(file_out);
			write_to_log("Error : WriteNver : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Write int array to binary file
	int WriteIntArrayBinary(char *file_name, vector<int> &arr)
	{
		FILE *file_out = NULL;
		try
		{
			int m = 1;

			if (!(file_out = fopen(file_name, "wb"))) return 1;

			for (int i = 0; i < arr.size(); i++)
			{
				m = arr[i];
				fwrite(&m, sizeof(int), 1, file_out);
			}
			fclose(file_out);

			return 0;
		}
		catch (exception ex)
		{
			if (file_out != NULL)
				fclose(file_out);
			write_to_log("Error : WriteNvkat : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Write double array to binary file
	int WriteDoubleArrayBinary(char *file_name, vector<double> &arr)
	{
		FILE *file_out = NULL;
		try
		{
			if (!(file_out = fopen(file_name, "wb"))) return 1;
			fwrite(arr.data(), sizeof(double), arr.size(), file_out);

			fclose(file_out);
			return 0;
		}
		catch (exception ex)
		{
			if (file_out != NULL)
				fclose(file_out);
			write_to_log("Error : WriteDoubleArrayBinary : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Write each cell neighbours to file
	int WriteElemNeib(char *file_name, int nx, int ny, int nz)
	{
		FILE *file_out = NULL;
		try
		{
			int m = 0;
			int neighbours[6];

			if (!(file_out = fopen(file_name, "wb"))) return 1;

			int cell = 0;
			for (int iz = 0; iz < nz - 1; iz++)
			{
				for (int iy = 0; iy < ny - 1; iy++)
				{
					for (int ix = 0; ix < nx - 1; ix++, cell++)
					{
						Mesh3D::GetNeighbours(nx, ny, nz, ix, iy, iz, cell, neighbours);
						for (int i = 0; i < 6; i++)
						{
							if (neighbours[i] == -1)
							{
								m = 0;
								fwrite(&m, sizeof(int), 1, file_out);
							}
							else
							{
								m = 1;
								fwrite(&m, sizeof(int), 1, file_out);
								fwrite(neighbours + i, sizeof(int), 1, file_out);
							}
						}
					}
				}
			}
			fclose(file_out);

			return 0;
		}
		catch (exception ex)
		{
			if (file_out != NULL)
				fclose(file_out);
			write_to_log("Error : WriteElemNeib : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Write material properties to file
	int WriteMtr3D2D(char *fileName, vector<Model::Material> &materials)
	{
		FILE *outputFile = NULL;
		try
		{
			int n;

			if (open_file_w(fileName, &outputFile) != 0)
				return 1;

			for (int i = 0; i < materials.size(); i++)
				fprintf(outputFile, "%d\t%d\n", i + 1, materials[i].numberInNormalTask);

			fclose(outputFile);
			return 0;
		}
		catch (exception ex)
		{
			if (outputFile != NULL)
				fclose(outputFile);
			write_to_log("Error : WriteProperty3D : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Write material properties to file
	int WriteProperty3D(char *fileName, vector<Model::Material> &materials, Model::PhysicalValueName prop1, Model::PhysicalValueName prop2)
	{
		FILE *outputFile = NULL;
		try
		{
			int n;

			if (open_file_w(fileName, &outputFile) != 0)
				return 1;

			for (int i = 0; i < materials.size(); i++)
				fprintf(outputFile, "%d\t%.7e\t%.7e\n", i + 1, materials[i].GetValue(prop1).constantValue, materials[i].GetValue(prop2).constantValue);

			fclose(outputFile);
			return 0;
		}
		catch (exception ex)
		{
			if (outputFile != NULL)
				fclose(outputFile);
			write_to_log("Error : WriteProperty3D : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Write material properties to file
	int WriteProperty2D(char *fileName, vector<Model::Material> &materials, Model::PhysicalValueName prop)
	{
		FILE *outputFile = NULL;
		try
		{
			int n;

			if (open_file_w(fileName, &outputFile) != 0)
				return 1;

			for (int i = 0; i < materials.size(); i++)
				fprintf(outputFile, "%d\t%.7e\n", i + 1, materials[i].GetValue(prop).constantValue);

			fclose(outputFile);
			return 0;
		}
		catch (exception ex)
		{
			if (outputFile != NULL)
				fclose(outputFile);
			write_to_log("Error : WriteProperty3D : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Write empty nodal T-Matrix size
	int WriteTSize3DNode(char *fileName)
	{
		FILE *outputFile = NULL;
		try
		{
			try
			{
				if (open_file_w(fileName, &outputFile) != 0)
					return 1;

				fprintf(outputFile, "%d\n", 0);

				fclose(outputFile);
				return 0;
			}
			catch (exception ex)
			{
				if (outputFile != NULL)
					fclose(outputFile);
				return 1;
			}
			return 0;
		}
		catch (exception ex)
		{
			if (outputFile != NULL)
				fclose(outputFile);
			write_to_log("Error : WriteTSize3DNode : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Write empty edge T-Matrix size
	int WriteTSize3DEdge(char *fileName, int edgesCount)
	{
		FILE *outputFile = NULL;
		try
		{
			if (open_file_w(fileName, &outputFile) != 0)
				return 1;

			fprintf(outputFile, "%d\n%d\n", 0, edgesCount);

			fclose(outputFile);
			return 0;
		}
		catch (exception ex)
		{
			if (outputFile != NULL)
				fclose(outputFile);
			write_to_log("Error : WriteTSize3DNode : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Write 1D meshes to file
	int Write3DMeshRegular(char *fileName, vector<double> *mesh)
	{
		FILE *outputFile = NULL;
		try
		{
			if (open_file_w(fileName, &outputFile) != 0)
				return 1;

			int nx, ny, nz;
			int edgesCount = 0;
			nx = mesh[0].size();
			ny = mesh[1].size();
			nz = mesh[2].size();
			edgesCount += (nx - 1) * ny * nz;
			edgesCount += (ny - 1) * nx * nz;
			edgesCount += (nz - 1) * nx * ny;
			for (int direction = 0; direction < 3; direction++)
			{
				fprintf(outputFile, "%d\n", mesh[direction].size());
				for (int meshIndex = 0; meshIndex < mesh[direction].size(); meshIndex++)
				{
					fprintf(outputFile, "%d\t%.14e\n", meshIndex + 1, mesh[direction][meshIndex]);
				}
			}

			fprintf(outputFile, "0\n0\n");

			fclose(outputFile);
			return 0;
		}
		catch (exception ex)
		{
			if (outputFile != NULL)
				fclose(outputFile);
			write_to_log("Error : Write3DMeshRegular : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Write regular mesh cells flags
	int WriteRegular(char *fileName, int nx, int ny, int nz)
	{
		FILE *outputFile = NULL;
		try
		{
			int m = 0;
			int neighbours[6];

			if (!(outputFile = fopen(fileName, "wb"))) return 1;

			int cellsCount = (nx - 1) * (ny - 1) * (nz - 1);
			fwrite(&cellsCount, sizeof(int), 1, outputFile);
			for (int cell = 1; cell <= cellsCount; cell++)
				fwrite(&cell, sizeof(int), 1, outputFile);
			fclose(outputFile);

			return 0;
		}
		catch (exception ex)
		{
			if (outputFile != NULL)
				fclose(outputFile);
			write_to_log("Error : WriteRegular : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Write boundary nodes
	int WriteL13D(char *fileName, int nx, int ny, int nz)
	{
		FILE *outputFile = NULL;
		try
		{
			int m = 0;
			int neighbours[6];

			if (!(outputFile = fopen(fileName, "wb"))) return 1;

			int node = 1;
			for (int iz = 0; iz < nz; iz++)
				for (int iy = 0; iy < ny; iy++)
					for (int ix = 0; ix < nx; ix++, node++)
						if (ix == 0 || ix == nx - 1 || iy == 0 || iy == ny - 1 || iz == 0 || iz == nz - 1)
							fwrite(&node, sizeof(int), 1, outputFile);
			fclose(outputFile);

			return 0;
		}
		catch (exception ex)
		{
			if (outputFile != NULL)
				fclose(outputFile);
			write_to_log("Error : WriteL13D : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Write node numbers for each edge
	int WriteNodesForEdges(char *fileName, set<pair<int, int>> &edges)
	{
		FILE *outputFile = NULL;
		try
		{
			if (!(outputFile = fopen(fileName, "wb"))) return 1;

			for (auto edge = edges.begin(); edge != edges.end(); ++edge)
			{
				fwrite(&edge->first, sizeof(int), 1, outputFile);
				fwrite(&edge->second, sizeof(int), 1, outputFile);
			}
			fclose(outputFile);

			return 0;
		}
		catch (exception ex)
		{
			if (outputFile != NULL)
				fclose(outputFile);
			write_to_log("Error : WriteNodesForEdges : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Write mesh cells described by edges
	int WriteEdges(char *fileName, vector<vector<int>> &elementsByEdges)
	{
		FILE *outputFile = NULL;
		try
		{
			int zero = 0;
			int one = 1;
			if (!(outputFile = fopen(fileName, "wb"))) return 1;

			for (auto element = elementsByEdges.begin(); element != elementsByEdges.end(); ++element)
			{
				for (int i = 0; i < element->size(); i++)
					fwrite(&((*element)[i]), sizeof(int), 1, outputFile);

				for (int i = 0; i < 12; i++)
					fwrite(&zero, sizeof(int), 1, outputFile);
				fwrite(&one, sizeof(int), 1, outputFile);
			}
			fclose(outputFile);

			return 0;
		}
		catch (exception ex)
		{
			if (outputFile != NULL)
				fclose(outputFile);
			write_to_log("Error : WriteNodesForEdges : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Write empty edge T-Matrix
	int WriteIg3D(char *fileName)
	{
		FILE *outputFile = NULL;
		try
		{
			if (!(outputFile = fopen(fileName, "wb"))) return 1;
			int m = 1;
			fwrite(&m, sizeof(int), 1, outputFile);
			fwrite(&m, sizeof(int), 1, outputFile);
			fclose(outputFile);

			return 0;
		}
		catch (exception ex)
		{
			if (outputFile != NULL)
				fclose(outputFile);
			write_to_log("Error : WriteIg3D : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}
	int WriteJg3D(char *fileName)
	{
		FILE *outputFile = NULL;
		try
		{
			if (!(outputFile = fopen(fileName, "wb"))) return 1;
			fclose(outputFile);

			return 0;
		}
		catch (exception ex)
		{
			if (outputFile != NULL)
				fclose(outputFile);
			write_to_log("Error : WriteJg3D : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}
	int WriteGg3D(char *fileName)
	{
		FILE *outputFile = NULL;
		try
		{
			if (!(outputFile = fopen(fileName, "wb"))) return 1;
			fclose(outputFile);

			return 0;
		}
		catch (exception ex)
		{
			if (outputFile != NULL)
				fclose(outputFile);
			write_to_log("Error : WriteGg3D : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Write receivers coordinates
	int WritePointRes(char *fileName, vector<Model::Receiver> &receivers)
	{
		FILE *outputFile = NULL;
		try
		{
			if (!(outputFile = fopen(fileName, "w"))) return 1;

			fprintf(outputFile, "%d\n", receivers.size());
			for (auto rec = receivers.begin(); rec != receivers.end(); ++rec)
				fprintf(outputFile, "%.13e\t%.13e\t%.13e\n", rec->coordinates[0], rec->coordinates[1], rec->coordinates[2]);

			fclose(outputFile);

			return 0;
		}
		catch (exception ex)
		{
			if (outputFile != NULL)
				fclose(outputFile);
			return 1;
		}
	}

	// Write receivers coordinates
	int WritexyzVectorE(char *fileName, vector<Model::Receiver> &receivers, int numberOfMNDipoles)
	{
		FILE *outputFile = NULL;
		try
		{
			if (!(outputFile = fopen(fileName, "w"))) return 1;

			fprintf(outputFile, "%d\n", receivers.size() * numberOfMNDipoles / 2);
			for (size_t recIndex = 0; recIndex < receivers.size(); recIndex += 2)
			{
				auto &rec1 = receivers[recIndex];
				auto &rec2 = receivers[recIndex + 1];

				for (size_t i = 0; i < numberOfMNDipoles; i++)
				{
					double d = ((double)i) / ((double)(numberOfMNDipoles - 1));

					fprintf(outputFile, "%.13e\t%.13e\t%.13e\n",
						rec1.coordinates[0] * (1.0 - d) + d * rec2.coordinates[0],
						rec1.coordinates[1] * (1.0 - d) + d * rec2.coordinates[1],
						rec1.coordinates[2] * (1.0 - d) + d * rec2.coordinates[2]);
				}
			}

			fclose(outputFile);

			return 0;
		}
		catch (exception ex)
		{
			if (outputFile != NULL)
				fclose(outputFile);
			return 1;
		}
	}

	// Write receivers coordinates count
	int WriteRecvsE(char *fileName, vector<Model::Receiver> &receivers, int numberOfMNDipoles, int generatorsCount)
	{
		FILE *outputFile = NULL;
		try
		{
			if (!(outputFile = fopen(fileName, "w"))) return 1;

			for (size_t i = 0; i < generatorsCount; i++)
				fprintf(outputFile, "%d\n", receivers.size() * numberOfMNDipoles / (2 * generatorsCount));

			fclose(outputFile);

			return 0;
		}
		catch (exception ex)
		{
			if (outputFile != NULL)
				fclose(outputFile);
			return 1;
		}
	}

	// Write group info
	int WriteGroup(char *fileName, vector<Model::Generator> &generators)
	{
		FILE *outputFile = NULL;
		try
		{
			if (!(outputFile = fopen(fileName, "w"))) return 1;

			fprintf(outputFile, "1\n1 1 %d\n", generators.size() / 2);

			fclose(outputFile);

			return 0;
		}
		catch (exception ex)
		{
			if (outputFile != NULL)
				fclose(outputFile);
			return 1;
		}
	}

	// Write gens
	int WriteGens(char *fileName, vector<Model::Generator> &generators)
	{
		FILE *outputFile = NULL;
		try
		{
			if (!(outputFile = fopen(fileName, "w"))) return 1;

			fprintf(outputFile, "%d\n", generators.size() / 2);
			for (size_t i = 0; i < generators.size() / 2; i++)
			{
				fprintf(outputFile, "%d\t%.13e\t%.13e\t%.13e\t%.13e\t%.13e\t%.13e\n"
					, i + 1
					, generators[i * 2].coordinates[0]
					, generators[i * 2].coordinates[1]
					, generators[i * 2].coordinates[2]
					, generators[i * 2 + 1].coordinates[0]
					, generators[i * 2 + 1].coordinates[1]
					, generators[i * 2 + 1].coordinates[2]);
			}

			fclose(outputFile);

			return 0;
		}
		catch (exception ex)
		{
			if (outputFile != NULL)
				fclose(outputFile);
			return 1;
		}
	}

	// Write 3D problem SLAE size
	int WriteKuSlau2(char *fileName, int edgesCount)
	{
		FILE *outputFile = NULL;
		try
		{
			if (!(outputFile = fopen(fileName, "w"))) return 1;
			fprintf(outputFile, "%d\t%d\t%d\n", edgesCount * 2, 0, 0);
			fclose(outputFile);

			return 0;
		}
		catch (exception ex)
		{
			if (outputFile != NULL)
				fclose(outputFile);
			return 1;
		}
	}

	// Write all mesh files procedure
	int WriteMesh3D(vector<double> *meshes1D, vector<double> *meshesTemplate1D, vector<int> &materialNumbers, set<pair<int, int>> &edges, vector<vector<int>> &elementsByEdges)
	{
		vector<double> mesh1D;
		int nx, ny, nz;
		int edgesCount;

		nx = meshes1D[0].size();
		ny = meshes1D[1].size();
		nz = meshes1D[2].size();
		edgesCount = Mesh3D::CalculateEdgesCount(meshes1D);

		if (WriteInftry("inftry.dat", nx, ny, nz) != 0) return 1;
		if (WriteXyz("xyz.dat", meshes1D) != 0) return 1;
		if (WriteXyz("xyzsh.dat", meshes1D) != 0) return 1;
		if (WriteXyz("xyz0.dat", meshesTemplate1D) != 0) return 1;
		if (WriteNver("nver.dat", nx, ny, nz) != 0) return 1;
		if (WriteIntArrayBinary("nvkat.dat", materialNumbers) != 0) return 1;
		if (WriteElemNeib("elem_neib", nx, ny, nz) != 0) return 1;
		if (WriteTSize3DNode("tsize3d.dat") != 0) return 1;
		if (WriteTSize3DEdge("tsize3d_.dat", edgesCount) != 0) return 1;
		if (Write3DMeshRegular("3dmeshregular", meshesTemplate1D) != 0) return 1;
		if (WriteRegular("regular", nx, ny, nz) != 0) return 1;
		if (WriteL13D("L13d.dat", nx, ny, nz) != 0) return 1;
		if (WriteNodesForEdges("nodesforedges.dat", edges) != 0) return 1;
		if (WriteEdges("edges.dat", elementsByEdges) != 0) return 1;
		if (WriteIg3D("ig3d_.dat") != 0) return 1;
		if (WriteJg3D("jg3d_.dat") != 0) return 1;
		if (WriteGg3D("gg3d_.dat") != 0) return 1;
		if (WriteKuSlau2("kuslau2", edgesCount) != 0) return 1;

		return 0;
	}



	
	// Write mesh file with nodes and cells count
	int WriteInf2tr(char *file_name, int nx, int ny, int bc1NodesCount)
	{
		FILE *file_out = NULL;
		try
		{
			if (!(file_out = fopen(file_name, "w"))) return 1;

			int nodesCount = nx * ny;
			int cellsCount = (nx - 1) * (ny - 1);

			fprintf(file_out, "islau= %d", 0);
			fprintf(file_out, " indku1= %d", 0);
			fprintf(file_out, " indfro= %d\n", 0);
			fprintf(file_out, "kuzlov= %d", nodesCount);
			fprintf(file_out, "   ktr= %d", cellsCount);
			fprintf(file_out, "    kt1= %d", bc1NodesCount);
			fprintf(file_out, "   kreb2= %d", 0);
			fprintf(file_out, "   kreb3= %d\n", 0);
			fprintf(file_out, "kisrr1= %d", 2);
			fprintf(file_out, " kisrr2= %d", 2);
			fprintf(file_out, " kisrr3= %d", 2);
			fprintf(file_out, "   kbrsr= %d\n", 8);
			fprintf(file_out, "kreb4= %d\n", 0);
			fclose(file_out);

			return 0;
		}
		catch (exception ex)
		{
			if (file_out != NULL)
				fclose(file_out);
			write_to_log("Error : WriteInf2tr : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Write mesh nodes coordinates
	int WriteRZDat(char *file_name, vector<double> *mesh)
	{
		FILE *file_out = NULL;
		try
		{
			double point[2];

			if (!(file_out = fopen(file_name, "wb"))) return 1;

			auto &meshR = mesh[0];
			auto &meshZ = mesh[1];

			for (int iz = 0; iz < meshZ.size(); iz++)
			{
				point[1] = meshZ[iz];
				for (int ir = 0; ir < meshR.size(); ir++)
				{
					point[0] = meshR[ir];
					fwrite(point, sizeof(double), 2, file_out);
				}
			}

			fclose(file_out);
			return 0;
		}
		catch (exception ex)
		{
			if (file_out != NULL)
				fclose(file_out);
			write_to_log("Error : WriteRZ : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Write mesh nodes count
	int WriteRZTxt(char *file_name, vector<double> *mesh)
	{
		FILE *file_out = NULL;
		try
		{
			if (!(file_out = fopen(file_name, "w"))) return 1;

			fprintf(file_out, "%d\t%d\n", mesh[0].size(), mesh[1].size());

			fclose(file_out);
			return 0;
		}
		catch (exception ex)
		{
			if (file_out != NULL)
				fclose(file_out);
			write_to_log("Error : WriteRZ : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Write mesh cells described with node numbers
	int WriteNvtr(char *file_name, int nr, int nz)
	{
		FILE *file_out = NULL;
		try
		{
			int element[6];
			element[4] = 0;
			element[5] = 1;

			if (!(file_out = fopen(file_name, "wb"))) return 1;

			int node;
			for (int iz = 0; iz < nz - 1; iz++)
			{
				node = iz * nr + 1;
				element[2] = node;
				element[3] = node + 1;
				element[0] = node + nr;
				element[1] = node + nr + 1;
				for (int ir = 0; ir < nr - 1; ir++)
				{
					fwrite(element, sizeof(int), 6, file_out);

					element[0]++;
					element[1]++;
					element[2]++;
					element[3]++;
				}
			}
			fclose(file_out);

			return 0;
		}
		catch (exception ex)
		{
			if (file_out != NULL)
				fclose(file_out);
			write_to_log("Error : WriteNvtr : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Write mesh file with nodes and cells count
	int WriteCurrentVal(char *file_name, double value)
	{
		FILE *file_out = NULL;
		try
		{
			if (!(file_out = fopen(file_name, "w"))) return 1;
			fprintf(file_out, "%.13e\n", value);
			fclose(file_out);

			return 0;
		}
		catch (exception ex)
		{
			if (file_out != NULL)
				fclose(file_out);
			write_to_log("Error : WriteCurrentVal : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	// Write all mesh files procedure
	int WriteMesh2D(vector<double> *meshes1D, vector<int> &materialNumbers, int index)
	{
		char fileName[2048];
		int nx, ny;

		nx = meshes1D[0].size();
		ny = meshes1D[1].size();

		auto boundaryNodesNumbers = Mesh2D::GetFirstBoundaryConditionNodesNumbers(nx, ny, index);

		sprintf(fileName, "inf2tr.dat%d", index);
		if (WriteInf2tr(fileName, nx, ny, boundaryNodesNumbers.size()) != 0) return 1;

		sprintf(fileName, "rz.dat%d", index);
		if (WriteRZDat(fileName, meshes1D) != 0) return 1;

		sprintf(fileName, "nvtr.dat%d", index);
		if (WriteNvtr(fileName, nx, ny) != 0) return 1;

		sprintf(fileName, "rz.txt%d", index);
		if (WriteRZTxt(fileName, meshes1D) != 0) return 1;

		sprintf(fileName, "nvkat2d.dat%d", index);
		if (WriteIntArrayBinary(fileName, materialNumbers) != 0) return 1;

		sprintf(fileName, "r.dat%d", index);
		if (WriteDoubleArrayBinary(fileName, meshes1D[0]) != 0) return 1;

		sprintf(fileName, "z.dat%d", index);
		if (WriteDoubleArrayBinary(fileName, meshes1D[1]) != 0) return 1;

		sprintf(fileName, "tsize.dat%d", index);
		if (WriteTSize3DNode(fileName) != 0) return 1;

		sprintf(fileName, "l1.dat%d", index);
		if (WriteIntArrayBinary(fileName, boundaryNodesNumbers) != 0) return 1;
		
		sprintf(fileName, "currentval%d", index);
		if (WriteCurrentVal(fileName, 1.0) != 0) return 1;
		
		
		return 0;
	}

	// Write material parameters files
	int WriteMaterials3D(vector<Model::Material> &materials)
	{
		if (WriteProperty3D("Sig3d", materials, Model::PhysicalValueName::SigmaH, Model::SigmaN) != 0)
		{
			write_to_log("Could not write 'Sig3d'\n");
			return 1;
		}

		if (WriteProperty3D("dpr3D", materials, Model::PhysicalValueName::Eps, Model::PhysicalValueName::Eps) != 0)
		{
			write_to_log("Could not write 'dpr3D'\n");
			return 1;
		}

		if (WriteProperty3D("mu3D", materials, Model::PhysicalValueName::Mu, Model::PhysicalValueName::Mu) != 0)
		{
			write_to_log("Could not write 'mu3D'\n");
			return 1;
		}

		if (WriteMtr3D2D("mtr3D2D", materials) != 0)
		{
			write_to_log("Could not write 'mtr3D2D'\n");
			return 1;
		}

		return 0;
	}

	// Write material parameters files
	int WriteMaterials2D(vector<Model::Material> &materials)
	{
		if (WriteProperty2D("Sigma", materials, Model::PhysicalValueName::SigmaH) != 0)
		{
			write_to_log("Could not write 'Sigma'\n");
			return 1;
		}

		if (WriteProperty2D("SigmaZ", materials, Model::PhysicalValueName::SigmaH) != 0)
		{
			write_to_log("Could not write 'SigmaZ'\n");
			return 1;
		}

		if (WriteProperty2D("mu", materials, Model::PhysicalValueName::Mu) != 0)
		{
			write_to_log("Could not write 'mu'\n");
			return 1;
		}

		if (WriteProperty2D("dpr", materials, Model::PhysicalValueName::Eps) != 0)
		{
			write_to_log("Could not write 'dpr'\n");
			return 1;
		}

		return 0;
	}


	int WriteCurrentFunction(char *file_name, vector<double> &mesh)
	{
		FILE *file_out = NULL;
		try
		{
			if (!(file_out = fopen(file_name, "w"))) return 1;

			for (size_t i = 0; i < mesh.size(); i++)
			fprintf(file_out, "%.13e\t%.13e\n", mesh[i], i == 0 ? 1.0 : 0.0);

			fclose(file_out);
			return 0;
		}
		catch (exception ex)
		{
			if (file_out != NULL)
				fclose(file_out);
			write_to_log("Error : WriteCurrentFunction : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	int WriteDeltaFunction(char *file_name, vector<double> &mesh)
	{
		FILE *file_out = NULL;
		try
		{
			if (!(file_out = fopen(file_name, "w"))) return 1;

			fprintf(file_out, "0\n%d\n", mesh.size());
			for (size_t i = 0; i < mesh.size(); i++)
				fprintf(file_out, "%.13e\t%.13e\n", mesh[i], i == 0 ? 1.0 : 0.0);

			fclose(file_out);
			return 0;
		}
		catch (exception ex)
		{
			if (file_out != NULL)
				fclose(file_out);
			write_to_log("Error : WriteDeltaFunction : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	int WriteInfite0(char *file_name, vector<double> &mesh)
	{
		FILE *file_out = NULL;
		try
		{
			if (!(file_out = fopen(file_name, "w"))) return 1;

			fprintf(file_out, "    ktime=  %d;  ntstop=  %d;   kiter=    1;   ntime=    1;   niter=    1;\n", mesh.size(), mesh.size());
			fprintf(file_out, "   kprogm=    1;  kpropi=    1;  kpropt=   -1;  kitrel=  250;              ;\n");
			fprintf(file_out, "       u0=   0.00000e+000\n");
			fprintf(file_out, " T I M E :\n");

			for (int i = 0; i < mesh.size(); i += 5)
				if (i < mesh.size() - 4)
				{
					for (int j = 0; j < 5; j++)
						fprintf(file_out, " %.13e;", mesh[i + j]);
					fprintf(file_out, "\n");
				}
				else
				{
					for (; i < mesh.size(); i++)
						fprintf(file_out, " %.13e;", mesh[i]);
					fprintf(file_out, "\n");
				}

			fclose(file_out);
			return 0;
		}
		catch (exception ex)
		{
			if (file_out != NULL)
				fclose(file_out);
			write_to_log("Error : WriteRZ : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}
}