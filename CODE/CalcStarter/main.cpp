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

#include <stdio.h>
#include <stdlib.h>
#include <filesystem>
#include <vector>
#include <string>
#include <Windows.h>
#include <direct.h>

using namespace std;

// Calculation settings
struct Settings
{
public:
	double steps3D[2];    // Mesh initial steps
	double sparse3D[2];   // Sparce coefficient
	double farBound3D[2]; // Distance from receivers bounding box+gap to the boundary of the calculation domain
	double gapFromReceivers3D[2]; // Gap between receivers bounding box and sparced mesh area
	int receiversRefinementsCount; // Mesh refinements around receivers

	int numberOfMNPoints;
	int numberOfABDipoles;

	double steps2D[2];    // Mesh initial steps
	double sparse2D[2];   // Sparce coefficient
	double farBound2D[2]; // Distance from receivers bounding box+gap to the boundary of the calculation domain
	double gapFromReceivers2D[2]; // Gap between receivers bounding box and sparced mesh area

	double timeStep = 1e-6;    // Time mesh initial steps
	double timeSparse = 1.1;   // Time mesh sparce coefficient
	double timeLast = 1.43e+4; // last time in time mesh

	bool directCalc = false; // calculate directly with no field decomposition (primary-secondary)

	int threadsCount = 1; // Threads count for parallel calculation

	Settings() { }
};

// Electrode position
struct Electrode
{
public:

	double coordinates[3]; // x, y, z
};

FILE *log_file = NULL;

// Write message to log
void write_to_log(const char *str)
{
	printf("%s", str);
	if (log_file != NULL)
	{
		fprintf(log_file, "%s", str);
		fflush(log_file);
	}
}

// Open file for reading
int open_file_r(const char *file_name, FILE **file_stream)
{
	char buf[2048];

	if (!((*file_stream) = fopen(file_name, "r")))
	{
		sprintf(buf, "Error : Could not read file '%s'\n", file_name);
		write_to_log(buf);
		return 1;
	}

	return 0;
}

// Open file for writing 
int open_file_w(const char *file_name, FILE **file_stream)
{
	char buf[2048];

	if (!((*file_stream) = fopen(file_name, "w")))
	{
		sprintf(buf, "Error : Could not write file '%s'\n", file_name);
		write_to_log(buf);
		return 1;
	}

	return 0;
}

// Open log file
int open_log(const char *file_name)
{
	if (open_file_w(file_name, &log_file) != 0)
	{
		printf("Error : could not open file '%s'\n", file_name);
		return 1;
	}

	return 0;
}


// Run Windows executable file
int ExecuteExe(const char *cmdline, const char *workdir)
{
	int retCode;
	STARTUPINFOA si;
	ZeroMemory(&si, sizeof(si));
	si.cb = sizeof(si);

	PROCESS_INFORMATION pi;
	ZeroMemory(&pi, sizeof(pi));

	CreateProcessA(NULL, (LPSTR)cmdline, NULL, NULL, FALSE, CREATE_NO_WINDOW, NULL, workdir, &si, &pi);
	WaitForSingleObject(pi.hProcess, INFINITE);
	GetExitCodeProcess(pi.hProcess, (LPDWORD)&retCode);

	CloseHandle(pi.hProcess);
	CloseHandle(pi.hThread);
	return retCode;
}


// Write single string file
int WriteSingleStringFile(const char *file_name, const char *str)
{
	FILE *file_out = NULL;
	if (open_file_w(file_name, &file_out) != 0) return 1;

	fprintf(file_out, "%s\n", str);

	fclose(file_out);
	return 0;

}

// Write single int file
int WriteSingleIntFile(const char *file_name, const int val)
{
	FILE *file_out = NULL;
	if (open_file_w(file_name, &file_out) != 0) return 1;

	fprintf(file_out, "%d\n", val);

	fclose(file_out);
	return 0;

}

// Write single string file
int WriteLcTxt(const char *file_name)
{
	FILE *file_out = NULL;
	if (open_file_w(file_name, &file_out) != 0) return 1;

	fprintf(file_out, "3\n");
	fprintf(file_out, "1e-3\n");
	fprintf(file_out, "1e-2\n");
	fprintf(file_out, "\n");
	fprintf(file_out, "3\n");
	fprintf(file_out, "1e-3\n");
	fprintf(file_out, "1e-2\n");
	fprintf(file_out, " \n");
	fprintf(file_out, "\n");
	fprintf(file_out, "\n");
	fprintf(file_out, "0\n");
	fprintf(file_out, "1e-3\n");
	fprintf(file_out, "1e-3\n");
	fprintf(file_out, "0.5\n");
	fprintf(file_out, "0.5\n");
	fprintf(file_out, "0.0\n");
	fprintf(file_out, "\n");
	fprintf(file_out, "0\n");
	fprintf(file_out, "1e-3\n");
	fprintf(file_out, "1e-3\n");
	fprintf(file_out, "0.5\n");
	fprintf(file_out, "0.5\n");
	fprintf(file_out, "0.5\n");
	fprintf(file_out, "\n");
	fprintf(file_out, "1\n");
	fprintf(file_out, "1e+30\n");
	fprintf(file_out, "1e-3\n");

	fclose(file_out);
	return 0;

}

// Write TimeIntervalForPrint
int WriteTimeIntervalForPrint(const char *file_name, Settings &settings)
{
	FILE *file_out = NULL;
	if (open_file_w(file_name, &file_out) != 0) return 1;

	fprintf(file_out, "%.13e\t%.13e\n", settings.timeStep * 5.0, settings.timeLast * 1.5);

	fclose(file_out);
	return 0;

}

// Write mes settings file
int WriteMeshSettings(const char *file_name, Settings &settings)
{
	FILE *file_out = NULL;
	if (open_file_w(file_name, &file_out) != 0) return 1;
	char c;

	fprintf(file_out, "%.13e\t%.13e\n", settings.steps3D[0], settings.steps3D[1]);
	fprintf(file_out, "%.13e\t%.13e\n", settings.sparse3D[0], settings.sparse3D[1]);
	fprintf(file_out, "%.13e\t%.13e\n", settings.gapFromReceivers3D[0], settings.gapFromReceivers3D[1]);
	fprintf(file_out, "%.13e\t%.13e\n", settings.farBound3D[0], settings.farBound3D[1]);
	fprintf(file_out, "%d\n", settings.receiversRefinementsCount);

	fprintf(file_out, "%.13e\t%.13e\n", settings.steps2D[0], settings.steps2D[1]);
	fprintf(file_out, "%.13e\t%.13e\n", settings.sparse2D[0], settings.sparse2D[1]);
	fprintf(file_out, "%.13e\t%.13e\n", settings.gapFromReceivers2D[0], settings.gapFromReceivers2D[1]);
	fprintf(file_out, "%.13e\t%.13e\n", settings.farBound2D[0], settings.farBound2D[1]);

	fprintf(file_out, "%.13e\n", settings.timeStep);
	fprintf(file_out, "%.13e\n", settings.timeSparse);
	fprintf(file_out, "%.13e\n", settings.timeLast);

	fclose(file_out);
	return 0;

}


// Read application settings file
int ReadSettings(const char *file_name, Settings &settings)
{
	FILE *file_in = NULL;
	if (open_file_r(file_name, &file_in) != 0) return 1;
	int tmpi;

	fscanf(file_in, "%lf%lf", settings.steps3D, settings.steps3D + 1);
	fscanf(file_in, "%lf%lf", settings.sparse3D, settings.sparse3D + 1);
	fscanf(file_in, "%lf%lf", settings.gapFromReceivers3D, settings.gapFromReceivers3D + 1);
	fscanf(file_in, "%lf%lf", settings.farBound3D, settings.farBound3D + 1);
	fscanf(file_in, "%d", &settings.receiversRefinementsCount);

	fscanf(file_in, "%d", &settings.numberOfMNPoints);
	fscanf(file_in, "%d", &settings.numberOfABDipoles);	

	fscanf(file_in, "%lf%lf", settings.steps2D, settings.steps2D + 1);
	fscanf(file_in, "%lf%lf", settings.sparse2D, settings.sparse2D + 1);
	fscanf(file_in, "%lf%lf", settings.gapFromReceivers2D, settings.gapFromReceivers2D + 1);
	fscanf(file_in, "%lf%lf", settings.farBound2D, settings.farBound2D + 1);

	fscanf(file_in, "%lf", &settings.timeStep);
	fscanf(file_in, "%lf", &settings.timeSparse);
	fscanf(file_in, "%lf", &settings.timeLast);

	fscanf(file_in, "%d", &tmpi);
	settings.directCalc = tmpi == 1;

	fscanf(file_in, "%d", &settings.threadsCount);
	fclose(file_in);

	return 0;
}


// Prepare directory for meshing
int PrepareDirectoryForMeshing(const char *path, Settings &settings)
{
	char fileName[2048];

	sprintf(fileName, "%s/objects", path); CopyFileA("objects", fileName, false);
	sprintf(fileName, "%s/z_sig_2d", path); CopyFileA("layers", fileName, false);
	sprintf(fileName, "%s/sours", path); CopyFileA("gen", fileName, false);
	sprintf(fileName, "%s/xyzmn", path); CopyFileA("rec", fileName, false);
	sprintf(fileName, "%s/lin", path); CopyFileA("recPerGen", fileName, false);

	sprintf(fileName, "%s/RegularMeshBuilderSettings.cfg", path);
	if (WriteMeshSettings(fileName, settings) != 0)
		return 1;

	return 0;
}

// Prepare calculations directory
int PrepareDirectoryCalculation(const char *path, const char *modulesPath, Settings &settings)
{
	char fileName[2048];
	char fileName2[2048];
	char dirName[2048];


	sprintf(fileName, "%s/vcomp90.dll", path); 
	sprintf(fileName2, "%s/vcomp90.dll", modulesPath);	
	CopyFileA(fileName2, fileName, false);

	sprintf(fileName, "%s/Microsoft.VC90.OpenMP.manifest", path);
	sprintf(fileName2, "%s/Microsoft.VC90.OpenMP.manifest", modulesPath);
	CopyFileA(fileName2, fileName, false);

	if (settings.directCalc)
	{
		sprintf(fileName, "%s/fdirect", path);
		if (WriteSingleStringFile(fileName, "1") != 0)
			return 1;

		sprintf(dirName, "%s/stationary", path);
		_mkdir(dirName);
		PrepareDirectoryForMeshing(dirName, settings);
	}

	sprintf(fileName, "%s/dminsrc", path);
	if (WriteSingleStringFile(fileName, "0.1") != 0)
		return 1;

	sprintf(fileName, "%s/igroup", path);
	if (WriteSingleStringFile(fileName, "1") != 0)
		return 1;

	sprintf(fileName, "%s/imp", path);
	if (WriteSingleStringFile(fileName, "20000\n20000\n0\n0") != 0)
		return 1;

	sprintf(fileName, "%s/inv_flags", path);
	if (WriteSingleStringFile(fileName, "0\n0\n1\n1\n0\n") != 0)
		return 1;
	
	sprintf(fileName, "%s/lc.txt", path);
	if (WriteLcTxt(fileName) != 0)
		return 1;

	sprintf(fileName, "%s/n_neib", path);
	if (WriteSingleStringFile(fileName, "3") != 0)
		return 1;

	sprintf(fileName, "%s/nSrsRazb", path);
	if (WriteSingleStringFile(fileName, "500") != 0)
		return 1;

	sprintf(fileName, "%s/numberofabdipoles", path);
	if (!filesystem::exists(std::string(fileName)))
		if (WriteSingleStringFile(fileName, std::to_string(settings.numberOfABDipoles).c_str()) != 0)
			return 1;

	sprintf(fileName, "%s/numberofmnpoints", path);
	if(!filesystem::exists(std::string(fileName)))
		if (WriteSingleStringFile(fileName, std::to_string(settings.numberOfMNPoints).c_str()) != 0)
			return 1;

	sprintf(fileName, "%s/recvsb", path);
	if (WriteSingleStringFile(fileName, "0") != 0)
		return 1;

	sprintf(fileName, "%s/scalingZ", path);
	if (WriteSingleStringFile(fileName, "1") != 0)
		return 1;

	sprintf(fileName, "%s/tcff", path);
	if (WriteSingleStringFile(fileName, "2.0") != 0)
		return 1;

	sprintf(fileName, "%s/xyzVectorB", path);
	if (WriteSingleStringFile(fileName, "0") != 0)
		return 1;

	sprintf(fileName, "%s/TimeIntervalForPrint", path);
	if (WriteTimeIntervalForPrint(fileName, settings) != 0)
		return 1;

	sprintf(fileName, "%s/nthreads.txt", path);
	if (WriteSingleIntFile(fileName, settings.threadsCount) != 0)
		return 1;

	return 0;
}

// Prepare modules directory
int PrepareDirectoryModules(const char *path)
{
	char fileName[2048];
	char dirName[2048];
	char fileName2[2048];

	sprintf(fileName, "%s/p.bat", path);
	
	GetCurrentDirectoryA(2048, dirName);
	sprintf(fileName2, "%s/%s/RegularMeshBuilder.exe", dirName, path);

	if (WriteSingleStringFile(fileName, fileName2) != 0)
		return 1;

	return 0;
}

// Prepare calculation directories
int PrepareDirectory(Settings &settings, const char *meshPath, const char *calculationPath, const char *modulesPath)
{
	_mkdir(calculationPath);

	PrepareDirectoryForMeshing(calculationPath, settings);
	PrepareDirectoryCalculation(calculationPath, modulesPath, settings);
	PrepareDirectoryModules(modulesPath);

	return 0;
}

// Run calculation
int RunCalculation(const char *meshPath, const char *calculationPath, const char *modulesPath)
{
	char fileName[2048];

	write_to_log("Starting RegularMeshBuilder.exe\n");
	sprintf(fileName, "%s/RegularMeshBuilder.exe", modulesPath);
	if (ExecuteExe(fileName, calculationPath) != 0) return 1;
	
	write_to_log("Starting GeoCalculatorSP.exe\n");
	sprintf(fileName, "%s/GeoCalculatorSP.exe", modulesPath);
	if (ExecuteExe(fileName, calculationPath) != 0) return 1;

	return 0;
}

// Calling common procedures
int MainProcedure()
{
	Settings settings;

	const char *meshPath = "Mesh";
	const char *calculationPath = "Calculations";
	const char *modulesPath = "Modules";

	write_to_log("Reading settings.cfg\n");
	if (ReadSettings("settings.cfg", settings) != 0)
		return 1;

	write_to_log("Preparing directories\n");
	if (PrepareDirectory(settings, meshPath, calculationPath, modulesPath) != 0)
		return 1;

	write_to_log("Running calculation\n");
	if (RunCalculation(meshPath, calculationPath, modulesPath) != 0)
		return 1;

	return 0;
}


int main()
{
	char buf[2048];

	if (open_log("CalcStarter.log") != 0)
	{
		printf("Cound not open CalcStarter.log\n");
		return 1;
	}

	int status = MainProcedure();
	sprintf(buf, "Status = %d\n", status);
	write_to_log(buf);

	fprintf(log_file, "All done\n");
	fclose(log_file);
	return status;
}