
//============================================ How to run calculation ===================================											  
To run calculation follow the steps (do not use long paths to the directory for calculation!!!):
1. Depending on the computational schemes you want to use: build all programs from folder 'CODE' and go to item 2 
   or use a ready-made folder for calculation  with already assembled programs (folder 'Folders for Calculations') and go to item 7.
2. Create empty directory for calculation running (let's denote it as Home directory).
3. Create 'Modules' directory inside of Home directory.
4. Put all calculation programs (BuildKUForHphi.exe, BuildNormVFor3D.exe, CalcNonstat2D.exe, CalcNonstat2dAr.exe, CalcNonstat2dAx.exe, CalcNonstat2dHphi.exe, CalcNonstat3DLine.exe, CalcStat2D.exe, CalcStat3D_3T_A.exe, CalcStat3D_3T_V.exe, CalcStat3D_A.exe, CalcStat3D_V.exe, CorrectTimeMesh.exe, GeoCalculatorSP.exe, OutputNonstat2D.exe, OutputNonstat2DOnlyRec.exe, OutputStat2D.exe, RegularMeshBuilder.exe, SrsMesh2D.exe, SummatorEDMN.exe, SumNonstat2D3D.exe, SumNonstat2DOnly.exe, SumNonstat3DOnly.exe, SumStat2D3D.exe, TimeDirectSolver.exe, UnloadAnomalNonstat.exe) to 'Modules'
5. (If needed) Put required dll files to ' Modules' directory.
6. Put CalcStarter.exe to Home directory.
7. Put all prepared text files (gen, layers, objects, rec, recPerGen, settings.cfg) to Home directory.
8. Run CalcStarter.exe and wait for it to complete.
9. Look for result file (edmnall0_N) in results directory 'Calculations' inside of Home directory.

All dimensions are in meters. In the calculations, the transmitter current is taken equal to 1 A.

//===================The following files contain model geoemtry and resistivity description. ====================================================================

1. layers
This file contains the resistivity of layers defined in the model. Do not include air in this file. The air layer will be added automatically

File format:
<Layers count>
For each layer:
<Upper bound Z> <Conductivity>

=============================================================================================
2. objects
This file contains 3D objects (heterogeneties) defined as a parallelepiped.

File format:
<Objects count>
For each polygon:
<X0> <X1> <Y0> <Y1> <Z0> <Z1> <Conductivity>

=============================================================================================
3. gen
This file contains transmitter (AB) lines placements

File format:
For each transmitter line:
<AX> <AY> <AZ> <BX> <BY> <BZ>

=============================================================================================
4. recPerGen
This file describes how many receiver (MN) linesfrom the common list belong to each transmitter line

File format:
For each transmitter line:
<Receiver lines count>

=============================================================================================
5. rec
This file contains common list of receiver (MN) lines placements

File format:
<Receiver lines count>
For each line:
<MX> <MY> <MZ> <NX> <NY> <NZ>

=============================================================================================
6. settings.cfg
This file contains mesh settings and some calculation settings

File format:

Mesh 3D section:
<Initial step X, Y> <Initial step Z>
<Sparce coefficient X, Y> <Sparce coefficient Z>
<Gap between receivers bounding box and sparced mesh area X, Y> <Gap between receivers bounding box and sparced mesh area Z>
<Distance from receivers bounding box+gap to the boundary of the calculation domain X, Y> <Distance from receivers bounding box+gap to the boundary of the calculation domain Z>
<Number of refinements around recevers>

<Approximation points count for receiver line>
<Approximation dipoles count for transmitter line>
    
Mesh 2D section:
<Initial step R> <Initial step Z>
<Sparce coefficient R> <Sparce coefficient Z>
<Gap between receivers bounding box and sparced mesh area R> <Gap between receivers bounding box and sparced mesh area Z>
<Distance from receivers bounding box+gap to the boundary of the calculation domain R> <Distance from receivers bounding box+gap to the boundary of the calculation domain Z>

Mesh time section:
<Initial step>
<Ssparce coefficient>
<Last time in time mesh>

Calculations section:
< 1 if calculate directly with no field decomposition (primary-secondary), 0 otherwise>
<Threads count for parallel calculation>


//============================================ Results file ==============================================
'edmnall0_N' file contains the signals in all receivers for transmitter line position with number N (according to ‘gen’ file).
File consists of small tables for each receiver (according to receivers from 'rec' file that belong to transmitter with number N as written in 'recPerGen' file).
The table format is followng:
<title line 1 (ignore it)> 
<title line 2 (ignore it)> 
<title line 3 (ignore it)>
<title line 4 (ignore it)>
<title line 5 (ignore it)>
Then a column with times (ms) and three columns with results: 
dVp (mV) — primary field
dVs (mV)  — secondary field
dV (mV)  — total field

To compare with the results presented in the article for DNME technology, for calculating the D2U signal, it is necessary to subtract the results obtained in the 2nd and 3rd receivers.

