# Simulation-Inhomogeneous-Ground

This build (as it stands) looks to visualise the landscape above the detector by comparing it to homogeneous ground. Parameters of the simulation can be modified by editing Filing.py. Relief data is fed to the simulation through the files Relief.cc and Relief.hh respectively. 

To run:

$ cd ./Source

$ python Filing.py

$ bash Compile


New Directories:

Source - contains the unaltered DetectorConstruction.cc file, Filing.py which determines the parameters of the simulation and alters DetectorConstruction.cc, and Triangles_logic.cc which was the first attempt at making a continuous landscape (The logic will be added to DetectorConstruction.cc once current issues are resolved).

Executables - contains all compiled executables in the run, and a script to generate the macro file with a random 8 digit seed

RowData - contains all RowData.out files for the run, as well as ClusterAnalysis.py and Plot.py which process the data and plots it respectively



