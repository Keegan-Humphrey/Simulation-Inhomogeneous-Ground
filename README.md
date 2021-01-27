# Simulation-Inhomogeneous-Ground

This build (as it stands) looks generate training data for Classify.py to recognize the Ubox being simulated.  Parameters of the simulation can be modified by editing Filing.py. Relief data is fed to the simulation through the files Relief.cc and Relief.hh. 

To run:

$ cd ./Source

$ python Filing.py

$ bash Compile

- the second command creates the BASH and GEANT scripts for the simulations
- the third command runs the BASH scripts to compile the source files, run the executables they create, and analyse the resulting data files

New Directories:

Source - contains the unaltered DetectorConstruction.cc file, Filing.py which determines the parameters of the simulation and alters DetectorConstruction.cc, and Triangles_logic.cc which was the first attempt at making a continuous landscape (The logic will be added to DetectorConstruction.cc once current issues are resolved).

Executables - contains all compiled executables in the run, and a script to generate the macro file with a random 8 digit seed.

RowData - contains all RowData.out files for the run, as well as ClusterAnalysis.py and Plot.py which process the data and plots it.

Training - end of the data flow. Visulaise the data using Plot.py and train a SGD algorithm to distinguish between runs with and without the Ubox present using Classify.py.




