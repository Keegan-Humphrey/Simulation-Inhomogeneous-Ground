# Simulation-Inhomogeneous-Ground

This build (as it stands) looks to visualise the landscape above the detector by comparing it to homogeneous ground. Parameters of the simulation can be modified by editing Filing.py. Relief data is fed to the simulation through the files Relief.cc and Relief.hh respectively. 

To run:

$ cd ./Source

$ python Filing.py

$ chmod 744 Compile

$ ./Compile


New Directories:

Source - contains altered DetectorConstruction.cc file, bash scripts for running the simulation, and Filing.py which determines the parameters of the simulation

Executables - contains all compiled executables in the run

RowData - contains all RowData.out files for the run, as well as ClusterAnalysis.py and Plot.py which process the data and plots it respectively



