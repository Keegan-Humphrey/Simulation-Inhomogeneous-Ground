# Simulation-Inhomogeneous-Ground


## Summary:

This application runs simulations, analysis, and visualisation for Muon Tomography using a portable plastic scintillator detector. Parameters of the simulation can be modified by editing Filing.py. Relief data is fed to the simulation through the files Relief.cc and Relief.hh (which are updated with Filing.py). Scripts are available to probe available information and process data in a variety of ways. 



## To run:
```
$ cd ./Source

$ python Filing.py

$ bash Compile

$ bash Run #(optional)
```
- The second command creates the BASH and GEANT files for the simulations
- The third command runs the BASH scripts to compile the source files (if chain == True - in Filing.py - it will also run the executables they create, and analyse the resulting data files)
- The fourth command is only needed if the chain boolean option is False in Filing.py (to visualise simulation geometry before running)



## Directory Contents:

###### Source 
- Contains the unaltered DetectorConstruction.cc file, Filing.py (which determines the parameters of the simulation and alters DetectorConstruction.cc), and contains the executables to run the simulations (Compile and Run).

###### Executables 
- Contains all compiled executables in the run, and a script to generate the macro file with a random 8 digit seed. It also contains visualisation tools (HepRApp.jar and vis.mac) which can be used to view the geometry of an executable using the following command.
```
$ cd ./Executables

$ ./date_and_time/executable_name vis.mac
```
- This creates .heprep files which can be visualised by running the HepRApp.jar executable and opening the generated files in the application.

###### RowData 
- Contains all RowData.out files for the run, as well as ClusterAnalysis.py and Plot.py which process the data and plots it. It also contains TestError.py which will introduce error to the data and show the change in the predicted positions of the muons.

###### Training 
- End of the data flow. Visulaise the data using Plot.py, train a SGD algorithm to distinguish between runs with and without the Ubox present using Classification.py, or run the tracking and clustering algorithms on a dataset to identify and visualise locations of objects.

