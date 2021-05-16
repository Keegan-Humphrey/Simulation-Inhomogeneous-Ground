#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 10:51:15 2020

@author: keegan
"""
import numpy as np
from joblib import dump, load
import random
import string
from time import time, asctime




#%%
# Free parameters for simulation (and landscape geometry)


alpha = [0] #[deg] x, y, z axis rotations
beta = [0] #[deg]
gamma = [0] #[deg]

Det_Position = [[0,0]] #[cm] (x,y)
Depth = [30] #[m] z

Seperation = [25] #[cm]

Events = int(1e4)
dump(Events, 'Events.joblib')

Ubox_Material_real = ['BMF', 'Air', 'GroundMat', 'Vacuum'][0]
Ubox_Material_sky = ['BMF', 'Air', 'GroundMat', 'Vacuum'][2]
Ubox_dims = [1, 1, 1] # [m] (x, y, z) half dimensions
Ubox_location = [[0,0,10]] #[m] (x,y,-z) from surface centre
        
Inhomogeneous = True # whether or not landscape is added
                     # -> depricated, no longer affects simulation
                     
Landscape_material = ['BMF', 'GroundMat','Air','Vacuum'][1]
GroundMaterial = ['BMF', 'GroundMat','Air','Vacuum'][1]

chain = True # whether or not bash scripts are chained
            # set as false to check visualisation before running
            # a large data set

transverse_scale = 10 #[m] the seperation between adjascent points in x and y

Elev = np.array([[715.77, 710.52, 708.86, 707.30, 703.07, 702.89, 703.36, 703.55], \
                 [714.05, 710.47, 704.03, 704.17, 703.42, 703.05, 702.86, 702.98], \
                 [711.74, 708.42, 703.88, 703.57, 702.96, 702.60, 702.17, 703.02], \
                 [710.26, 705.41, 704.11, 703.72, 703.17, 702.85, 702.29, 701.65], \
                 [708.63, 704.29, 703.65, 703.62, 703.05, 702.74, 702.18, 701.39], \
                 [706.67, 703.89, 703.28, 703.30, 702.91, 702.79, 702.05, 701.17], \
                 [707.68, 704.84, 704.07, 704.16, 703.40, 702.68, 702.17, 701.52]])

Building_height = 712.00

#Elev[0:2,5:] = Building_height

Z = Elev - np.amin(Elev)

Z = np.round(Z, 1)

Shape_Z = np.shape(Z)


#%%


def WriteSource(fName, Seperation, Depth, X, Y, Real, Alpha, Beta, Gamma, boxloc): 
    '''
    Takes DetectorConstruction.cc file and creates new file with altered parameters.
    Writes in landscape geometry to be read from Relief.cc for Real files only
    
    input:
        fName [str] - new file name
        Seperation [float] - [cm] seperation of detector planes
        Depth [float] - [m] depth of detector
        X [float] - [cm] x position of detector
        Y [float] - [cm] y position of detector
        Real [bool] - True ==> Real, False ==> Sky
        Alpha, Beta, Gamma [float] - [deg] Euler angles of detector
        boxloc [float] - [m] location of the object
    
    ouput:
        None
    '''
    Source = 'DetectorConstruction.cc'
    
    g = open(Source)
    f = open(fName,'w+')
    
    g_lines = g.readlines()
    
    g_lines[51] = 'G4double UD_Det_planes_seperation = {}*cm;\n'.format(Seperation) ####check line numbers 
    g_lines[59] = 'G4double DetectorDepth = {}*m;\n'.format(Depth)
    g_lines[60] = 'G4double UD_Det_centerX = {}*cm;\n'.format(X)
    g_lines[61] = 'G4double UD_Det_centerY = {}*cm;\n'.format(Y)
    
    g_lines[63] = 'G4double UD_Det_alpha = {}*deg; \n'.format(Alpha)
    g_lines[64] = 'G4double UD_Det_beta = {}*deg; \n'.format(Beta)
    g_lines[65] = 'G4double UD_Det_gamma = {}*deg; \n'.format(Gamma)
    
    g_lines[216] = '\tG4LogicalVolume* logicGround = new G4LogicalVolume(solidGround, {}, "HomoGround"); \n'.format(GroundMaterial)
    
    g_lines[227] = '\tG4double UBoxHight = %s*m; \n'%(Ubox_dims[2])
    g_lines[228] = '\tG4VSolid* UBox = new G4Box("UBox", %s*m, %s*m, UBoxHight); \n'%(Ubox_dims[0], Ubox_dims[1])
    
    if Real:
        g_lines[230] = '\tG4LogicalVolume* UBoxLogicVolume = new G4LogicalVolume(UBox, %s, "UBox");\n'%(Ubox_Material_real)
        g_lines[231] = '\t//G4LogicalVolume* UBoxLogicVolume = new G4LogicalVolume(UBox, %s, "UBox");\n'%(Ubox_Material_sky)
    
    else:
        g_lines[230] = '\t//G4LogicalVolume* UBoxLogicVolume = new G4LogicalVolume(UBox, %s, "UBox");\n'%(Ubox_Material_real)
        g_lines[231] = '\tG4LogicalVolume* UBoxLogicVolume = new G4LogicalVolume(UBox, %s, "UBox");\n'%(Ubox_Material_sky)
    
    #g_lines[224] = '\tG4double UboxDepth = {} * m;'.format(boxloc[2])
    g_lines[235] = '\tG4ThreeVector({} * m, {} * m, DetectorDepth / 2 - UBoxHight - {} * m),'.format(*boxloc)
    
    g_lines.insert(43, '#include "Relief.hh" \n \n')
    g_lines.insert(43, '#include "G4Tet.hh" \n \n')
    
    f.writelines(g_lines)
    
    f.close()
    g.close()



def WriteHeader(Data):
    '''
    Write a Geant readable source file to describe the landscape geometry described by Data
    Assumes there already exists an 'include' file called Relief.hh

    input:
        Data [float] - 2 dim array of heights above homogeneous ground
    
    output:
        None
    '''
    Shape = np.shape(Data)
    
    Mins = np.array([[np.amin([Data[i,j],Data[i,j+1],Data[i+1,j],Data[i+1,j+1]]) for j in range(Shape[1]-1)] for i in range(Shape[0]-1)])
    
    g = open('Relief.cc','w+')
    
    # Write Relief array
    g_lines = ['#include <array> \n', \
               '#include "G4SystemOfUnits.hh" \n', \
               '#include "G4UnitsTable.hh" \n \n'] #, \
               #'namespace RData \n', \
               #'{ \n \n']
               
    g_lines.append('G4double Pnt_Sep = {} * m; //seperation of data points \n \n'.format(transverse_scale))
    
    g_lines.append('int iaxis = {}; \n'.format(Shape[0]))
    g_lines.append('int jaxis = {}; \n \n'.format(Shape[1]))
    g_lines.append('G4double max_landscape_height = {} * m; \n \n'.format(np.amax(Z)))
    
    
    
    #paste data ----------------------------------------------------------#
    g_lines.append('\n std::array<std::array<G4double, %s>, %s> Relief = {{ \n'%(Shape[1],Shape[0]))
    
    for i in range(Shape[0]-1):
        Str = '{{'

        for j in range(Shape[1]-1):
            Str += str(Data[i,j]) + ' * m, '

        Str += str(Data[i,Shape[1]-1]) + ' * m}}, \n'

        g_lines.append(Str)

    Str = '{{'

    for j in range(Shape[1]-1):
        Str += str(Data[i,j]) + ' * m, '

    Str += str(Data[i,Shape[1]-1]) + ' * m}} \n'

    g_lines.append(Str)

    g_lines.append('}}; \n \n')



    #paste minima ----------------------------------------------------------#
    g_lines.append('\n std::array<std::array<G4double, %s>, %s> RMins = {{ \n'%(Shape[1]-1,Shape[0]-1))

    for i in range(Shape[0]-2):
        Str = '{{'

        for j in range(Shape[1]-2):
            Str += str(Mins[i,j]) + ' * m, '

        Str += str(Mins[i,Shape[1]-2]) + ' * m}}, \n'

        g_lines.append(Str)

    Str = '{{'

    for j in range(Shape[1]-2):
        Str += str(Mins[i,j]) + ' * m, '
	
    Str += str(Mins[i,Shape[1]-2]) + ' * m}} \n'

    g_lines.append(Str)

    g_lines.append('}}; \n \n')
    #g_lines.append('} \n \n')

    g.writelines(g_lines)

    g.close()
    
    # ----------------------------------------------------------#
    
    f = open('Relief.hh','w+')
    
    f_lines = ['#include <array> \n', \
                '#include "G4SystemOfUnits.hh"  \n', \
                '#include "G4UnitsTable.hh" \n \n', \
                'extern G4double Pnt_Sep; //seperation of data points \n', \
                'extern int iaxis; \n', \
                'extern int jaxis; \n \n', \
                'extern G4double max_landscape_height; \n \n', \
                'extern std::array<std::array<G4double, %s>, %s> Relief; \n \n'%(Shape[1],Shape[0]), \
                'extern std::array<std::array<G4double, %s>, %s> RMins; \n \n'%(Shape[1]-1,Shape[0]-1)]
    
    f.writelines(f_lines)
    
    f.close()
    
    

def WriteCommands(src_name, ex_name, dir_name, row_name, date_dir): 
    '''
    Takes file parameters and generates scripts to compile and 
    to run them on the cluster. 
    
    input:
        src_name [str] - source file name
        ex_name [str] - desired name of executable
        dir_name [str] - name of source file parent directory
        row_name [str] - name of RowData.out file
        date_dir [str] - name of sub directory
        
    output:
        None
    '''
    f = open('Compile','a') #'a' ==> append (so the fn can be iterated)
    g = open('Run','a')
    
    Lines_f = ['mv {} ../../src \n'.format(src_name), \
            'cd ../../src \n', \
            'mv {} DetectorConstruction.cc \n'.format(src_name), \
            'cd ../ \n', \
            'rm -r build \n', \
            'mkdir build \n', \
            'cd build \n', \
            'rm -rf * \n', \
            'cmake ../ \n', \
            'make \n', \
            'mv Muons {} \n'.format(ex_name), \
            'mv {} ../Executables \n'.format(ex_name), \
            'cd ../src \n', \
            'mv DetectorConstruction.cc {} \n'.format(src_name), \
            'mv {} ../{}/{} \n \n'.format(src_name, dir_name,date_dir), \
            'cd ../{}/{} \n'.format(dir_name,date_dir)]
            
    Lines_g = ['python Macro_Seed.py \n', \
               './{} \n'.format(ex_name), \
               'mv RowData.out {} \n'.format(row_name), \
               'mv {} ../RowData/{} \n \n'.format(row_name,date_dir)]
    
    f.writelines(Lines_f)
    g.writelines(Lines_g)
    
    f.close()
    g.close()



def Write_Description(cwd):
    '''
    Writes a .txt file to describe the data generated by the run.
    To be put in the Training folder corresponding to the data. 
    
    input:
        cwd [str] - name of current working directory
        
    output:
        None
    '''
    f = open('DESCRIPTION.txt','w+')
    
    f.writelines(['The Following describes the data in ../{}/ \n\n\n'.format(cwd),
                  'The depth of the detector from the surface minimum is: \t {} m \n\n'.format(Depth),
                  'The x-y position of the detector is: \t {} rad \n\n'.format(Det_Position),
                  'The seperation of the detector planes is: \t {} cm \n\n'.format(Seperation),
                  'The euler angles of the detector are: \t {}, {}, {} deg \n\n'.format(alpha,beta,gamma),
                  'The number of events each run is: \t {} \n\n'.format(Events),
                  'The material of the box is: \t {} \n\n'.format(Ubox_Material_real),
                  'The dimensions of the box are: \t ({}, {}, {}) m \n\n'.format(*Ubox_dims),
                  'The location of the box is: \t {} m \n\n'.format(Ubox_location),
                  'The inhomogeneous landscape {} simulated \n\n'.format(["wasn't",'was'][Inhomogeneous]),
                  'The material of the landscape {} {} \n\n'.format(['would have been','is'][Inhomogeneous],Landscape_material)])
        
    f.close()
    
    

def CreateSource(det_position, depth, seperation):
    '''
    Write detector construction files for all combinations of the elements in
    radius theta and seperation for both real and sky geometries
    Write macro file with random 10 digit seed
    Write bash scripts to compile and run the simulations
    Output the parameters used to to make them available for later analysis
    
    input:
        radius [float] - list of distances from the surface level origin
        theta [float] - list of angles relative to vertical for placement of Ubox
        seperation [float] - list of seperations of detector planes
        
    output:
        None
    '''
    
    current_date_time = asctime().replace(' ','_').replace(':','-')
    
    dump(current_date_time,'time.joblib')

    Write_Description(current_date_time)
    
    f = open('Compile','w+')
    g = open('Run','w+')
    
    f.writelines(['#!/bin/bash \n \n', \
                  'mv Relief.cc ../src \n', \
                  'mv Relief.hh ../include \n \n', \
                  'mkdir ../Training/{} \n'.format(current_date_time), \
                  'mkdir ../Executables/{} \n'.format(current_date_time), \
                  'mkdir ../RowData/{} \n'.format(current_date_time), \
                  'mkdir {} \n \n'.format(current_date_time), \
                  'cp DESCRIPTION.txt ../RowData/{} \n \n'.format(current_date_time), \
                  'mv DESCRIPTION.txt ../Training/{} \n \n'.format(current_date_time), \
                  'mv Parameters.joblib ../RowData \n \n', \
                  'mv Events.joblib ../Executables \n \n', \
                  'mv time.joblib ../RowData \n \n',
                  'mv DC*.cc ./{} \n'.format(current_date_time), \
                  'cd ./{} \n'.format(current_date_time)])
                  
    g.writelines(['#!/bin/bash \n \n', \
                  'cd ../Executables \n \n']) #, \
                  
    f.close()
    g.close()
    
    Parameters = []
    
    for r in det_position:
        for d in depth:
            for sep in seperation:
                for a in alpha:
                    for b in beta:
                        for c in gamma:
                            for loc in Ubox_location:
                                
                                Parameters.append(np.array([r,d,sep,a,b,c,loc],dtype=object))
                                
                                WriteSource('DCR_{}cm_{}m_{}cm_{}a_{}b_{}c_{}m.cc'.format(*Parameters[-1]).replace(' ',''),
                                            sep, 
                                            d,
                                            r[0],
                                            r[1],
                                            True,
                                            a,b,c,loc)
                                
                                WriteSource('DCS_{}cm_{}m_{}cm_{}a_{}b_{}c_{}m.cc'.format(*Parameters[-1]).replace(' ',''),
                                            sep, 
                                            d,
                                            r[0],
                                            r[1],
                                            False,
                                            a,b,c,loc)
                                
                                WriteCommands('DCR_{}cm_{}m_{}cm_{}a_{}b_{}c_{}m.cc'.format(*Parameters[-1]).replace(' ',''),
                                              'MuonsR_{}cm_{}m_{}cm_{}a_{}b_{}c_{}m'.format(*Parameters[-1]).replace(' ',''),
                                              'Source', 
                                              'RDR_{}cm_{}m_{}cm_{}a_{}b_{}c_{}m.out'.format(*Parameters[-1]).replace(' ',''),
                                              current_date_time)
                                              
                                
                                WriteCommands('DCS_{}cm_{}m_{}cm_{}a_{}b_{}c_{}m.cc'.format(*Parameters[-1]).replace(' ',''),
                                              'MuonsS_{}cm_{}m_{}cm_{}a_{}b_{}c_{}m'.format(*Parameters[-1]).replace(' ',''),
                                              'Source', 
                                              'RDS_{}cm_{}m_{}cm_{}a_{}b_{}c_{}m.out'.format(*Parameters[-1]).replace(' ',''),
                                              current_date_time)
       
    f = open('Compile','a')
    g = open('Run','a')
    
    if chain:
        f.writelines(['cd ../ \n', \
                      'cd ../Source \n',
                      'bash Run']) #Chains together the command scripts, useful for batch jobs
    else:
        f.writelines(['cd ../ \n', \
                      'cd ../Source \n'])
        
    g.writelines(['mv Muons* ./{} \n'.format(current_date_time), \
                  'cd ../RowData \n \n', \
                  'python ClusterAnalysis.py \n \n', \
                  'mv ./{}/RD*.joblib ../Training/{} \n \n'.format(current_date_time,current_date_time), \
                  'rm ../RowData/Parameters.joblib ../Executables/Events.joblib \n \n', \
                  'rm ../Executables/inputCurrent ../Executables/r2.in \n \n', \
                  'mv ../Executables/Cts_dx_dy.out ../Executables/{} \n \n'.format(current_date_time), \
                  'rm -rf ../build \n \n', \
                  'cp time.joblib ../RowData/{} \n \n'.format(current_date_time), \
                  'mv time.joblib ../Training/{} \n \n'.format(current_date_time)])
                  #'python Plot.py'])
    
    f.close()    
    g.close()
    
    dump(Parameters,'Parameters.joblib')
    
    
    

CreateSource(Det_Position, Depth, Seperation)
WriteHeader(Z)


