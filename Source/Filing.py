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
beta = np.round([180 * mult * np.arccos(0.25) / np.pi for mult in [-1,1]], 1) #[deg]
#beta = [0] #[deg]
gamma = [0] #[deg]

Radius = [20] #[m] ~ distance from origin

#Theta = np.round(np.pi / 2 * np.array([0, 1, 2, 3, 4]) / 10, 2) #[rad] polar angle
#Theta = np.round([0], 2)
Theta = np.round(np.pi * beta / 180,1) #[rad]

Seperation = [30] #[cm]
#Seperation = [40] #[cm] seperation of the detector planes

Events = int(5e5)
#Events = int(5e6)
dump(Events, 'Events.joblib')

Ubox_Material = ['BMF', 'Air', 'GroundMat', 'Vacuum'][0]
Ubox_dims = [1, 100, 0.5] # [m] (x, y, z) half dimensions
Ubox_location = [[5,0,1],[-5,0,1],[10,0,1],[-10,0,1]] #[m] (x,y,z)

Inhomogeneous = True # whether or not landscape is added
Landscape_material = ['BMF', 'GroundMat'][1]



#%%
# Relief Data

#dim = 60
#
#x = np.arange(1,dim+1).reshape(dim,1)
#y = np.arange(1,dim+1).reshape(1,dim)
#X,Y = np.meshgrid(x,y)
#
#Z = 25 * np.exp(-((X-25)**2 + (Y-25)**2)/5) + 25 * np.exp(-(((X-35)**2 + (Y-35)**2)/5))
#
#Z = np.zeros((dim,dim))
#Z[30:35,30:35] = 5 * np.ones((5,5))



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
    
    
    g_lines[217] = '\tG4double UBoxHight = %s*m; \n'%(Ubox_dims[2])
    g_lines[218] = '\tG4VSolid* UBox = new G4Box("UBox", %s*m, %s*m, UBoxHight); \n'%(Ubox_dims[0], Ubox_dims[1])
    
    if Inhomogeneous:
        LandscapeList = ['\t \t /* Landscape Volume */ \n ', \
                        '\t G4VSolid* solidLandscape = new G4Box("Landscape", UD_World_sizeX, UD_World_sizeY, (UD_World_sizeZ - DetectorDepth) / 2.0); \n ', \
                        '\t G4LogicalVolume* logicLandscape = new G4LogicalVolume(solidLandscape, Vacuum, "Landscape"); \n ', \
                        '\t  \n ', \
                        '\t new G4PVPlacement(0,               // no rotation \n ', \
                        '\t     G4ThreeVector(0.0, 0.0, DetectorDepth + (UD_World_sizeZ - DetectorDepth) / 2.0), //above Homogeneous ground \n ', \
                        '\t     logicLandscape, \n ', \
                        '\t     "Landscape",         // its name \n ', \
                        '\t     logicWorld,  // its mother  volume \n ', \
                        '\t     false,           // no boolean operations \n ', \
                        '\t     0); \n ', \
                        '\t  \n ', \
                        '\t  \n ', \
                        '\t G4double RBox_1_height = 1 * m; \n ', \
                        '\t G4VSolid* RBox_1 = new G4Box("RBox_1", Pnt_Sep/2.0, Pnt_Sep/2.0, RBox_1_height/2.0); \n ', \
                        '\t G4LogicalVolume* RBox_1_LV = new G4LogicalVolume(RBox_1, %s, "RBox_1"); \n '%(Landscape_material), \
                        '\t  \n ', \
                        '\t int RBox_copy = 0; \n ', \
                        '\t  \n ', \
                        '\t for (int i = 0; i < iaxis - 1; i++) \n ', \
                        '\t { \n ', \
                        '\t     for (int j = 0; j < jaxis - 1; j++) \n ', \
                        '\t     { \n ', \
                        '\t         if (RMins[i][j] >= 1.0) \n ', \
                        '\t         { \n ', \
                        '\t             for (int k = 0; k < RMins[i][j]; k++) \n ', \
                        '\t             { \n ', \
                        '\t                 G4RotationMatrix* rotBox_null = new G4RotationMatrix(); \n', \
                        '\t                 rotBox_null->rotateX(0.0*deg); \n', \
                        '\t                 G4ThreeVector RBox_Centre = G4ThreeVector((1-iaxis) * Pnt_Sep /2.0 + Pnt_Sep * i, (1-jaxis) * Pnt_Sep /2.0 + Pnt_Sep * j, (0.5 + k) * m - (UD_World_sizeZ - DetectorDepth) / 2.0); \n ', \
                        '\t                 new G4PVPlacement(rotBox_null,               // no rotation \n ', \
                        '\t                         RBox_Centre, \n', \
                        '\t                         RBox_1_LV, \n ', \
                        '\t                         "RBox_1",         // its name \n ', \
                        '\t                         logicLandscape,  // its mother  volume \n ', \
                        '\t                         false,           // no boolean operations \n ', \
                        '\t                         RBox_copy); \n ', \
                        '\t  \n ', \
                        '\t                 RBox_copy++; \n ', \
                        '\t            } \n ', \
                        '\t            G4cout << G4endl << RMins[i][j] * m - (UD_World_sizeZ - DetectorDepth) / 2.0 << " mm is local height " << (1-iaxis) * Pnt_Sep /2.0 + Pnt_Sep * i << ", " << (1-jaxis) * Pnt_Sep /2.0 + Pnt_Sep * j << " mm are the x,y positions " << G4endl; \n ', \
                        '\t  \n ', \
                        '\t         } \n ', \
                        '\t     } \n ', \
                        '\t } \n ', \
                        '\t  \n ', \
                        '\t G4cout << G4endl << Pnt_Sep/2.0 << "  " << Pnt_Sep/2.0 << "  " << RBox_1_height /2.0 << G4endl << G4endl; \n ', \
                        '\t //G4cout << G4endl << RBox_copy << " boxes were placed to create the Landscape" << G4endl << G4endl; \n ', \
                        '\t G4cout << G4endl << (UD_World_sizeZ - DetectorDepth) / 2.0 << " mm is the half height of the sky" << G4endl << G4endl; \n ' ]
        
        LandscapeList.reverse()
        [g_lines.insert(381,Str) for Str in LandscapeList]
    
    if Real:
        g_lines[219] = '\tG4LogicalVolume* UBoxLogicVolume = new G4LogicalVolume(UBox, %s, "UBox");\n'%(Ubox_Material)
        g_lines[220] = '\t//G4LogicalVolume* UBoxLogicVolume = new G4LogicalVolume(UBox, GroundMat, "UBox");\n'
    
    else:
        g_lines[219] = '\t//G4LogicalVolume* UBoxLogicVolume = new G4LogicalVolume(UBox, %s, "UBox");\n'%(Ubox_Material)
        g_lines[220] = '\tG4LogicalVolume* UBoxLogicVolume = new G4LogicalVolume(UBox, GroundMat, "UBox");\n'
    
    g_lines[224] = '\tG4double UboxDepth = {} * m;'.format(boxloc[2])
    g_lines[227] = '\tG4ThreeVector({} * m, {} * m, DetectorDepth / 2 - UBoxHight / 2 - UboxDepth),'.format(boxloc[0],boxloc[1])
    
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
    g_lines.append('G4double max_landscape_height = {}; \n \n'.format(np.amax(Z)))
    
    #paste data
    g_lines.append('\n std::array<std::array<G4double, %s>, %s> Relief = {{ \n'%(Shape[1],Shape[0]))
    
    for i in range(Shape[0]-1):
        Str = '{{'

        for j in range(Shape[1]-1):
            Str += str(Data[i,j]) + ', '

        Str += str(Data[i,Shape[1]-1]) + '}}, \n'

        g_lines.append(Str)

    Str = '{{'

    for j in range(Shape[1]-1):
        Str += str(Data[i,j]) + ', '

    Str += str(Data[i,Shape[1]-1]) + '}} \n'

    g_lines.append(Str)

    g_lines.append('}}; \n \n')

    #paste minima
    g_lines.append('\n std::array<std::array<G4double, %s>, %s> RMins = {{ \n'%(Shape[1]-1,Shape[0]-1))

    for i in range(Shape[0]-2):
        Str = '{{'

        for j in range(Shape[1]-2):
            Str += str(Mins[i,j]) + ', '

        Str += str(Mins[i,Shape[1]-2]) + '}}, \n'

        g_lines.append(Str)

    Str = '{{'

    for j in range(Shape[1]-2):
        Str += str(Mins[i,j]) + ', '
	
    Str += str(Mins[i,Shape[1]-2]) + '}} \n'

    g_lines.append(Str)

    g_lines.append('}}; \n \n')
    #g_lines.append('} \n \n')

    g.writelines(g_lines)

    g.close()

    
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
                  'The distance from the origin is: \t {} m \n\n'.format(Radius),
                  'The angle of the detector relative to the vertical is: \t {} rad \n\n'.format(*Theta),
                  'The seperation of the detector planes is: \t {} cm \n\n'.format(Seperation),
                  'The euler angles of the detector are: \t {}, {}, {} deg \n\n'.format(alpha,beta,gamma),
                  'The number of events each run is: \t {} \n\n'.format(Events),
                  'The material of the box is: \t {} \n\n'.format(Ubox_Material),
                  'The dimensions of the box are: \t ({}, {}, {}) m \n\n'.format(*Ubox_dims),
                  'The depth of the box is: \t {} m \n\n'.format(Ubox_location),
                  'The inhomogeneous landscape {} simulated \n\n'.format(["wasn't",'was'][Inhomogeneous]),
                  'The material of the landscape {} {} \n\n'.format(['would have been','is'][Inhomogeneous],Landscape_material)])
        
    f.close()
    
    

def CreateSource(radius, theta, seperation):
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
#    r = open('r2.in','w+') #beware!!! using the same seed for real and sky
#    
#    r.writelines(['/run/initialize \n \n', \
#                  '/SetSeed/set {} \n \n'.format(''.join(random.choices(string.digits, k=10))), \
#                  '/gun/particle mu- \n \n', \
#                  '#/vis/scene/endOfEventAction accumulate \n \n', \
#                  '/run/beamOn {} \n \n'.format(Events)])
#    
#    r.close()
    
    
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
                  'mv DESCRIPTION.txt ../Training/{} \n \n'.format(current_date_time), \
                  'mv DC*.cc ./{} \n'.format(current_date_time), \
                  'cd ./{} \n'.format(current_date_time)])
                  
    g.writelines(['#!/bin/bash \n \n', \
                  'mv Parameters.joblib ../RowData \n \n', \
                  'mv Events.joblib ../Executables \n \n', \
                  'mv time.joblib ../RowData \n \n', \
                  'cd ../Executables \n \n']) #, \
                  #"alias python='/tomerv/SENSEI1/Keegan/Python-3.9.0/python' \n \n"])
    
    f.close()
    g.close()
    
    Parameters = []
    
    for r in radius:
        for t in theta:
            for sep in seperation:
                for a in alpha:
                    for b in beta:
                        for c in gamma:
                            for loc in Ubox_location:
                                
                                Parameters.append(np.array([r,t,sep,a,b,c,loc],dtype=object))
                                
                                WriteSource('DCR_{}m_{}rad_{}cm_{}a_{}b_{}c_{}m.cc'.format(*Parameters[-1]).replace(' ',''), 
                                            sep, 
                                            np.round(np.cos(t)*r,2),
                                            np.round(np.sin(t)*r*100),
                                            0,
                                            True,
                                            a,b,c,loc)
                                
                                WriteSource('DCS_{}m_{}rad_{}cm_{}a_{}b_{}c_{}m.cc'.format(*Parameters[-1]).replace(' ',''),
                                            sep, 
                                            np.round(np.cos(t)*r,2),
                                            np.round(np.sin(t)*r*100),
                                            0,
                                            False,
                                            a,b,c,loc)
                                
                                WriteCommands('DCR_{}m_{}rad_{}cm_{}a_{}b_{}c_{}m.cc'.format(*Parameters[-1]).replace(' ',''),
                                              'MuonsR_{}m_{}rad_{}cm_{}a_{}b_{}c_{}m'.format(*Parameters[-1]).replace(' ',''),
                                              'Source', 
                                              'RDR_{}m_{}rad_{}cm_{}a_{}b_{}c_{}m.out'.format(*Parameters[-1]).replace(' ',''),
                                              current_date_time)
                                              
                                
                                WriteCommands('DCS_{}m_{}rad_{}cm_{}a_{}b_{}c_{}m.cc'.format(*Parameters[-1]).replace(' ',''), 
                                              'MuonsS_{}m_{}rad_{}cm_{}a_{}b_{}c_{}m'.format(*Parameters[-1]).replace(' ',''), 
                                              'Source', 
                                              'RDS_{}m_{}rad_{}cm_{}a_{}b_{}c_{}m.out'.format(*Parameters[-1]).replace(' ',''),
                                              current_date_time)
                                
                                #r,t,sep,a,b,c,dep),
       
    f = open('Compile','a')
    g = open('Run','a')
    
    f.writelines(['cd ../ \n', \
                  #'cd ../Source \n',
                  'bash Run']) #Chains together the command scripts, useful for batch jobs
        
    g.writelines(['mv Muons* ./{} \n'.format(current_date_time), \
                  'cd ../RowData \n \n', \
                  'python ClusterAnalysis.py \n \n', \
                  'mv ./{}/RD*.joblib ../Training/{} \n \n'.format(current_date_time,current_date_time)])
                  #'python Plot.py'])
    
    f.close()    
    g.close()
    
    dump(Parameters,'Parameters.joblib')
    
    
    

CreateSource(Radius, Theta, Seperation)
WriteHeader(Z)


