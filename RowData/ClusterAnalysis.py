#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 18 13:29:30 2020

@author: keegan
"""



# import ROOT
import numpy as np
from joblib import dump, load
from time import time
import math

import glob
import os



print("Hello Viewer!")

begin_time = np.round(time(),1)



#---------------------------------- Global Variables -----------------------------------------


#---------------------
# Detector Information
#---------------------
TriggerSize = 41.6                                              #[cm]
BarWidth = 3.2                                                  #[cm]
BarHight = 1.7                                                  #[cm]
NumOfBars = 27
PlanesSeperation = 25                                           #[cm]
TriggerWidth = 1                                                #[cm]

#---------------------
# Volume Specs
#---------------------
ImageLayerSize = [1500, 1500]                                   #[cm] Dimensions of the image layers


#---------------------
# Imaging 
#---------------------
ProjectionPixel = [200, 200, 5]                                #Resolution of images, ProjectionPixel[2] is number of image layers
ClusterLayer = 16
ImageLayers = [500]                                            #[cm] height above bottom of detector for image slices

#---------------------
# Clustering 
#---------------------
Divide = [3, 3]                                               #[X divisions, Y divisions] of image layers in LocalMaxIndices() and LayerCluster()
SigFrac = 0.2                                                   #Percent of low fourier frequencies kept 
                                                                

#---------------------
# Run Options
#---------------------          
Iterate = False                                               #True ==> analysis will be done on every layer given by ProjectionPixel[2]


#---------------------
# Readout Error
#---------------------     
STD = 0 #[MeV]                                             # Standard deviation of added gaussian

#%%

# --------------------------------- Translation of Gilad's Code -----------------------------------------



def ReadRowDataFileFastest(FileName,DetectorPos): #must use 'path/Filename' (including quoatations) in function call 
    
    #--------------------------------------------------------------------
    # Reads RowData.out files and produces a Dictionary with the 
    # information contained in the file
    #--------------------------------------------------------------------
    
    begin_time = time()
    
    G = np.loadtxt(FileName, dtype = str, delimiter = '~') #Delimiter '~' is chosen in order to avoid default whitespace delimiter
    Gsize = len(G)
    
    EventWordInd = []
    n = 0
    
    for i in range(Gsize):    
        if G[i] == '*Event*':
            n += 1
            EventWordInd.append(i)
    
    RowData = {'FileName':FileName,
               'NumberOfEvents':n,
               'DateAndTime':[],
               'DetectorPos':DetectorPos, #Coordinates (in x,y plane) of the Detector
               'BarsReadout': [[],[]],
               'UpperTrigPos':[[],[]],
               'LowerTrigPos':[[],[]]} #Trig Positions are left empty, can be filled if needed
    
    for i in range(n):
        RowData['DateAndTime'].append(G[EventWordInd[i]+2])   
       
    #Error = np.random.normal(scale=STD, size=(n,4))
    #err = np.random.normal(scale)
    
    for i in range(n-1): # -1 to avoid error from length of last event 
        Bar = []
        Length = []
        
        for l in range(12):
            if 12 <= len(G[EventWordInd[i]+l]) <= 13: #Bounds depend delicately on precision of path lengths
                Bar.append(G[EventWordInd[i]+l][0:3])
                Length.append(G[EventWordInd[i]+l][4:])
    
        BarFloat = np.float_(Bar).tolist() #Converts string elements to float
        LengthFloat = (np.float_(np.array(Length))+np.random.normal(scale=STD,size=len(Length))).tolist()
        #LengthFloat = (np.float_(np.array(Length))).tolist()
        
        RowData['BarsReadout'][0].append(BarFloat)
        RowData['BarsReadout'][1].append(LengthFloat)
        
    print('There were ', n,' events in ', FileName,', the simulation ran between ', \
          RowData['DateAndTime'][0],' - ',RowData['DateAndTime'][n-1],'.')
    
    tictoc = np.round(time() - begin_time,1)
    print('It took ', tictoc,' s to read the file.') #',FileName)

    return RowData



def CalcLocalPos(Bar,Length): #Arguments are BarsReadout elements from one event
    # LocalPos == [LocalX or LocalY, LocalZ]
    
    #--------------------------------------------------------------------
    # Determines position of the muon through each layer of scintillators
    # as well as outputs list of which bars produced the signal for each
    # event
    #--------------------------------------------------------------------
    
    LengthLocal = [[],[],[],[]]
    BarLocal = [[],[],[],[]]
    LocalPos = [[],[],[],[]]
    
    for n in range(len(Bar)): #Sorts Bar and Length data into nested lists corresponding to each layer on the detector
         LengthLocal[math.floor(Bar[n]/100)-1].append(Length[n])
         BarLocal[math.floor(Bar[n]/100)-1].append(Bar[n])
         
    a = np.sqrt(BarHight**2+(BarWidth/2)**2)
    Alpha = np.arctan(2*BarHight/BarWidth)
    
    for i in range(len(LengthLocal)):
        NumOfBarsLocal = len(LengthLocal[i])
        
        if NumOfBarsLocal == 0: # No signal, return error
            X,Z = -9999,-9999
        
        elif NumOfBarsLocal == 1: 
            X,Z = BarWidth/2, BarHight/2
            '''
            if BarLocal[i][0]%2 == 0: #The first bar's vertex is facing down
                X,Z = 0, 0
            
            else:
                X,Z = 0, BarHight
            '''
            #Takes tip instead of middle of bar
            
            
        elif NumOfBarsLocal >= 2:
            Readout = []
            mxind = LengthLocal[i].index(max(LengthLocal[i]))
            
            if mxind == 0: #The first bar has the max readout so we take the first and second bars
                Readout = [LengthLocal[i][mxind],LengthLocal[i][mxind+1]]
            
            elif mxind == len(LengthLocal[i])-1: #The last bar has the max readout so we take the last bar and the one before
                Readout = [LengthLocal[i][mxind],LengthLocal[i][mxind-1]]
            
            else: #The max readout is somewhere at the middle, so we take this bar and it's highest neighbor
                Readout = [LengthLocal[i][mxind], np.amax([LengthLocal[i][mxind+1], LengthLocal[i][mxind-1]])]

            if BarLocal[i][0]%2 == 0: #The first bar's vertex is facing down
                X = BarWidth/2 - (a*Readout[0]/(Readout[0]+Readout[1]))*math.cos(Alpha)
                Z = BarHight/2 - (a*Readout[0]/(Readout[0]+Readout[1]))*math.sin(Alpha)
            
            else: #The first bar's vertex is facing up
                X = -BarWidth/2 + (a*Readout[1]/(Readout[0]+Readout[1]))*math.cos(Alpha)
                Z = BarHight/2 - (a*Readout[1]/(Readout[0]+Readout[1]))*math.sin(Alpha)
            
        LocalPos[i] = [X,Z] #[X/2 + np.random.random()*X, Z/2 + np.random.random()*Z] #randomize 50% (increases computing time)
            
    return LocalPos, BarLocal 



def CalcAbsPos(LocalPos,BarLocal, Seperation):
    #AbsPos[i] == [AbsX or AbsY, AbsZ] 
    
    #--------------------------------------------------------------------
    # Deterimines position relative to centre of the detector layer, 
    # and height measured from the bottom of the detector
    #--------------------------------------------------------------------
    
    AbsPos = [[0,0],[0,0],[0,0],[0,0]]
    
    if [-9999,-9999] in LocalPos:
        return -9999
        
    else:    
        for l in range(len(LocalPos)):  
            FirstBarIndex = BarLocal[l][0]
            i = math.floor(FirstBarIndex/100)
            
            if i == 1: #XUp
                AbsPos[l][1] = LocalPos[l][1] + TriggerWidth + Seperation + 3.5 * BarHight
                FirstBarIndex = FirstBarIndex - 100
           
            if i == 2: #YUp
                AbsPos[l][1] = LocalPos[l][1] + TriggerWidth + Seperation + 2.5 * BarHight
                FirstBarIndex = FirstBarIndex - 200
            
            if i == 3: #XDown
                AbsPos[l][1] = LocalPos[l][1] + TriggerWidth + 1.5 * BarHight
                FirstBarIndex = FirstBarIndex - 300
            
            if i == 4: #YDown
                AbsPos[l][1] = LocalPos[l][1] + TriggerWidth + 0.5 * BarHight
                FirstBarIndex = FirstBarIndex - 400
            
            
            if FirstBarIndex%2 == 0 : #The first bar's vertex is facing down
                AbsPos[l][0] = LocalPos[l][0] - (NumOfBars / 4 - 0.25) * BarWidth + BarWidth / 2 * FirstBarIndex
                
            else: #The first bar's vertex is facing up
                AbsPos[l][0] = LocalPos[l][0] - (NumOfBars / 4 - 0.25) * BarWidth + BarWidth / 2 * (FirstBarIndex + 1) 
                
        return AbsPos 



def CalcEventHittingPoints(Bar, Length, ZImage, DetectorPos, Seperation): #( <RowData>['BarsReadout'][0][i], <RowData>['BarsReadout'][1][i], ...)
    
    #--------------------------------------------------------------------
    # Determines where on each plane (z = const.) the muon passed through 
    #--------------------------------------------------------------------
    
    HittingPoints = np.zeros(9).reshape((3,3))
    
    [LocalPos,BarLocal] = CalcLocalPos(Bar,Length)
    
    AbsPos = CalcAbsPos(LocalPos,BarLocal, Seperation)
    
    if AbsPos == -9999:
        HittingPoints = [[-9999]*3]*3
    
    else:    
        AbsXUp = AbsPos[0]
        AbsYUp = AbsPos[1]
        AbsXDown = AbsPos[2]
        AbsYDown = AbsPos[3]
        
        ZUp = 2 * TriggerWidth + 4 * BarHight + Seperation
        ZSurf = DetectorPos[2] + ZUp
        
        dZx = AbsXUp[1] - AbsXDown[1]           
        dX = AbsXUp[0] - AbsXDown[0]            
        
        dZy = AbsYUp[1] - AbsYDown[1]
        dY = AbsYUp[0] - AbsYDown[0]
        
        if dX != 0 and dY != 0:
            ax = dZx / dX  
            ay = dZy / dY
        
            bx = AbsXUp[1] - ax * AbsXUp[0]
            by = AbsYUp[1] - ay * AbsYUp[0]
            
            HittingPoints[0][0] = (ZImage - bx) / ax + DetectorPos[0] #XImage
            HittingPoints[0][1] = (ZImage - by) / ay + DetectorPos[1] #YImage
            HittingPoints[0][2] = ZImage #ZImage
            
            HittingPoints[1][0] = (ZUp - bx) / ax + DetectorPos[0] #XUp
            HittingPoints[1][1] = (ZUp - by) / ay + DetectorPos[1] #YUp
            HittingPoints[1][2] = ZUp #ZUp
            
            HittingPoints[2][0] = (ZSurf - bx) / ax + DetectorPos[0] #XSurf
            HittingPoints[2][1] = (ZSurf - by) / ay + DetectorPos[1] #YSurf
            HittingPoints[2][2] = ZSurf #ZSurf
            
        else:
            HittingPoints = [[-9999]*3]*3
    
    return HittingPoints 



def PazAnalysis(RowData, Seperation, Iterate, TopDepth, ImgLayers):
    
    #--------------------------------------------------------------------
    # Creates a 3D array counting the number of trajectories passsing 
    # through each pixel in the image layers specified by the global
    # variables 
    #--------------------------------------------------------------------
    
    begin_time = time()
    
    N = RowData['NumberOfEvents']
    
    if Iterate == True:
        DetectorCounts = np.zeros((ProjectionPixel[0], ProjectionPixel[1], ProjectionPixel[2])) 
        
        for k in range(N-1):
            for i in range(ProjectionPixel[2]):
                ZImage = 2 * TriggerWidth + 4 * BarHight + Seperation + i*TopDepth/ProjectionPixel[2] #Z coordinate of Image Layer
                HittingPoints = CalcEventHittingPoints(RowData['BarsReadout'][0][k], RowData['BarsReadout'][1][k], ZImage, RowData['DetectorPos'], Seperation)
                
                if np.all(HittingPoints == -9999):
                    DetectorCounts = -9999
                    break
                
                else:  
                    Iind = int(round((HittingPoints[0][0] + ImageLayerSize[0]/2)*(ProjectionPixel[0]-1)/ImageLayerSize[0])) 
                    Jind = int(round((HittingPoints[0][1] + ImageLayerSize[1]/2)*(ProjectionPixel[1]-1)/ImageLayerSize[1]))
                    
                    if Iind>0 and Jind>0 and Iind<len(DetectorCounts[0]) and Jind<len(DetectorCounts[1]): 
                            DetectorCounts[Iind,Jind,i] += 1 
        
    else:
        DetectorCounts = np.zeros((ProjectionPixel[0], ProjectionPixel[1],len(ImageLayers))) 
        
        for k in range(N-1):
            for depth in ImgLayers:
                ZImage = depth
                #ZImage = 2 * TriggerWidth + 4 * BarHight + Seperation + ClusterLayer*TopDepth/ProjectionPixel[2] #Z coordinate of Image Layer
        
                HittingPoints = CalcEventHittingPoints(RowData['BarsReadout'][0][k], RowData['BarsReadout'][1][k], ZImage, RowData['DetectorPos'], Seperation)
                
                if np.all(HittingPoints == -9999):
                    DetectorCounts = -9999
                    break
                
                else:  
                    Iind = int(round((HittingPoints[0][0] + ImageLayerSize[0]/2)*(ProjectionPixel[0]-1)/ImageLayerSize[0])) 
                    Jind = int(round((HittingPoints[0][1] + ImageLayerSize[1]/2)*(ProjectionPixel[1]-1)/ImageLayerSize[1]))
                    
                    if Iind>0 and Jind>0 and Iind<len(DetectorCounts[0]) and Jind<len(DetectorCounts[1]): 
                            DetectorCounts[Iind,Jind] += 1 
    
    tictoc = np.round(time() - begin_time,1) 
    print('It took ', tictoc,' s to analyse the RowData from ', RowData['FileName'],'\n\n')
    
    return DetectorCounts 



#%%
#---------------------------------- Beam Manipulation Functions --------------------------------------



def ObjectView(Data, Resolution, TopDepth, Seperation):
    
    #-----------------------------------------------------------------
    # Takes Data hitting points and counts hits in the ImageVolume at 
    # ObjectZ above the detector with resolution Resolution
    #-----------------------------------------------------------------
    
    DetectorCounts = np.zeros((Resolution[0],Resolution[1],Resolution[2]))
    
    for i in range(len(Data)):
        AbsXUp = Data[i][0][2][0]
        AbsYUp = Data[i][0][2][1]  
        AbsXDown = Data[i][0][1][0]
        AbsYDown = Data[i][0][1][1]
        
        ZUp = Data[i][0][2][2]             
        ZDown = Data[i][0][1][2]
        dZ = ZUp - ZDown
        
        if dZ == 0:
            continue                    
        
        dX = AbsXUp - AbsXDown     
        ax = dZ / dX                           
        bx = ZUp - ax * AbsXUp         
        
        dY = AbsYUp - AbsYDown
        ay = dZ / dY
        by = ZUp - ay * AbsYUp
        
        for j in range(Resolution[2]):
            ZImage = 2 * TriggerWidth + 4 * BarHight + Seperation + j*TopDepth/Resolution[2]
            
            HittingPoints = np.zeros(3)
            
            HittingPoints[0] = (ZImage - bx) / ax #XImage  
            HittingPoints[1] = (ZImage - by) / ay #YImage
            HittingPoints[2] = ZImage #ZImage
            
            Iind = int(round((HittingPoints[0] + ImageLayerSize[0]/2)*(Resolution[0]-1)/ImageLayerSize[0]))
            Jind = int(round((HittingPoints[1] + ImageLayerSize[1]/2)*(Resolution[1]-1)/ImageLayerSize[1]))
        
            if Iind > 0 and Jind > 0 and Iind < len(DetectorCounts[0]) and Jind < len(DetectorCounts[1]):
                DetectorCounts[Iind][Jind][j] += 1
                    
    return DetectorCounts 


        
def GroupOverlaps(Data):
    #Expects 3D array for Data
    
    #--------------------------------------------------------------------
    # Checks to see if each pair of object data arrays overlap by seeing 
    # if subtracting one from the other changes any of the values in the 
    # first array.  
    # Then forms groups of overlapping sets of data (each of which is 
    # interpretted as an object)
    #--------------------------------------------------------------------
    
    Data = np.array(Data)       
           
    N = len(Data)
    OverlapLists = []
    
    for i in range(N):
        for j in range(N):
            if i != j: 
                Minus = (Data[i] - Data[j])
                Minus[Minus < 0] = 0 # sets negative pixels to zero
    
                if np.any(Data[i] != Minus):
                    OverlapLists.append([i,j])    
    
    ObjectGroups = []
    
    for i in range(len(OverlapLists)):
        if len(ObjectGroups) == 0:
            ObjectGroups.append([OverlapLists[i][0],OverlapLists[i][1]])
    
        Appended = 0 #Used to keep track of whether element is already in ObjectGroups across if gates
    
        for j in range(len(ObjectGroups)):
            if not OverlapLists[i][0] in ObjectGroups[j]:
                if OverlapLists[i][1] in ObjectGroups[j]:
                    ObjectGroups[j].append(OverlapLists[i][0])
                    Appended += 1                
    
            elif not OverlapLists[i][1] in ObjectGroups[j]:
                if OverlapLists[i][0] in ObjectGroups[j]:
                    ObjectGroups[j].append(OverlapLists[i][1])
                    Appended += 1
            
            if OverlapLists[i][0] in ObjectGroups[j]:
                if OverlapLists[i][1] in ObjectGroups[j]:
                    Appended += 1
                    
        if Appended == 0:
            ObjectGroups.append([OverlapLists[i][0],OverlapLists[i][1]])
            
    return ObjectGroups, OverlapLists
        


def ScaleGroups(Counts, Groups, Maxima, OverlapCutoff):
    
    #--------------------------------------------------------------------
    # Scales all beams identically so regions of interference can be            
    # identified and decides whether pairs of 3D arrays overlap to 
    # create Target. Removes the need for the classification algorithm 
    # with this step.
    #--------------------------------------------------------------------
    
    DetectorCounts = np.zeros((len(Groups),ProjectionPixel[0],ProjectionPixel[1],ProjectionPixel[2]))
    AddedMaxima = []
    Targets = np.zeros((len(Groups),ProjectionPixel[2]))
    
    for i in range(len(Groups)):
        Temp = []
        
        for j in Groups[i]: # Decides how many pairs of beams overlap at each layer by permuting them within their groups
            for k in Groups[i]:
                if j != k: 
                    for m in range(ProjectionPixel[2]):
                        Minus = (Counts[j,:,:,m] - Counts[k,:,:,m])
                        Minus[Minus < 0] = 0
            
                        #if np.any(Counts[j,:,:,m] != Minus):
                        if abs(np.sum(Counts[j,:,:,m] - Minus)) > np.max(Counts[j,:,:,m]) * OverlapCutoff: 
                            Targets[i,m] += 1
                        
                        else:
                            Targets[i,m] += 2
                        
            if Maxima[j] != 0 and not Maxima[j] in Temp: # assumes that if maxima are the same then clustered array is a duplicate
                ScaledData = Counts[j] * np.max(Maxima)/Maxima[j] 
            
                DetectorCounts[i] += ScaledData
    
                Temp.append(Maxima[j])
        
        AddedMaxima.append(Temp)
        
    return DetectorCounts, AddedMaxima, Targets



def ScaleLayers(Counts, Scale):
    
    if Scale == True:
        Shape = np.shape(Counts)
        TempCounts = np.zeros(Shape)
        '''
        Counts[Counts > 0] = 100
        
        TempCounts = Counts
        '''
        for i in range(Shape[0]):
            Max = np.max(Counts[i])
            
            for j in range(ProjectionPixel[2]):
                if np.max(Counts[i,:,:,j]) != 0:
                    TempCounts[i,:,:,j] = Counts[i,:,:,j] * Max / np.max(Counts[i,:,:,j])
        
    else:
        TempCounts = Counts
    
    return TempCounts
    
#%%
#---------------------------------- Study Resolution --------------------------------------



def ResolveError(ClArray):
    '''
    Finds standard deviation in clustered image for standard object
    compares nonzero pixels in expected and clustered image based on
    depth
    
    *** neglects small differences in position ~ 0.1m
    
    input:
        ClArray [float] - [2D] clustered array
        
    output:
        Sigma [float] - Percent of incorrect pixels 
    '''
    Shape = np.array(np.shape(ClArray))
    
    Inds = np.array(Shape * 100 / np.array(ImageLayerSize) / 2, dtype=int) #number of indices occupied
                                                     # by half the object
    ExpArray = np.zeros(Shape)
    
    ExpArray[int(Shape[0]/2)-Inds[0]:int(Shape[0]/2)+Inds[0],int(Shape[1]/2)-Inds[1]:int(Shape[1]/2)+Inds[1]] = 1
    
    ClArray[ClArray > 0] = 1
    
    Sigma = np.count_nonzero(ClArray - ExpArray)
    SigmaPer = 100 * np.count_nonzero(ClArray - ExpArray) / (Shape[0] * Shape[1])
    
    ### How should size of object be affected by height?
    
    return np.array([SigmaPer,Sigma])
    

#%%
#---------------------------------- Get Data --------------------------------------



def ReadDataFiles():
    
    Parameters = load('Parameters.joblib')
    Dir_name = load('time.joblib')
    
#    RealFiles = [Dir_name+'/'+'RDR_{}cm_{}m_{}cm_{}a_{}b_{}c_{}m.out'.format(*par).replace(' ','') for par in Parameters]
#    SkyFiles = [Dir_name+'/'+'RDS_{}cm_{}m_{}cm_{}a_{}b_{}c_{}m.out'.format(*par).replace(' ','') for par in Parameters] #int(par[0]),par[1],int(par[2])
#    ZPositions = [par[1] for par in Parameters] #[m]
#    Seperations = [par[2] for par in Parameters]
#    Shape = np.shape(Parameters)
#
    global names
    names = [filename[filename.find('/')+4:] for filename in glob.iglob('{}/RDR*.out'.format(Dir_name), recursive=True)]
    
    RealFiles = [Dir_name+'/'+'RDR'+name for name in names]
    SkyFiles = [Dir_name+'/'+'RDS'+name for name in names]
    ZPositions = [float(names[0].split('m_')[1]) for i in range(len(names))]
    Seperations = [float(names[i].split('cm')[1][-2:]) for i in range(len(names))]
    Shape = [len(names)]
        
    
    XPositions = np.zeros(Shape[0]) # centres images
    #XPositions = [round(np.sin(par[1])*par[0]*100) for par in Parameters] #[cm]
    YPositions = np.zeros(Shape[0])
    
    SkyCountList = []
    RealCountList = []
    CountList = []
    RowSkyList = []
    RowRealList = []
    
    #sign = input('Input True if you expect a dense object \n \n')
    sign = False
    #sign = 'True'
    
    for i in range(len(RealFiles)):
        RDSky = ReadRowDataFileFastest(RealFiles[i], [XPositions[i],YPositions[i],ZPositions[i]])
        RDReal = ReadRowDataFileFastest(SkyFiles[i], [XPositions[i],YPositions[i],ZPositions[i]])
    
        DCS = PazAnalysis(RDSky,Seperations[i],Iterate,ZPositions[i],ImageLayers) 
        DCR = PazAnalysis(RDReal,Seperations[i],Iterate,ZPositions[i],ImageLayers)
        
        if sign == 'True':
            DCdatPlus = (DCS-DCR) 
            
        else:
            DCdatPlus = (DCR-DCS) 
        
        DCdatPlus[DCdatPlus < 0] = 0
        
        dump(DCS,SkyFiles[i][:SkyFiles[i].find('.out')]+'.joblib')
        dump(DCR,RealFiles[i][:RealFiles[i].find('.out')]+'.joblib')
        
        RowSkyList.append(RDSky)
        RowRealList.append(RDReal)
        SkyCountList.append(DCS) 
        RealCountList.append(DCR)
        
        CountList.append(DCdatPlus)    
            
    ReadDict = {'Length':Shape[0], \
                'Iterate':Iterate, \
                'Parameters':Parameters, \
                #'Row Sky List':RowSkyList, \
                #'Row Real List':RowRealList, \
                'Sky Count List':SkyCountList, \
                'Real Count List':RealCountList, \
                'Subtracted Count List':CountList, \
                'Seperations':Seperations}
    
    RowDict = {'Row Sky List':RowSkyList, \
                'Row Real List':RowRealList}
    
    return ReadDict, RowDict



def ReadFile(RealFile, SkyFile, Position, Seperation):
    
#        RDSky = ReadRowDataFileFastest(RealFile, Position)
#        DCS = PazAnalysis(RDSky,Seperation,Iterate,Position[2], ImageLayers)
#        dump(DCS,SkyFile[:SkyFile.find('.out')]+'.joblib')
#
        RDReal = ReadRowDataFileFastest(SkyFile, Position)
        DCR = PazAnalysis(RDReal,Seperation,Iterate,Position[2],ImageLayers)
        dump(DCR,RealFile[:RealFile.find('.out')]+'.joblib')

#%%           
#---------------------------------- Call Functions --------------------------------------



#-----------------------------------------------------------------
# Runs the imaging procedure on data specified by the user on
# the command line. Creates 2D images as well as tracked beams
# and estimates of the object sizes and shapes visually
#-----------------------------------------------------------------

if __name__ == '__main__':
    ReadDictionary, RowDictionary = ReadDataFiles()

#    ReadFile('Real','Sky',[0,0,3000], 25)

    #dump(ReadDictionary,'ReadDictionary.joblib')
    #dump(RowDictionary,'RowDictionary.joblib')
    #ReadDictionary = load('ReadDictionary.joblib')

    #ClusterDictionary  = Cluster(ReadDictionary)

    #ClustError = [ResolveError(ClusterDictionary['Cluster'][i]['Filled Array']) for i in range(ReadDictionary['Length'])]
    #FFTError = [ResolveError(ClusterDictionary['FFT Cluster'][i]['Filled Array']) for i in range(ReadDictionary['Length'])]

    #dump(ReadDictionary['Sky Count List'],'Sky_Arrays.joblib')
    #dump(ReadDictionary['Real Count List'],'Real_Arrays.joblib')
    #dump(ClusterDictionary,'ClusterDict.joblib')
    #dump(SubCounts,'SubCounts.joblib')
    #dump(np.array([ClustError,FFTError]),'Error.joblib')

    ###change to only run on cluster layer specified by objdepth etc..

#--------------------------------------------------------------------------------------------



tictoc = time() - begin_time

print('The total runtime of the analysis was: {:.02g} s'.format(tictoc))

#print('The total runtime of the simulation was: {} s'.format(load('time.joblib')))
