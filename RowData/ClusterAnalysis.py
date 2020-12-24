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
ImageLayerSize = [1000, 1000]                                   #[cm] Dimensions of the image layers 


#---------------------
# Imaging 
#---------------------
ProjectionPixel = [200, 200, 5]                                #Resolution of images, ProjectionPixel[2] is number of image layers
ClusterLayer = 16

#---------------------
# Clustering 
#---------------------
Divide = [3, 3]                                               #[X divisions, Y divisions] of image layers in LocalMaxIndices() and LayerCluster()
SigFrac = 0.2                                                   #Percent of low fourier frequencies kept 
                                                                

#---------------------
# Run Options
#---------------------          
Iterate = True                                                #True ==> analysis will be done on every layer given by ProjectionPixel[2]


# --------------------------------- Translation of Gilad's Code -----------------------------------------



def ReadRowDataFileFastest(FileName,DetectorPos): #must use 'path/Filename' (including quoatations) in function call 
    
    #--------------------------------------------------------------------
    # Reads RowData.out files and produces a Dictionary with the 
    # information contained in the file
    #--------------------------------------------------------------------
    
    begin_time = np.round(time(),1)
    
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
       
    for i in range(n-1): # -1 to avoid error from length of last event 
        Bar = []
        Length = []
        
        for l in range(12):
            if 12 <= len(G[EventWordInd[i]+l]) <= 13: #Bounds depend delicately on precision of path lengths
                Bar.append(G[EventWordInd[i]+l][0:3])
                Length.append(G[EventWordInd[i]+l][4:])
    
        BarFloat = np.float_(Bar).tolist() #Converts string elements to float
        LengthFloat = np.float_(Length).tolist()
        
        RowData['BarsReadout'][0].append(BarFloat)
        RowData['BarsReadout'][1].append(LengthFloat)
        
    print('There were ', n,' events in ', FileName,', the simulation ran between ', \
          RowData['DateAndTime'][0],' - ',RowData['DateAndTime'][n-1],'.')
    
    tictoc = np.round(time(),1) - begin_time
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



def PazAnalysis(RowData, Seperation, Iterate, TopDepth):
    
    #--------------------------------------------------------------------
    # Creates a 3D array counting the number of trajectories passsing 
    # through each pixel in the image layers specified by the global
    # variables 
    #--------------------------------------------------------------------
    
    begin_time = np.round(time(),1)
    
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
        DetectorCounts = np.zeros((ProjectionPixel[0], ProjectionPixel[1])) 
        
        for k in range(N-1):
            ZImage = 2 * TriggerWidth + 4 * BarHight + Seperation + ClusterLayer*TopDepth/ProjectionPixel[2] #Z coordinate of Image Layer
            HittingPoints = CalcEventHittingPoints(RowData['BarsReadout'][0][k], RowData['BarsReadout'][1][k], ZImage, RowData['DetectorPos'], Seperation)
            
            if np.all(HittingPoints == -9999):
                DetectorCounts = -9999
                break
            
            else:  
                Iind = int(round((HittingPoints[0][0] + ImageLayerSize[0]/2)*(ProjectionPixel[0]-1)/ImageLayerSize[0])) 
                Jind = int(round((HittingPoints[0][1] + ImageLayerSize[1]/2)*(ProjectionPixel[1]-1)/ImageLayerSize[1]))
                
                if Iind>0 and Jind>0 and Iind<len(DetectorCounts[0]) and Jind<len(DetectorCounts[1]): 
                        DetectorCounts[Iind,Jind] += 1 
    
    tictoc = np.round(time(),1) - begin_time
    print('It took ', tictoc,' s to analyse the RowData from ', RowData['FileName'])
    
    return DetectorCounts 



# -------------------------- Clustering Algorithms -----------------------------------------



def ClusterAlgorithm(Data, Threshold, StartIndex):
    #Expects 2D array
    
    #--------------------------------------------------------------------
    # Isolates regions in Data above the threshold given the StartIndex
    #--------------------------------------------------------------------
    
    def Check(ActiveIndices, Shift, Zeros, Temp): # Checks to see if the pixel Shift-ed relative to ActiveIndices[j][k] is above Threshold
        
        if 0 <= ActiveIndices[j][k][0] + Shift[0] < Shape[0] and 0 <= ActiveIndices[j][k][1] + Shift[1] < Shape[1]:
            if Data[ActiveIndices[j][k][0] + Shift[0], ActiveIndices[j][k][1] + Shift[1]] < Threshold:   
                if Data[ActiveIndices[j][k][0] + Shift[0], ActiveIndices[j][k][1] + Shift[1]] == 0:
                    Zeros += 1
                
            elif not [ActiveIndices[j][k][0] + Shift[0], ActiveIndices[j][k][1] + Shift[1]] in ActiveIndices[j-1]:
                if not [ActiveIndices[j][k][0] + Shift[0], ActiveIndices[j][k][1] + Shift[1]] in Temp:
                    Temp.append([ActiveIndices[j][k][0] + Shift[0], ActiveIndices[j][k][1] + Shift[1]])
            
        return Zeros, Temp
    
    ActiveIndices = [[StartIndex]]  
    Shape = np.shape(Data)
    ClusteredData = np.zeros((Shape[0], Shape[1]))
    
    i = 0
    j = 0
    while i == 0:
        Temp = []
        LayerZeros = 0
        
        for k in range(len(ActiveIndices[j])):   
                                             
            Zeros = 0
            
            Zeros, Temp = Check(ActiveIndices, [1,0], Zeros, Temp)
            Zeros, Temp = Check(ActiveIndices, [0,1], Zeros, Temp)
            Zeros, Temp = Check(ActiveIndices, [-1,0], Zeros, Temp)
            Zeros, Temp = Check(ActiveIndices, [0,-1], Zeros, Temp)
            
            if Zeros >= 2:
                LayerZeros += 1
        
            if len(Temp) > 1000 or k > 1000:
                i += 1
                break
        
        if len(Temp) == 0:
            i += 1
            
        elif LayerZeros <= len(ActiveIndices[j]):
            ActiveIndices.append(Temp)
            j += 1
        
        else: 
            i += 1
        
    for layer in ActiveIndices:
        for pixel in layer:
            ClusteredData[pixel[0], pixel[1]] = Data[pixel[0], pixel[1]]
    
    ClusterDict = {'Active Indices':ActiveIndices, \
                   'Clustered Array':ClusteredData, \
                   'Start Index':StartIndex, \
                   'Threshold Value':Threshold}
        
    return ClusterDict



def LocalMaxIndices(Data, LocalCutoff, Divide):
    #Expects 2D arrays for Data
        
    #--------------------------------------------------------------------
    # Finds extrema above the LocalCutoff in the regions or Data 
    # specified by Divide and outputs the corresponding Index and Value
    #--------------------------------------------------------------------
    
    Shape = np.shape(Data)
    
    DivX = Divide[0]
    DivY = Divide[1]
    
    IndListX = np.linspace(0, Shape[0], num=DivX).astype(int)
    IndListY = np.linspace(0, Shape[1], num=DivY).astype(int)
    
    Value = []
    Index = []
    LayerMaxima = np.max(Data)
    
    for j in range(DivX -1):
        for k in range(DivY - 1):
            Max = np.max(Data[IndListX[j]:IndListX[j+1],IndListY[k]:IndListY[k+1]])  
            
            if Max > LayerMaxima*LocalCutoff: #LocalCutoff:
                Temp = []
                for l in range(len(np.where(Data[IndListX[j]:IndListX[j+1],IndListY[k]:IndListY[k+1]] == Max)[0])): 
                    TempInd = [IndListX[j] + np.where(Data[IndListX[j]:IndListX[j+1],IndListY[k]:IndListY[k+1]] == Max)[0][l], \
                               IndListY[k] + np.where(Data[IndListX[j]:IndListX[j+1],IndListY[k]:IndListY[k+1]] == Max)[1][l]] 
                                                                                                    
                    if not TempInd in Temp:
                        Temp.append(TempInd) 
        
                Value.append(Max)
                Index.append(Temp)
    
    return Value, Index



def ClusterMaxima(Data, Value, Index, PercentCutoff):
    #Expects 2D array
    
    #--------------------------------------------------------------------
    # Runs the ClusterAlgorithm at each element of Index and provided the 
    # results have more than five nonzero pixels, then combines them to  
    # form one image
    #--------------------------------------------------------------------
    
    AllClusters = []
    Areas = []
    ActiveIndices = []
    
    LayerMaxima = np.max(Data)
    
    #PercentCutoff = np.std(Data) / LayerMaxima
    
    for i in range(len(Index)):
        
        for l in range(len(Index[i])):
            Dict = ClusterAlgorithm(Data, LayerMaxima*PercentCutoff, [Index[i][l][0],Index[i][l][1]]) 
            
            if np.count_nonzero(Dict['Clustered Array']) > 5: # Requires two full layers of clustered Pixels 
                AllClusters.append(Dict['Clustered Array']) # ( 1 start pixel + 4 adjacent to it) to be considered significant
                n = np.count_nonzero(Dict['Clustered Array'])
                
                for layer in Dict['Active Indices']:
                    for ind in layer:
                        if not ind in ActiveIndices and len(Dict['Active Indices']) > 1:
                            ActiveIndices.append(ind)
                
                Areas.append(n)
        
    AllClusters = np.transpose(AllClusters)
    
    Shape = np.shape(Data)
    
    Maximized = np.zeros((Shape[1],Shape[0]))
    
    if len(AllClusters) != 0:
        for k in range(Shape[0]):
            for l in range(Shape[1]):  
                Maximized[k][l] = np.max(AllClusters[k][l]) 
                    
    Maximized = np.transpose(Maximized) 
    
    LayerDict = {'Active Indices':ActiveIndices, \
                 'Clustered Array':Maximized, \
                 'Start Value':Value, \
                 'Start Index':Index, \
                 'Index Cluster Areas':Areas, \
                 'Final Area':len(ActiveIndices)}
    
    return LayerDict



def FillCluster(Array,ClArray):
    '''
    takes clustered array and fills gaps in the solid object 
    with corresponding data from the original array
    
    input:
        Array [float] - [2D] unclustered array
        ClArray [float] - [2D] clustered array
    
    output:
        Fill [float] - [2D] filled array
    '''
    Shape = np.shape(ClArray)
    Fill = np.copy(ClArray)
    
    ind = np.array([[0,1],[0,-1],[1,0],[-1,0]])
    
    for i in range(1,Shape[0]-1):
        for j in range(1,Shape[1]-1):
            if ClArray[i,j] == 0:
                Non0 = [0 if Fill[i+k[0],j+k[1]] == 0 else 1 for k in ind]
                
                if np.count_nonzero(Non0) >= 3:
                    if Array[i,j] > 0:
                        Fill[i,j] += Array[i,j]
    
    return Fill



def TruncateFourier(DCounts,SigFrac,ThreeD):
    '''
    Truncates high frequency signals in DetectorCounts to reduce noise
    Uses real fast fourier transform and truncates the array, then
    performs inverse transform to recover data with reduced noise.
    
    input:
        DetectorCounts [float] - real, 2D or 3D array
        SigFrac [float] - in (0,1) fraction of frequencies kept
        ThreeD [bool] - True if DCounts is 3D
        
    output:
        Trunc [float] - 2D array of truncated signals
    '''
    def Truncate(DetectorCounts):
        
        Shape = np.shape(DetectorCounts)
        
        FFT2 = np.fft.rfft2(DetectorCounts)
        
        FFT2[int(np.shape(FFT2)[0] * SigFrac):,int(np.shape(FFT2)[1] * SigFrac):] = 0
        
        Trunc = np.fft.irfft2(FFT2) # loses a column (changes dimensions)
        
        #Trunc = np.concatenate((Trunc,np.zeros((Shape[0],1))), axis=1) #restores original dimensions
                                                            #by concatenating column of zeros
        return Trunc
    
    if ThreeD == True:
        Trunc = np.transpose(([Truncate(DCounts[:,:,i]) for i in range(np.shape(DCounts)[2])]),[1,2,0])
    
    else:
        Trunc = Truncate(DCounts)
        
    return Trunc




def ClusterArray(Array):
    '''
    Clusteres a 2D array by both FFT and standard methods (filling both)
    
    input:
        Array [float] - array to be clustered
        
    output:
        ClusterDict [Dict] - cluster info by standard method
        FFTClusterDict [Dict] - cluster info by fft method
    '''
    def Fill(Original, Cluster):
    
        gate, FilledCounts = np.zeros(np.shape(Original)), np.ones(np.shape(Original))
    
        while np.any(gate != FilledCounts): #Runs FillCluster until it doesn't do anything
            gate = FilledCounts
    
            FilledCounts = FillCluster(Original, Cluster)
    
        return FilledCounts
    
    Cutoff = np.std(Array) / np.max(Array)
    
    Value, Index = LocalMaxIndices(Array, Cutoff, Divide)
    
    ClusterDict = ClusterMaxima(Array, Value, Index, Cutoff)
    
    FFTCounts = TruncateFourier(Array, SigFrac, False)
    
    FFTClusterDict = ClusterMaxima(FFTCounts, Value, Index, 4 * Cutoff)
    #remove 4 as a free parameter some how?
    
    FilledCounts = Fill(Array, ClusterDict['Clustered Array'])
    FilledFFT = Fill(FFTCounts, FFTClusterDict['Clustered Array'])
    
    ClusterDict['Filled Array'] = FilledCounts
    FFTClusterDict['Filled Array'] = FilledFFT
    
    return ClusterDict, FFTClusterDict



def Cluster(ReadDict):
    '''
    Finds local maxima in the array and runs clustering algorithm at each for
    both the original array and the truncated fourier version
    
    input:
        ReadDict [dict] - dictionary of data read and processed from external file
        detdepth [float] - [m] depth of the bottom of detector
        
    output:
        AnalysisDict [dict] - clustering data with and without fft smoothing
    '''
    ClusterList = []
    FFTClusterList = [] 
    OriginalList = []

    for i in range(len(ReadDict['Subtracted Count List'])):
        objdepth = 2.5 #[m] depth of bottom of object
        detdepth = ReadDict['Parameters'][:,0] * np.cos(ReadDict['Parameters'][:,1])
        k = np.array([objdepth * ProjectionPixel[2] / dep for dep in detdepth],dtype=int)[i]
        
        ClusterDict, FFTClusterDict = ClusterArray(ReadDict['Subtracted Count List'][i][:,:,k])
        
        OriginalList.append(ReadDict['Subtracted Count List'][i][:,:,k])        
        ClusterList.append(ClusterDict)
        FFTClusterList.append(FFTClusterDict)
        
    AnalysisDict = {'Original':OriginalList, \
                    'Cluster':ClusterList, \
                    'FFT Cluster':FFTClusterList}
    
    return AnalysisDict



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
    


#---------------------------------- Get Data --------------------------------------



def ReadDataFiles():
    
    Parameters = load('Parameters.joblib')
    
    Shape = np.shape(Parameters)
    
    XPositions = [round(np.sin(par[1])*par[0]*100) for par in Parameters] #[cm]
    YPositions = np.zeros(Shape[0])
    ZPositions = [round(np.cos(par[1])*par[0]*100,2) for par in Parameters] #[cm]
    
    Seperations = [par[2] for par in Parameters]
    
    RealFiles = ['RDR_{}m_{}rad_{}cm.out'.format(par[0],par[1],par[2]) for par in Parameters]
    SkyFiles = ['RDS_{}m_{}rad_{}cm.out'.format(par[0],par[1],par[2]) for par in Parameters]
    
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
    
        DCS = PazAnalysis(RDSky,Seperations[i],Iterate,ZPositions[i]) 
        DCR = PazAnalysis(RDReal,Seperations[i],Iterate,ZPositions[i])
        
        if sign == 'True':
            DCdatPlus = (DCS-DCR) 
            
        else:
            DCdatPlus = (DCR-DCS) 
        
        DCdatPlus[DCdatPlus < 0] = 0
        
        RowSkyList.append(RDSky)
        RowRealList.append(RDReal)
        SkyCountList.append(DCS) 
        RealCountList.append(DCR)
        
        CountList.append(DCdatPlus)    
            
    ReadDict = {'Length':Shape[0], \
                'Iterate':Iterate, \
                'Parameters':Parameters, \
                'Row Sky List':RowSkyList, \
                'Row Real List':RowRealList, \
                'Sky Count List':SkyCountList, \
                'Real Count List':RealCountList, \
                'Subtracted Count List':CountList, \
                'Seperations':Seperations}
    
    return ReadDict


                    
#---------------------------------- Call Functions --------------------------------------



#-----------------------------------------------------------------
# Runs the imaging procedure on data specified by the user on
# the command line. Creates 2D images as well as tracked beams
# and estimates of the object sizes and shapes visually
#-----------------------------------------------------------------


ReadDictionary = ReadDataFiles()

dump(ReadDictionary,'ReadDictionary.joblib')
#ReadDictionary = load('ReadDictionary.joblib')

#ClusterDictionary  = Cluster(ReadDictionary) 

#ClustError = [ResolveError(ClusterDictionary['Cluster'][i]['Filled Array']) for i in range(ReadDictionary['Length'])] 
#FFTError = [ResolveError(ClusterDictionary['FFT Cluster'][i]['Filled Array']) for i in range(ReadDictionary['Length'])]

dump(ReadDictionary['Sky Count List'],'Sky_Arrays.joblib')
dump(ReadDictionary['Real Count List'],'Real_Arrays.joblib')
#dump(ClusterDictionary,'ClusterDict.joblib')
#dump(SubCounts,'SubCounts.joblib')
#dump(np.array([ClustError,FFTError]),'Error.joblib')

###change to only run on cluster layer specified by objdepth etc..

#--------------------------------------------------------------------------------------------



tictoc = np.round(time(),1) - begin_time

print('The total runtime of the program was: {} s'.format(tictoc))
