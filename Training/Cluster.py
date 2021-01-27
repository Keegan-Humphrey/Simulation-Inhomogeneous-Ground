#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  9 20:22:07 2021

@author: keegan
"""
import numpy as np


#%%
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


#%%


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

