#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 00:16:15 2021

@author: keegan
"""

import numpy as np
from joblib import dump, load
from time import time
import math
import matplotlib.pyplot as plt

import glob
import os

from scipy import ndimage

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
ImageLayers = [500]                                            #[cm] height above bottom of detector for image slices
Iterate = False                                                #True ==> analysis will be done on ProjectionPixel[2] regularly spaced layers in decay volume
                                                               #False ==> analysis will be done only on ImageLayers

#---------------------
# Clustering 
#---------------------
Divide = [3, 3]                                               #[X divisions, Y divisions] of image layers in LocalMaxIndices() and LayerCluster()
SigFrac = 0.2                                                   #Percent of low fourier frequencies kept 
                                                                

#---------------------
# Run Options
#---------------------          


#---------------------
# Readout Error
#---------------------     
#STD = 0 #[MeV]                                             # Standard deviations of added gaussian errors
#key = 'Calibration'                                         # type of error



#%%
# --------------------------------- Translation of Gilad's Code -----------------------------------------



def ReadRowDataFileFastest(FileName,DetectorPos, STD=0, spread=0, Threshold=0): #must use 'path/Filename' (including quoatations) in function call 
    
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

    Decalib = 2 * spread * np.random.ranf([4,NumOfBars]) - spread # [%] random callibration error for each bar
    threshold = np.random.ranf([4,NumOfBars]) * Threshold # [MeV] Random threshold 
    
    for i in range(n-1): # -1 to avoid error from length of last event 
        Bar = []
        Length = []
        
        for l in range(12):
            if 12 <= len(G[EventWordInd[i]+l]) <= 13: #Bounds depend delicately on precision of path lengths
                Bar.append(G[EventWordInd[i]+l][0:3])
                Length.append(G[EventWordInd[i]+l][4:])
    
        
        BarFloat = np.float_(Bar).tolist() #Converts string elements to float
        
        '''Introduce errors to data'''
        
        Add = [np.random.normal(scale=STD,size=len(Length)),
               [np.float_(np.array(Length))[i] * \
                              Decalib[int(BarFloat[i]/100) - 1,int(BarFloat[i]%100)] \
                              for i in range(len(BarFloat))]]
        
        LengthFloat = np.float_(np.array(Length))+Add[0]+Add[1] # Introduce random gaussian and calibration erro
        
        cur_bars_thresh = [threshold[int(BarFloat[i]/100) - 1,int(BarFloat[i]%100)] \
                              for i in range(len(BarFloat))]
        
#        LengthToList = [LengthFloat[i] for i in range(len(LengthFloat)) if LengthFloat[i] >= cur_bars_thresh[i]] # introduce random cutoff threshold
        
        LengthToList = []
        rm_indices = []
    
        for i in range(len(LengthFloat)):
            if LengthFloat[i] >= cur_bars_thresh[i]:
                LengthToList.append(LengthFloat[i])
            
            else:
                rm_indices.append(i)
                
        BarFloat = list(np.delete(BarFloat,rm_indices))
        
        ''''''
              
        RowData['BarsReadout'][0].append(BarFloat)
        RowData['BarsReadout'][1].append(LengthToList)
        
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
                    continue
                
                else:  
                    try:
                        Iind = int(round((HittingPoints[0][0] + ImageLayerSize[0]/2)*(ProjectionPixel[0]-1)/ImageLayerSize[0])) 
                        Jind = int(round((HittingPoints[0][1] + ImageLayerSize[1]/2)*(ProjectionPixel[1]-1)/ImageLayerSize[1]))
                        
                        if Iind>0 and Jind>0 and Iind<len(DetectorCounts[0]) and Jind<len(DetectorCounts[1]): 
                                DetectorCounts[Iind,Jind,i] += 1 
            
                    except ValueError:
                        print(HittingPoints,'\n',RowData['BarsReadout'][0][k], RowData['BarsReadout'][1][k])
                        continue
                    
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
                    try:
                        Iind = int(round((HittingPoints[0][0] + ImageLayerSize[0]/2)*(ProjectionPixel[0]-1)/ImageLayerSize[0])) 
                        Jind = int(round((HittingPoints[0][1] + ImageLayerSize[1]/2)*(ProjectionPixel[1]-1)/ImageLayerSize[1]))
                        
                        if Iind>0 and Jind>0 and Iind<len(DetectorCounts[0]) and Jind<len(DetectorCounts[1]): 
                                DetectorCounts[Iind,Jind] += 1 
                            
                    except ValueError:
                        print(HittingPoints,'\n',RowData['BarsReadout'][0][k], RowData['BarsReadout'][1][k])
                        continue
    
    tictoc = np.round(time() - begin_time,1) 
    print('It took ', tictoc,' s to analyse the RowData from ', RowData['FileName'],'\n\n')
    
    return DetectorCounts 



#%%
# error analysis
    
class Error():
    
    def __init__(self,Fname):
        self.fname = Fname
    
    def analyse(self,Position,Seperation):
        
        self.Position = Position
        self.Seperation = Seperation
#        self.spread = spread
        
        self.RowData = ReadRowDataFileFastest(self.fname, Position, 0, 0, 0)
        self.Counts = PazAnalysis(self.RowData,Seperation,Iterate,Position[2],ImageLayers) 
        self.shape = np.shape(self.Counts)


    def propagate(self,STDs=[0],spread=[0],Threshold=[0]):
        
        self.spread = spread
        self.stds = STDs
        self.threshold = Threshold
        
        self.errs = [self.Counts]
        
        for STD in STDs: 
            for sprd in spread:
                for thrs in Threshold:
                    RD = ReadRowDataFileFastest(self.fname, self.Position, STD, sprd, thrs)
                    Counts = PazAnalysis(RD,self.Seperation,Iterate,self.Position[2],ImageLayers)
                    self.errs.append(Counts)
        
    
    def show_counts(self):
        
        try:
        
            for arr in self.errs:
                for i in range(np.shape(arr)[2]):
                    fig,ax = plt.subplots(dpi=100,figsize=(8,5))
                    
                    img = ndimage.gaussian_filter(arr[:,:,i],0.5)
                
                    ax.imshow(img)
                    plt.show()
             
        except Exception:   
            for arr in self.errs:
                fig,ax = plt.subplots(dpi=100,figsize=(8,5))
                
                img = ndimage.gaussian_filter(arr,0.5)
                
                ax.imshow(img)
                plt.show()
            
            
#    def make_dist(self):
#        
#        self.x = []
#        self.y = []
#        
#        
#        for arr in self.errs:
#            for i in range(np.shape(arr)[2]):
#                self.x.append(np.sum(arr[:,:,i],axis=0))
#                self.y.append(np.sum(arr[:,:,i],axis=1))
    
    
    def show_dist_position(self,labl=None):
        
        self.x = []
        self.y = []
        
#        z = self.Position[2]
        
        anglex = (np.arange(ProjectionPixel[0])-ProjectionPixel[0]/2) * ImageLayerSize[0] /(100 * ProjectionPixel[0]/2)
        angley = (np.arange(ProjectionPixel[1])-ProjectionPixel[1]/2) * ImageLayerSize[1] /(100 * ProjectionPixel[1]/2)
        
        for arr in self.errs:
            for i in range(np.shape(arr)[2]):
                self.x.append(np.sum(arr[:,:,i],axis=0))
                self.y.append(np.sum(arr[:,:,i],axis=1))
        
        #self.stds = [0.1,0.25,0.5] #[0.5,1,2]
        
        plt.rcParams.update({'font.size': 16})

        fig, ax = plt.subplots(dpi=150,nrows=len(err.errs),ncols=2,figsize=[16,5*len(err.errs)])

        [ax[i,0].bar(anglex,self.x[i],label=labl) for i in range(len(self.errs))]
        [ax[i,1].bar(angley,self.y[i],label=labl) for i in range(len(self.errs))]
        
        #np.hist
        
        for i in range(len(self.errs)):
            ax[i,0].text(np.amin(anglex),np.max(self.x[i])*0.05,\
            'Mean: {:.02g} \n STD: {:.02g}'.format(np.average(anglex,weights=self.x[i]), np.sqrt(np.average(anglex**2,weights=self.x[i]))))
            ax[i,1].text(np.amin(angley),np.max(self.y[i])*0.05,\
            'Mean: {:.02g} \n STD: {:.02g}'.format(np.average(angley,weights=self.y[i]), np.sqrt(np.average(angley**2,weights=self.y[i]))))
            
        [a.semilogy() for a in ax[:,0]]
        [a.semilogy() for a in ax[:,1]]
        
        ax[0,0].set_title('Predicted X Position')
        ax[0,1].set_title('Predicted Y Position')
        
        [ax[i,0].set_ylabel('Error STD = {} MeV'.format([0,*self.stds][i])) for i in range(len(self.errs))]
        [[ax[j,i].set_xlabel('Position [m]') for i in range(2)] for j in range(len(self.errs))]
        
        if labl != None:
            ax.legend()
        
        plt.show()
        
        
        
    def show_dist_slice(self,labl=None,title=None):
        
        self.x = []
        
        z = self.Position[2]
        
        anglex = np.arctan((np.arange(ProjectionPixel[0])-ProjectionPixel[0]/2)[100:] * ImageLayerSize[0] /(100 * ProjectionPixel[0]/2) / z)
        
        for arr in self.errs:
            for i in range(np.shape(arr)[2]):
                self.x.append(np.sum(arr[99:101,100:,i],axis=0))
        
        plt.rcParams.update({'font.size': 16})
        
        fig, ax = plt.subplots(dpi=100,nrows=len(self.errs),figsize=[8,5*len(self.errs)])
        
#        [ax[i].bar(anglex,self.x[i],label=labl) for i in range(len(self.errs))]
        [ax[i].plot(anglex,self.x[i],label=labl,marker='.') for i in range(len(self.errs))]
        
        for i in range(len(self.errs)):
            mean = np.average(anglex,weights=self.x[i])
            ax[i].text(np.amin(anglex),np.max(self.x[i])*0.15,\
            'Mean: {:.04g} \n STD: {:.04g}'.format(mean, np.sqrt(np.average((anglex-mean)**2,weights=self.x[i]))))
               
#        [a.semilogy() for a in ax[:]]
        
        ax[0].set_title('Predicted Position')
        
        n = 1
        ax[0].set_ylabel('Original')
        for i in range(len(self.stds)):
            for j in range(len(self.spread)):
                for k in range(len(self.threshold)):
                    ax[n].set_ylabel('STD = {} MeV, {} %, > {} Mev'.format([*self.stds][i],[*self.spread][j],[*self.threshold][k]))
                    n += 1
                
        [ax[j].set_xlabel('Angle [rad]') for j in range(len(self.errs))]
        
        if labl != None:
            ax.legend()
        plt.show()
            
        
    
    def show_dist_angle(self,length):
    
        
        self.bins = []
        self.hist = []
        
        for arr in self.errs:

            shape = np.shape(arr)    
        
            r_max = np.sqrt(shape[0]**2 + shape[1]**2)
        
            bins = np.linspace(0,r_max,length)
            hist = np.zeros(length)
        
            for i in range(shape[0]):
                for j in range(shape[1]):
                    r = np.sqrt(i**2 + j**2)
                    
                    ind = np.where(bins<r)[-1]-1           
        
                    hist[ind] += arr[i][j] 
        
                    
            bins = np.array(bins)
            hist = np.array(hist).flatten()
            
            fig, ax = plt.subplots(dpi=150,figsize=[8,5])
            
            inds = np.argsort(bins)
                
            hist = np.take_along_axis(hist,inds,axis=0)
            bins = np.take_along_axis(bins,inds,axis=0)
            
            ax.bar(bins,hist)
            plt.show()
            
            self.bins.append(bins)
            self.hist.append(hist)
        
    
    
'''

'''

if __name__ == '__main__':
    
    fname_list = ['Sun_Mar_28_11-44-41_2021/RDR_[0,0]cm_5m_25cm_0a_0b_0c_[0,0,0]m.out',
                  'Fri_Mar_19_23-48-16_2021/RDR_[0,0]cm_5m_25cm_0a_0b_0c_[0,0,0]m.out',
                  'Sun_Mar_28_15-39-44_2021/RDR_[0,0]cm_5m_25cm_0a_0b_0c_[0,0,0]m.out', 
                  'Sat_Apr__3_23-46-44_2021/RDR_[0,0]cm_5m_25cm_0a_0b_0c_[0,0,0]m.out',
                  'Sat_Mar_13_09-50-39_2021/RDR_[0,0]cm_5m_40cm_0a_0b_0c_[0,0,0]m.out',
                  'Mon_Apr__5_22-49-14_2021/RDR_[0,0]cm_5m_25cm_0a_0b_0c_[0,0,0]m.out']
    
    ob_list = []
    
    for name in fname_list:
        err = Error(name)
        err.analyse([0,0,5],40)
#        err.propagate([0.25,0.5,1,2])
        err.propagate([0],[0],[0.25,0.5,1])
        err.show_counts()
        err.show_dist_slice()
        
        ob_list.append(err)


