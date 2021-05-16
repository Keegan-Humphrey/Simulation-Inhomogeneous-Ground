#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 15:19:38 2021

@author: keegan
"""

import numpy as np
import matplotlib.pyplot as plt
import datetime
import glob
from time import time
import os
import numpy.random as rnd
#import sys
import warnings

#from Plot import PlotQuick

from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier
from sklearn.linear_model import SGDRegressor, SGDClassifier
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import AffinityPropagation
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import train_test_split
from skimage.feature import hog
from skimage.io import imread
from skimage.transform import rescale, resize
from skimage.color import rgb2gray
 
from joblib import dump, load


from sklearn.base import BaseEstimator, TransformerMixin
        
from sklearn.linear_model import SGDClassifier
from sklearn.model_selection import cross_val_predict
from sklearn.preprocessing import StandardScaler
import skimage

from sklearn.pipeline import Pipeline
from sklearn import svm

from sklearn.model_selection import GridSearchCV

from scipy import ndimage
 

#import matplotlib.pyplot as plt
   
#%%
# Noise Manipulation

def TruncateFourier(DCounts, SigFrac, ThreeD):
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
        
        #Shape = np.shape(DetectorCounts)
        
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



def LearnPreProcess(Data, Copies, Percent_Error, std, ThreeD): 
    #Expects 3D array for Data    
    
    '''
    Adds noise and error to Data
    
    input:
        std [float] - (must be > 0) standard deviation of noise introduced
        Copies [int] - number of copies introduced to the data set 
        Percent_Error [float] - (must be > 0) fractional error introduced
        ThreeD [bool] - True ==> data is three dimensional
    '''
    #begin_time = time()
    
    Shape = np.shape(Data)

    #LearningData = np.zeros((Shape[0],Shape[1],Copies*Shape[2]))
    
    
    if ThreeD:
        Noisy_Data = [np.copy(Data[:,:,i]) for i in range(Shape[2]) for j in range(Copies)]
        
    else:
        Noisy_Data = [np.copy(Data) for i in range(Copies)]
    
    rnd.shuffle(Noisy_Data)
    
    Noisy_Data = np.transpose(Noisy_Data, axes=[1,2,0])
    
    Noisy_Data += np.array(np.abs(rnd.normal(0, std, \
                        (Shape[0],Shape[1],Shape[2]*Copies))), dtype=int) #adds noise
    
    Noisy_Data += np.array(rnd.choice([-1, 1],size=(Shape[0],Shape[1],Shape[2]*Copies)) \
                           * rnd.random((Shape[0],Shape[1],Shape[2]*Copies)) * \
                           Noisy_Data * Percent_Error, dtype=int) #adds a uniformly distributed error of PercentError to data
        
    #tictoc = np.round(time() - begin_time, 1)
    #print('It took ', tictoc,' to preprocess')
    
    return Noisy_Data #LearningData



def PrepareData(Dir_names, SigFrac, Copies, Percent_Error, std):
    '''
    input:
        # DetectorCounts [float] - real, 2D or 3D array
        SigFrac [float] - in (0,1) (list) fraction of frequencies kept
        std [float] - (must be > 0) standard deviation of noise introduced (mu=0)
        Copies [int] - number of copies introduced to the data set 
        Percent_Error [float] - (must be > 0) fractional error introduced
        # ThreeD [bool] - True ==> data is three dimensional
    '''
    def Prep(DC,percent,std_dev):
        
        for i in range(np.shape(DC)[0]):
            #Data = load(Files[i])
            Data = DC[i]
            
            try:
                Data = np.transpose([TruncateFourier(Data[:,:,j], frac, False) for frac in SigFrac for j in range(np.shape(Data)[2])],[1,2,0])
            
            except IndexError:
                Data = np.transpose([TruncateFourier(Data, frac, False) for frac in SigFrac],[1,2,0])
                
            Data = LearnPreProcess(Data, Copies, percent, std_dev, len(np.shape(Data))>2)
            
            if i == 0:
                Prepped = Data
                
            else:
                Prepped = np.concatenate([Prepped, Data], axis=2)
        
        return Prepped
    
    
    Directories = [Dir for Dir in Dir_names]
    
#    print(Directories)
    
    RealFiles = []
    SkyFiles = []
       
    for Directory in Directories:
        RealFiles.extend([filename for filename in glob.iglob(Directory+'RDR*.joblib', recursive=True)])
        SkyFiles.extend([filename for filename in glob.iglob(Directory+'RDS*.joblib', recursive=True)])

    
    
    R_Data = np.array([load(RealFiles[i]) for i in range(len(RealFiles))])
    S_Data = np.array([load(SkyFiles[i]) for i in range(len(SkyFiles))])
    
    #indices = {(i,j) for i in range(len(names)) for j in range(len(names)) if i != j}
    
    #print([np.any(S_Data[i]==S_Data[j]) for i, j in indices])
    
    #DCdat = R_Data - S_Data
    DCdat = S_Data - R_Data
    DCdat[DCdat < 0] = 0
        
#    PlotQuick(DCdat[0],False,'',False)
    
    Prepped_R = Prep(R_Data, Percent_Error[0], std[0])
    Prepped_S = Prep(S_Data, Percent_Error[1], std[1])
    Prepped_D = Prep(DCdat, Percent_Error[0], std[0])    
    
    return Prepped_R, Prepped_S, Prepped_D


        
#%%
#this cell (and much of what follows) is courtesy of:
#https://kapernikov.com/tutorial-image-classification-with-scikit-learn/

class RGB2GrayTransformer(BaseEstimator, TransformerMixin):
    """
    Convert an array of RGB images to grayscale
    """
 
    def __init__(self):
        pass
 
    def fit(self, X, y=None):
        """returns itself"""
        return self
 
    def transform(self, X, y=None):
        """perform the transformation and return an array"""
        return np.array([skimage.color.rgb2gray(img) for img in X])
 
 
class HogTransformer(BaseEstimator, TransformerMixin):
    """
    Expects an array of 2d arrays (1 channel images)
    Calculates hog features for each img
    """
 
    def __init__(self, y=None, orientations=9,
                 pixels_per_cell=(8, 8),
                 cells_per_block=(3, 3), block_norm='L2-Hys'):
        self.y = y
        self.orientations = orientations
        self.pixels_per_cell = pixels_per_cell
        self.cells_per_block = cells_per_block
        self.block_norm = block_norm
 
    def fit(self, X, y=None):
        return self
 
    def transform(self, X, y=None):
 
        def local_hog(X):
            return hog(X,
                       orientations=self.orientations,
                       pixels_per_cell=self.pixels_per_cell,
                       cells_per_block=self.cells_per_block,
                       block_norm=self.block_norm)
 
        try: # parallel
            return np.array([local_hog(img) for img in X])
        except:
            return np.array([local_hog(img) for img in X])
        




#%%
# get data



### try again with prep stuff outside again
## try without randomizing?? (preprocessing)


class make_dist:
    
    def __init__(self,dir_names):
        
        warnings.filterwarnings("ignore")
        
        self.dir = dir_names
        
        self.LearningData = None
        self.Target = None
        
        self.grayify = RGB2GrayTransformer()
        self.hogify = HogTransformer(pixels_per_cell=(8, 8),
                                     cells_per_block=(2,2),
                                     orientations=9,
                                     block_norm='L2-Hys')
        self.scalify = StandardScaler()
    
        self.Scores = None
        self.grid_res = None
        
        self.seed = 42 # "the answer to life, the universe, and everything"
        self.sgd = SGDClassifier(random_state=self.seed, max_iter=1000, tol=1e-3)
        self.clf = Pipeline([('grayify', self.grayify),
                             ('hogify',self.hogify),
                             ('scalify', self.scalify),
                             ('classify', self.sgd)])
    
        
        
    def get_data(self):
        
#        Prepped_R, Prepped_S, Prepped_D = PrepareData(self.dir, [1], 2, [0.05,0.05], [0.4,0.4])
        #Prepped_R, Prepped_S, Prepped_D = PrepareData(dir_names, [1], 2, [0.05,0.05], [0.4,0.4])
        Prepped_R, Prepped_S, Prepped_D = PrepareData(self.dir, [1], 1, [0.0,0.0], [0.0,0.0])
        
        NumObj = 1 # number of objects in the real images
        
        TargetSky = np.zeros(np.shape(Prepped_S)[2]) 
        TargetReal = np.ones(np.shape(Prepped_R)[2]) * NumObj
        #TargetSub = np.ones(np.shape(Prepped_D)[2]) * NumObj 
        
        LearningData = np.concatenate((Prepped_S, Prepped_R),axis=2)#, Prepped_D), axis=2)
        self.Target = np.concatenate((TargetSky, TargetReal))#, TargetSub))
        
        LearningData = rescale(LearningData, 0.25, mode='reflect')
        
        self.LearningData = np.transpose(LearningData, [2,0,1])
        
        return self.LearningData, self.Target
    
    
    
    def optimize(self):
        
        param_grid = [{'classify__tol':np.power(np.ones(3)*10,[-1,-3,-5]),
                       'classify__penalty':['l2', 'l1', 'elasticnet']}] 
                       
#                       'hogify__orientations': [5,10,15],
#                       'hogify__cells_per_block': [(3, 3),(2, 2),(1, 1)],
#                       'hogify__pixels_per_cell': [(7,7),(8, 8),(9, 9)]}]
                    
        
                             
        grid_search = GridSearchCV(self.clf,
                       param_grid,
                       cv=None,
                       n_jobs=1,
                       scoring='accuracy',
                       verbose=1,
                       return_train_score=True)
        
        grid_res = grid_search.fit(*self.split(self.LearningData,self.Target)[::2])
        
        #use mean values of the array data for grid search?
        
        #print(grid_res.best_params_)
        
        self.grid_res = grid_res
        
        self.clf = self.grid_res.estimator
        

    def preview(self):
        
        for i in range(0,np.shape(self.LearningData)[0],4):
#            for j in range(2):
            
                j = 0    
            
                plt.figure()
                
                if j == 1:
#                    img = ndimage.gaussian_filter(self.LearningData[i],0.5)
                    plt.imshow(self.LearningData[i])
                    print('denoised')
                
                else:
                    plt.imshow(self.LearningData[i])
                
                plt.show()
                    
                print(self.Target[i])
        
        
        
    def split(self,X,y):
        
        return train_test_split(X,
                                y,
                                test_size=0.2,
                                shuffle=True,
                                random_state=self.seed)
        
        
        
    def get_scores(self, verbose=True, reps=100):    
    
        Scores = []
        
        p_sky = []
        p_real = []
        
        real = np.loadtxt('real_access.txt',skiprows=1)
        sky = np.loadtxt('sky_access.txt',skiprows=1)
        
        real = resize(real,np.shape(dist.LearningData[0]))
        sky = resize(sky,np.shape(dist.LearningData[0]))

        
        for i in range(reps):
            
            X, y = self.get_data()
            
            self.seed = np.random.randint(314159)
            
            X_train, X_test, y_train, y_test = self.split(X, y)
            
            X_train_prepared = X_train
            
            # call fit_transform on each transform converting X_train step by step
#            X_train_gray = self.grayify.fit_transform(X_train)
#            X_train_hog = self.hogify.fit_transform(X_train_gray)
#            X_train_prepared = self.scalify.fit_transform(X_train_hog)
#             
#            try:
#                self.clf.set_params(random_state=self.seed)
#                
#            except ValueError:
            
            self.clf.set_params(classify__random_state=self.seed)
        
            self.clf.fit(X_train_prepared, y_train)

            scores = cross_val_score(self.clf, X_train_prepared, y_train, cv=5)
            
            Scores.append(scores)
            
            p_real.append(self.clf.predict([rgb2gray(real)]))
            p_sky.append(self.clf.predict([rgb2gray(sky)]))
             
            
            if verbose:
                if i % int(reps / 10) == 0:
                    print(str(int(i * 100 / reps)) + ' %')
            
        self.Scores = np.array(Scores).flatten()
        
        print('Mean sky score: ',np.mean(p_sky))
        print('Mean real score: ',np.mean(p_real))
        
        return Scores
     
      
    def show_dist(self,labl=None):
        
        fig, ax = plt.subplots(dpi=100,nrows=1,figsize=[8,5])
        ax.hist(self.Scores,bins=15,range=(0,1),label=labl)
        
        #ax[1].hist(Scores_0,bins=15,range=(0,1),label='Large Angles (beta=0)')
        #ax[2].hist(Scores_45,bins=15,range=(0,1),label='Large Angles (beta=-45)')
        
        ax.text(0,len(self.Scores)*0.05,'Mean: {:.02g} \n STD: {:.02g}'.format(np.mean(self.Scores), np.std(self.Scores)))
        
        ax.set_title('Estimator Accuracy')
        ax.set_xlabel('Accuracy [%]')
        
        if labl != None:
            ax.legend()
        plt.show()


Directory_names = [['/Users/keegan/Desktop/Research/visualisation/useful_data/Sat_Apr__3_06-53-14_2021/', \
                    '/Users/keegan/Desktop/Research/visualisation/useful_data/Sun_Feb_14_20-16-48_2021/b=0/',
                    '/Users/keegan/Desktop/Research/visualisation/useful_data/Tue_Feb__9_09-30-18_2021/']]


#Directory_names = [['/Users/keegan/Desktop/Research/visualisation/useful_data/Sun_Feb_14_20-16-48_2021/b=0/', \
#'/Users/keegan/Desktop/Research/visualisation/useful_data/Tue_Feb__9_09-30-18_2021/'], \
#['/Users/keegan/Desktop/Research/visualisation/useful_data/Sun_Feb_14_20-16-48_2021/b=-45/', \
#'/Users/keegan/Desktop/Research/visualisation/versions/simulation_1.0.2/Training/Sat_Feb_27_21-57-18_2021/']]



if __name__ == '__main__':
    
#    Directory_names = [['/Users/keegan/Desktop/Research/visualisation/Gilad_Builds_and_data/data/']] 

    dists = []
    
    for Directory in Directory_names:
        dist = make_dist(Directory)
        
        try:
            dist.get_data()
        
        except UnboundLocalError:
            print('Error, no files in that directory')
            continue
            
#        dist.preview()
        dist.optimize()
        dist.get_scores(verbose=True)
        dist.show_dist()
        
        dists.append(dist)

    
    
    


    
    
    
    
    
    
    
    
    
    
    
    
    

