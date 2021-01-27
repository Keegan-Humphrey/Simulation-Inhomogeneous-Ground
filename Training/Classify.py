#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 15:19:38 2021

@author: keegan
"""

import numpy as np
import datetime
import glob
from time import time
import os
import numpy.random as rnd

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
from skimage.transform import rescale
 
from joblib import dump, load

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
    begin_time = time()
    
    Shape = np.shape(Data)

    #LearningData = np.zeros((Shape[0],Shape[1],Copies*Shape[2]))
    
    
    if ThreeD:
        Noisy_Data = [np.copy(Data[:,:,i]) for i in range(Shape[2]) for j in range(Copies)]
    
    else:
        Noisy_Data = [np.copy(Data) for i in range(Copies)]
    
    rnd.shuffle(Noisy_Data)
    
    Noisy_Data = np.transpose(Noisy_Data, axes=[1,2,0])
    
    Noisy_Data += np.array(np.abs(rnd.normal(0, std, Shape)), dtype=int)
    
    Noisy_Data += np.array(rnd.choice([-1, 1]) * rnd.random(Shape) * Noisy_Data * Percent_Error, dtype=int) #adds a uniformly distributed error of PercentError to data
        
    #tictoc = np.round(time() - begin_time, 1)
    #print('It took ', tictoc,' to preprocess')
    
    return Noisy_Data #LearningData



def PrepareData(SigFrac, Copies, Percent_Error, std):
    '''
    input:
        DetectorCounts [float] - real, 2D or 3D array
        SigFrac [float] - in (0,1) (list) fraction of frequencies kept
        std [float] - (must be > 0) standard deviation of noise introduced (mu=0)
        Copies [int] - number of copies introduced to the data set 
        Percent_Error [float] - (must be > 0) fractional error introduced
        ThreeD [bool] - True ==> data is three dimensional
    '''
    Directory = os.getcwd() + '/Air_1_0.5_100_m_True/'

    names = [filename[filename.find('RDR')+3:] for filename in glob.iglob(Directory+'RDR*.joblib', recursive=True)]
    
    RealFiles = [Directory+'RDR'+name for name in names] 
    SkyFiles = [Directory+'RDS'+name for name in names] 
    
    def Prep(DC,percent,std_dev):
        
        for i in range(len(names)):
            #Data = load(Files[i])
            Data = DC[i]
            
            Data = np.transpose([TruncateFourier(Data[:,:,j], frac, False) for frac in SigFrac for j in range(np.shape(Data)[2])],[1,2,0])
            
            Data = LearnPreProcess(Data, Copies, percent, std_dev, len(np.shape(Data))>2)
            
            if i == 0:
                Prepped = Data
                
            else:
                Prepped = np.concatenate([Prepped, Data], axis=2)
        
        return Prepped
    
    R_Data = np.array([load(RealFiles[i]) for i in range(len(names))])
    S_Data = np.array([load(SkyFiles[i]) for i in range(len(names))])
    
    #indices = {(i,j) for i in range(len(names)) for j in range(len(names)) if i != j}
    
    #print([np.any(S_Data[i]==S_Data[j]) for i, j in indices])
    
    #DCdat = R_Data - S_Data
    DCdat = S_Data - R_Data
    DCdat[DCdat < 0] = 0
        
    Prepped_R = Prep(R_Data, Percent_Error[0], std[0])
    Prepped_S = Prep(S_Data, Percent_Error[1], std[1])
    Prepped_D = Prep(DCdat, Percent_Error[0], std[0])    
    
    return Prepped_R, Prepped_S, Prepped_D
        


#%%


Prepped_R, Prepped_S, Prepped_D = PrepareData([1], 1, [0,0.5], [0.1,0.1])

NumObj = 1 # number of objects in the real images

TargetSky = np.zeros(np.shape(Prepped_S)[2])
TargetReal = np.ones(np.shape(Prepped_R)[2]) * NumObj
TargetSub = np.ones(np.shape(Prepped_R)[2]) * NumObj

LearningData = np.concatenate((Prepped_S, Prepped_R, Prepped_D), axis=2)
Target = np.concatenate((TargetSky, TargetReal, TargetSub))

LearningData = rescale(LearningData, 1/4, mode='reflect')

LearningData = np.transpose(LearningData, [2,0,1])


X_train, X_test, y_train, y_test = train_test_split(LearningData,
                                                    Target,
                                                    test_size=0.1,
                                                    shuffle=True,
                                                    random_state=42)


#Train with dense objects only
#use both dc and real to train on, hopefully it will work
# for very weak signal than training on the weak signal


#%%
#what follows is courtesy of:
#https://kapernikov.com/tutorial-image-classification-with-scikit-learn/

#%%
#Transform

from sklearn.base import BaseEstimator, TransformerMixin
 
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
        
        
        
from sklearn.linear_model import SGDClassifier
from sklearn.model_selection import cross_val_predict
from sklearn.preprocessing import StandardScaler
import skimage
 
# create an instance of each transformer
grayify = RGB2GrayTransformer()
hogify = HogTransformer(
    pixels_per_cell=(8, 8),
    cells_per_block=(2,2),
    orientations=9,
    block_norm='L2-Hys'
)
scalify = StandardScaler()
 
# call fit_transform on each transform converting X_train step by step
X_train_gray = grayify.fit_transform(X_train)
X_train_hog = hogify.fit_transform(X_train_gray)
X_train_prepared = scalify.fit_transform(X_train_hog)
 
print(X_train_prepared.shape)


#%% 
#Train and test

sgd_clf = SGDClassifier(random_state=42, max_iter=1000, tol=1e-3)
sgd_clf.fit(X_train_prepared, y_train)


#big_test = load(os.getcwd() + '/large run/big_real.joblib')
#big_test = rescale(big_test, 1/4, mode='reflect')
#big_test = np.transpose(big_test, [2,0,1])
#
#X_test = big_test
#y_test = np.ones(np.shape(big_test)[0])


X_test_gray = grayify.transform(X_test)
X_test_hog = hogify.transform(X_test_gray)
X_test_prepared = scalify.transform(X_test_hog)

y_pred = sgd_clf.predict(X_test_prepared)
#print(np.array(y_pred == y_test)[:25])
print('')
print('Percentage correct: ', 100*np.sum(y_pred == y_test)/len(y_test))


#%%
#Diagnose

from sklearn.metrics import confusion_matrix
 
cmx = confusion_matrix(y_test, y_pred)

print(cmx)


#%%

from sklearn.pipeline import Pipeline
from sklearn import svm
 
HOG_pipeline = Pipeline([
    ('grayify', RGB2GrayTransformer()),
    ('hogify', HogTransformer(
        pixels_per_cell=(8, 8),
        cells_per_block=(2,2),
        orientations=9,
        block_norm='L2-Hys')
    ),
    ('scalify', StandardScaler()),
    ('classify', SGDClassifier(random_state=42, max_iter=1000, tol=1e-3))
])
 
clf = HOG_pipeline.fit(X_train, y_train)
print('Percentage correct: ', 100*np.sum(clf.predict(X_test) == y_test)/len(y_test))

#%%
from sklearn.model_selection import GridSearchCV
 
param_grid = [
    {'hogify__orientations': [9],
    'hogify__cells_per_block': [(3, 3)],
    'hogify__pixels_per_cell': [(8, 8), (14, 14)]},
    {'hogify__orientations': [9],
     'hogify__cells_per_block': [(3, 3)],
     'hogify__pixels_per_cell': [(14, 14)],
     'classify': [
         SGDClassifier(random_state=42, max_iter=1000, tol=1e-3),
         svm.SVC(kernel='linear')]}
]

#%%

grid_search = GridSearchCV(HOG_pipeline,
                           param_grid,
                           cv=3,
                           n_jobs=-1,
                           scoring='accuracy',
                           verbose=1,
                           return_train_score=True)
 
grid_res = grid_search.fit(X_train, y_train)


#%%
# description of the best performing object, a pipeline in our case.
grid_res.best_estimator_

grid_res.best_score_

pp.pprint(grid_res.best_params_)

#%%

from sklearn.model_selection import RandomizedSearchCV

best_pred = grid_res.predict(X_test)
print('Percentage correct: ', 100*np.sum(best_pred == y_test)/len(y_test))

cmx_svm = confusion_matrix(y_test, best_pred)
plot_confusion_matrix(cmx, vmax1=225, vmax2=100, vmax3=12)

plot_confusion_matrix(cmx_svm, vmax1=225, vmax2=100, vmax3=12)
