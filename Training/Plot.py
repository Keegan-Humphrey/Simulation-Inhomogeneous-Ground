import matplotlib.pyplot as plt
from joblib import load
import numpy as np
import glob
import os

def PlotQuick(Data, ThreeD, Title):   
    #Expects Boolean value for ThreeD
    
    #--------------------------------------------------------------------
    # Plots each image layer in Data as a list of colormapped images
    #--------------------------------------------------------------------    
    
    def PlotQuick2D(Data):
        
        fig = plt.figure(figsize=(15,15))
        ax = fig.add_subplot(1,1,1)
        
        ax.set_title(Title)
        plt.imshow(Data, alpha = 0.5)
        
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.patch.set_alpha(0)
        ax.set_frame_on(False)
        plt.colorbar(orientation='vertical')
        
        plt.show()
   
    
    if ThreeD == True:
        Shape = np.shape(Data)
        
        for i in range(Shape[2]):
            PlotQuick2D(Data[:,:,i])
    
    else:
        PlotQuick2D(Data)

#Dict = load('ReadDictionary.joblib')

Directory = os.getcwd() + '/Tue_Jan_26_09-44-54_2021/' #load('time.joblib')

names = [filename[filename.find('RDR')+3:] for filename in glob.iglob(Directory+'RDR*.joblib', recursive=True)]
    

#names = [filename[3:] for filename in glob.iglob('RDR*.joblib', recursive=True)]
#SkyFiles = [filename for filename in glob.iglob('RDS*.joblib', recursive=True)]

#names = ['_30m_0.63rad_25cm.joblib']
   
#sky = load('Sky_Arrays.joblib')
#real = load('Real_Arrays.joblib')

#PlotQuick(Dict['Subtracted Count List'][0], Dict['Iterate'])

for i in range(len(names)):
    sky = np.array(load(Directory + 'RDR' + names[i]), dtype=int)
    real = np.array(load(Directory + 'RDS' + names[i]), dtype=int)
    
    #DCdatPlus = (real-sky)
    DCdatPlus = (sky-real)
    DCdatPlus[DCdatPlus < 0] = 0
    
    for j in range(np.shape(sky)[2]):
        PlotQuick(sky[:,:,j], False, 'RDS'+names[i])
        PlotQuick(real[:,:,j], False, 'RDR'+names[i])
        PlotQuick(DCdatPlus[:,:,j], False, 'DCR'+names[i])
