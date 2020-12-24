import matplotlib.pyplot as plt
from joblib import load
import numpy as np

def PlotQuick(Data, ThreeD):   
    #Expects Boolean value for ThreeD
    
    #--------------------------------------------------------------------
    # Plots each image layer in Data as a list of colormapped images
    #--------------------------------------------------------------------    
    
    def PlotQuick2D(Data):
        
        fig = plt.figure(figsize=(15,15))
        ax = fig.add_subplot(1,1,1)
        
        #ax.set_title('PazVsReal')
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

sky = load('Sky_Arrays.joblib')
real = load('Real_Arrays.joblib')

#PlotQuick(Dict['Subtracted Count List'][0], Dict['Iterate'])

DCdatPlus = (sky[0]-real[0])
#DCdatPlus = (real[0]-sky[0])
DCdatPlus[DCdatPlus < 0] = 0

PlotQuick(sky[0], True)
PlotQuick(real[0], True)
PlotQuick(DCdatPlus, True)
