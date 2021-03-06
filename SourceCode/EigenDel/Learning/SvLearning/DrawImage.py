import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import spline
import math

def WaveletTransfer(arryY):    
    MexicoHatBaseData = []
    ibaseNum = 16
    arryWTY = []
    iWindowSize = 60  # real one is iWindowSize*2 = 72 
    for i in range(iWindowSize*-1, iWindowSize + 1):  #(-36, 37): # get the value from -17<--->+17    
        MexicoHatBaseData.append( (1-float(i*i)/ibaseNum) / math.exp(float(i*i)/(2*ibaseNum)) );
      
    refData = [0]*(iWindowSize * 2 + 1)                
    valueY = 0.0;
    
    iFirstMaxPos = 0;
    iEndMinPos = len(arryY)-1;
    
    #print(len(arryY))
    for i in range(0, len(arryY)):    
        for j in range(iWindowSize*-1, iWindowSize + 1): #(-36, 37):        
            if (i+j)<iFirstMaxPos or (i+j)>iEndMinPos:
                #print(i-j)            
                refData[j+iWindowSize] = arryY[i-j]           
            else:
                refData[j+iWindowSize] = arryY[i+j]
                        
        # get the y value from sum product methord
        for j in range(0, iWindowSize*2+1):        
            valueY += MexicoHatBaseData[j]*refData[j];
    
        arryWTY.append(valueY)                                
        valueY = 0.0; # reset this temp summation values
    nparryY = np.array(arryWTY)
    return nparryY    

def DrawPics(x, y, filePath, status):    
    xnew = np.linspace(x.min(), x.max(), 100) #300 represents number of points to make between T.min and T.max 
    y_smooth = spline(x, y, xnew)
    #f2 = interp1d(x, y, kind='cubic')
    
          
    fig_size = plt.rcParams["figure.figsize"]
    fig_size[0] = 12
    fig_size[1] = 9
      
    # plotting the points
    plt.subplot(211)
    plt.ylim(0, y.max() + 1) 
    plt.xlim(0,x.max())
    plt.plot(x, y, color='green', linestyle='dashed', linewidth = 2, 
             marker='o', markerfacecolor='blue', markersize=8)
      
    # setting x and y axis range 
    plt.ylim(0, y.max() + 1) 
    plt.xlim(x.min()*1.1, x.max()*1.1) 
    # naming the x axis 
    #plt.xlabel('x - axis') 
    # naming the y axis 
    #plt.ylabel('y - axis')
    
    plt.subplot(212)    
    plt.plot(xnew, y_smooth)
    
    
#     y_wt = WaveletTransfer(y_smooth)    
#     plt.subplot(313)
#     plt.ylim(y_wt.min()-1, y_wt.max() + 1) 
#     plt.xlim(x.min()*.85, x.max()*1.1) 
#     print(len(xnew), len(y_wt))
#     plt.plot(xnew, y_wt)
    
    # giving a title to my graph 
    plt.subplot(211)
    plt.title(y) 
      
    # function to show the plot 
    #plt.show()
    if status == 0:
        plt.savefig(filePath)
    else:
        #plt.show();
        plt.savefig(filePath)  
    plt.clf()
    plt.cla()
    plt.close()

def ReadFile(filePath):
    ofs = open(filePath, "r");
    lines = ofs.readlines();
    features = [];
    for line in lines:
        features.append(line.split(","))
    nparryFeatures = np.array(features)
    nparryFeatures = nparryFeatures.astype(np.float)       
    return nparryFeatures
            
#Main 
# x axis values 
x1 = [1,2,3,4,5] 
x = np.array(x1)
print(x.shape)      

# corresponding y axis values 
yPositive = ReadFile("./Positive/PositiveSamplesDepthInfo.txt");

rootPath = "./Positive/"
DrawPics(x, yPositive[29], rootPath + str(0) + ".jpg", 0)
#print(yPositive[0].shape)
 
for i in range (0, yPositive.shape[0]):
    if i < yPositive.shape[0] - 1: 
        DrawPics(x, yPositive[i], rootPath + str(i) + ".jpg", 0)
        #DrawPics(x, yPositive[i], rootPath + "sum" + ".jpg", 0)
    else:
        DrawPics(x, yPositive[i], rootPath + str(i) + ".jpg", 0)
        #DrawPics(x, yPositive[i], rootPath + "sum" + ".jpg", 1)
  
#For Negative
yNegative = ReadFile("./Negative/NegativeSamplesDepthInfo.txt");
rootPath = "./Negative/"
print(yNegative[0].shape)
for i in range (0, yNegative.shape[0]):
    if i < yNegative.shape[0] - 1: 
        DrawPics(x, yNegative[i], rootPath + str(i) + ".jpg", 0)
        #DrawPics(x, yNegative[i], rootPath + "sum" + ".jpg", 0)
    else:
        DrawPics(x, yNegative[i], rootPath + str(i) + ".jpg", 1)
        #DrawPics(x, yNegative[i], rootPath + "sum" + ".jpg", 0)
      
