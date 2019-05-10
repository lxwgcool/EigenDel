import pandas as pd 
import numpy as np
import sys
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA, KernelPCA
from mpl_toolkits.mplot3d import Axes3D
from sklearn.manifold import TSNE
from sklearn import datasets
from sklearn.cluster import DBSCAN
import scipy.cluster.hierarchy as shc
from sklearn.cluster import AgglomerativeClustering

def Draw2DPics(samples, labels, strChromName, strSampleName):    
    x = []
    y = []
    n = []    
    for i in range(0, len(samples)):
        x.append(samples[i][0])
        y.append(samples[i][1])
        n.append(i)          
                 
    npX = np.array(x)
    npY = np.array(y)    
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('PC_1', fontsize = 14)
    ax.set_ylabel('PC_2', fontsize = 14)
    ax.set_title('PCA ' + strSampleName + " " + strChromName, fontsize = 19)
    
    print(len(labels[0]))
    print(npX.shape, npY.shape)

#     for i in range(0, len(labels[0])):
#         if labels[0][i] == 1:
#             color = 'c'
#         elif labels[0][i] == 0:
#             color = 'k'
#         else:
#             color = 'b'
#         ax.scatter(npX[i], npY[i], c=color)

    #1: Draw 0 first
    for i in range(0, len(labels[0])):
        if labels[0][i] == 0:
            color = 'k'
            ax.scatter(npX[i], npY[i], c=color) 
         
    #2: Draw 1 second
    for i in range(0, len(labels[0])):
        if labels[0][i] == 1:
            color = 'c'
            ax.scatter(npX[i], npY[i], c=color)
        
#     for i, txt in enumerate(n):
#         ax.annotate(txt, (npX[i], npY[i]))
    plt.savefig("./Pics/PCA_" + strSampleName + "_Chrom_" + strChromName + ".png")                    
    #plt.show()    


def Draw3DPics(samples, labels, strChromName):
#    fig = plt.figure(figsize = (8,8))
    
    #ax = fig.add_subplot(1,1,1) 
    #ax.set_xlabel('Principal Component 1', fontsize = 15)
    #ax.set_ylabel('Principal Component 2', fontsize = 15)
    #ax.set_title('2 Component PCA', fontsize = 20)
    
    x = []
    y = []
    z = []
    for i in range(0, len(samples)):
        x.append(samples[i][0])
        y.append(samples[i][1])
        z.append(samples[i][2])
    
    print(x)
    print(y)
    print(z)             
    npX = np.array(x)
    npY = np.array(y)
    npZ = np.array(z)
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('PC_1', fontsize = 10)
    ax.set_ylabel('PC_2', fontsize = 10)
    ax.set_zlabel('PC_3', fontsize = 10)
    ax.set_title('3 Component PCA', fontsize = 15)
    
    for i in range(0, len(labels[0])):
        if labels[0][i] == 1:
            color = 'r'
        elif labels[0][i] == 0:
            color = 'g'
        else:
            color = 'b'
        ax.scatter(npX[i], npY[i], npZ[i], c=color)                
    plt.show()
        
#     targets = ['Iris-setosa', 'Iris-versicolor', 'Iris-virginica']
#     colors = ['r', 'g', 'b']
#     for target, color in zip(targets,colors):
#         indicesToKeep = finalDf['target'] == target
#         ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
#                    , finalDf.loc[indicesToKeep, 'principal component 2']
#                    , c = color
#                    , s = 50)
#     ax.legend(targets)
#     ax.grid()

def RunPCA3D(samples):    
    pca = PCA(n_components=3)
    #pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(samples)
    principalDf = pd.DataFrame(data = principalComponents
                 , columns = ['PC_1', 'PC_2', 'PC_3']) 
    #, columns = ['PC_1', 'PC_2'])
    print(principalDf.head(5))
    return principalComponents    

def RunPCA(samples):    
    #pca = PCA(n_components=3)
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(samples)
#     principalDf = pd.DataFrame(data = principalComponents
#     #             , columns = ['PC_1', 'PC_2', 'PC_3']) 
#     , columns = ['PC_1', 'PC_2'])
#     print(principalDf.head(5))
    return principalComponents

def TSN(samples, labels):
    model = TSNE(n_components=2, learning_rate=100)
    transformed = model.fit_transform(samples)
    
    x_axis = transformed[:, 0]
    y_axis = transformed[:, 1]
    #z_axis = transformed[:, 2]
    
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # ax.set_xlabel('PC_1', fontsize = 10)
    # ax.set_ylabel('PC_2', fontsize = 10)
    # ax.set_zlabel('PC_3', fontsize = 10)
    # ax.set_title('3 Component PCA', fontsize = 15)
    # 
    # for i in range(0, len(labels[0])):
    #     color = 'g'
    #     if i < 59:
    #         color = 'r'
    #     ax.scatter(x_axis[i], y_axis[i], z_axis[i], c=color)                
    # plt.show()
    
    plt.scatter(x_axis, y_axis, c=labels[0])
    plt.show()

def RunDBSCAN(samples, labels):
    clustering = DBSCAN(eps=0.3, min_samples=2).fit_predict(samples)
    #print(clustering.labels_)
    #print(clustering)
    plt.scatter(samples[:, 0], samples[:, 1], c=clustering)
    plt.show()
    
def DrawSHC(samples, labels):
    plt.title("Customer Dendograms")  
    dend = shc.dendrogram(shc.linkage(samples, method='ward'))
    plt.show()

def DrawAGGC(samples, labels, strChromName, strSampleName):
    cluster = AgglomerativeClustering(n_clusters=4, affinity='euclidean', linkage='ward')  
    cluster.fit_predict(samples)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title('HC ' + strSampleName + " " + strChromName, fontsize = 19)
                 
    #plt.scatter(samples[:,0], samples[:,1], c=cluster.labels_, cmap='rainbow')
    
    x = []
    y = []
    n = []    
    for i in range(0, len(samples)):
        x.append(samples[i][0])
        y.append(samples[i][1])        
        n.append(i)
    
    labels = cluster.labels_
#     print(labels)
#     print(len(labels))
    iBad = 0;
    for i in range(0, len(labels)):
        if labels[i] == 0:
            color = 'r'            
        elif labels[i] == 1:
            color = 'g'
            iBad = iBad + 1;
        elif labels[i] == 2:
            color = 'b'
        else:
            color = 'y'
            
        ax.scatter(x[i], y[i], c=color)  
                    
#     for i, txt in enumerate(n):
#         ax.annotate(txt, (x[i], y[i]))
     
    print("iBad: ", iBad)    
    #plt.show()
    plt.savefig("./Pics/AGGC_" + strSampleName + "_Chrom_" + strChromName + ".png")
    return labels
    
def ReadFile(filePath):
    ofs = open(filePath, "r");
    lines = ofs.readlines();
    features = [];
    for line in lines:
        features.append(line.split(","))
    nparryFeatures = np.array(features)
    nparryFeatures = nparryFeatures.astype(np.int)       
    return nparryFeatures

def ReadFileFloat(filePath):
    ofs = open(filePath, "r");
    lines = ofs.readlines();
    features = [];
    for line in lines:
        features.append(line.split(","))
    nparryFeatures = np.array(features)
    nparryFeatures = nparryFeatures.astype(np.float)       
    return nparryFeatures

def SaveLables(lables, strChromName):
    filePath = "./Pics/ClusterLables_Chrom_" + strChromName + ".txt"
    ifs = open(filePath, 'w')
    #Write the labels into this file
    for i in range(0, lables.shape[0]):
        ifs.write(str(lables[i]))        
    ifs.close()
    
#1: Print the input samples
#samples = ReadFile("/home/xin/lxwg/WorkStudio/Prototype/SV_Draft/SvDelPainter/Learning/SvLearning/Samples/NA12776/chrom11/CluterPureZero.txt");
#print(sys.argv[1])
print(sys.version)

# samples = ReadFile("./Samples/NA12776/chrom11/CluterPureZero.txt");
# labels = ReadFile("./Samples/NA12776/chrom11/CluterLabel.txt");
# strChromName = "ChromX"
# strSampleName = "NAxxxxx"

samples = ReadFile(sys.argv[1]) #ReadFile("./Samples/NA12776/chrom11/CluterPureZero.txt");
labels = ReadFile(sys.argv[2]) #ReadFile("./Samples/NA12776/chrom11/CluterLabel.txt");
strChromName = sys.argv[3]
strSampleName = sys.argv[4]

#print(labels[0])
print(labels[0].shape)
#features = ['No_Sv_Left', 'DC_Left', 'DC_Core', 'DC_Right', 'No_Sv_Right']
# features = ['DC_Left_Zero_Pure', 'DC_Left_Zero_Week', 'DC_Left_Zero_Mix', 'DC_Left_Rep_Week', 'DC_Left_Rep_Strong', 
#             'DC_Core_Zero_Pure', 'DC_Core_Zero_Week', 'DC_Core_Zero_Mix', 'DC_Core_Rep_Week', 'DC_Core_Rep_Strong',
#             'DC_Right_Zero_Pure', 'DC_Right_Zero_Week', 'DC_Right_Zero_Mix', 'DC_Right_Rep_Week', 'DC_Right_Rep_Strong']
# print(pd.DataFrame(data = samples, columns = features).head())

# principalComponents = RunPCA(samples)
#TSN(samples, labels)

#RunDBSCAN(samples, labels)
principalComponents = RunPCA(samples)
#DrawSHC(principalComponents, labels)

realClusterlables = DrawAGGC(principalComponents, labels, strChromName, strSampleName)
SaveLables(realClusterlables, strChromName)

#2: Do PCA and print the PCA result (the first 5 samples)
#3: Draw Pics
# kpca = KernelPCA(kernel="rbf", n_components = 2, fit_inverse_transform=True, gamma=10)
# X_kpca = kpca.fit_transform(samples)
# X_back = kpca.inverse_transform(X_kpca)

principalComponents = RunPCA(samples)
Draw2DPics(principalComponents, labels, strChromName, strSampleName)


# principalComponents = RunPCA3D(samples)
# Draw3DPics(principalComponents, labels)
