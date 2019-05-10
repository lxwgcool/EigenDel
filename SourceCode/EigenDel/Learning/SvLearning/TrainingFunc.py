from sklearn import svm
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import make_classification

from sklearn.tree import DecisionTreeRegressor
from sklearn.tree import DecisionTreeClassifier

from sklearn.ensemble import AdaBoostClassifier
from sklearn.ensemble import AdaBoostRegressor

def TrainingSVM(arrayPos, arrayNeg, arrayPosTesting, arrayNegTesting): #  Do tomorrow
    #1: Prepare data
    numPosSample, numPosFeature = arrayPos.shape
    
    numNegSample, numNegFeature = arrayNeg.shape
    print(numPosSample, numNegSample)
        
    arryLabel = [];
    for i in range(0, numPosSample): # add label for positive samples
        arryLabel.append(1)
    for i in range(0, numNegSample): # add label for negative samples
        arryLabel.append(0)
    print(len(arryLabel))
        
    #Get Sample sum 
    arraySampleSum = np.concatenate((arrayPos, arrayNeg), axis=0)    
     
    #2: SVM Training --> Now let's for training --> Go!!
    #clf = svm.SVC(gamma='scale')
    
#     clf = RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1)

    clf = AdaBoostClassifier(DecisionTreeClassifier(max_depth=1),
                         algorithm="SAMME", n_estimators=200)
    clf.fit(arraySampleSum, arryLabel)     
    
    print(clf.predict(arrayPosTesting))    
    return;