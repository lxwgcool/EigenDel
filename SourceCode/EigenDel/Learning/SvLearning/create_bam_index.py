import os

#1: Read all of files in current folder
#   Get the File List that need to create bai file (index file)    
def LoadFileInDir():
    arryFile = []
    strCmd = 'ls ./ > ./tmp.txt'
    os.system(strCmd)
    f = open('./tmp.txt', 'r')
    lines = f.readlines()        
    for line in lines:        
        if(line.find('REF') != -1):
            arryFile.append(line[0:len(line) - 1])
    
    arryRealFile = []
    for file in arryFile:
        if(file[(len(file) - 3) : len(file)] == "bam"):
            arryRealFile.append( "./" + file)                            
    return arryRealFile        
 
#2: user samtool to create index one by one
def CreateBamIndex(arryFile):
    for file in arryFile:
        strCmd = "samtools index " + file + " " + file + ".bai"
        print(strCmd)
        os.system(strCmd)
        
#### Main ####
arryFile = LoadFileInDir()
print(arryFile) 
CreateBamIndex(arryFile)

   