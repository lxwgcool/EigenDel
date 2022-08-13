#!/bin/bash

strRef="/data/COVID_ADHOC/Sequencing/Debug/lxwg/Data/1000Genome/Ref/hs37d5.fa"
strRefDir="/data/COVID_ADHOC/Sequencing/Debug/lxwg/Data/1000Genome/Ref"

function func_callSV()
{
   echo "Sample ID    : $1"
   echo "CurSample Dir: $2"
   echo "Bam File     : $3"
   echo
   strSampleID=$1
   strCurSampleDir=$2
   strBAM=$3

   BINSIZE=500

   # 1: Extract read mapping
   echo "Step 1: Extract read mapping"
   CMD="cnvnator -root file.root -tree ${strBAM}"
   echo ${CMD}
   eval ${CMD}
   echo
   
   # 2: Generate histogram
   echo "Step 2: Generate histogram"
   CMD="cnvnator -root file.root -his $BINSIZE -chrom 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"
   echo ${CMD}
   eval ${CMD}
   echo

   # 3: Calculate statistics
   echo "Step 3: Calculate statistics"
   CMD="cnvnator -root file.root -stat $BINSIZE -d ${strRefDir}/"
   echo ${CMD}
   eval ${CMD}
   echo

   # 4: Partition
   echo "Step 4: Partition"
   CMD="cnvnator -root file.root -partition $BINSIZE"
   echo ${CMD}
   eval ${CMD}
   echo

   # 5:Call CNVs
   echo "Step 5: Call CNVs"
   CMD="cnvnator -root file.root -call $BINSIZE > ${strCurSampleDir}/call_sv_cnvnator_${strSampleID}.txt"
   echo ${CMD}
   eval ${CMD}
   echo
   
   echo "============="
   echo
}

echo "==============="
date
echo "==============="
echo

strRoot="/data/COVID_ADHOC/Sequencing/Debug/lxwg/Data/1000Genome"
strBaseDir=$(dirname $(realpath $0))

for bamFile in $(find ${strRoot} -iname '*.bam')
do
   sampleID=$(basename $(dirname ${bamFile}))
   strCurSampleDir=${strBaseDir}/${sampleID}

   # Create Folder 
   mkdir -p ${strCurSampleDir}   

   # Call
   SECONDS=0 	
   
   cd ${strCurSampleDir} && func_callSV "${sampleID}" "${strCurSampleDir}" "${bamFile}"
   
   result=$?
   duration=$SECONDS
   echo "Running Time: $(($duration / 3600))hrs $((($duration / 60) % 60))min $(($duration % 60))sec"
   if [[ ${result} -ne 0 ]]; then
      echo "Error: $(date) run was failed!"
      exit 1
   else
      echo "$(date) run was finished successfully!"
   fi
done
