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

   GASVPATH="/gpfs/gsfs6/users/COVID_ADHOC/Sequencing/Debug/lxwg/SV/Callers/GASVpro/bin"

   #1: Using BAMToGASV to Preprocess BAM les
   echo "Step 1: Using BAMToGASV to Preprocess BAM files"
   CMD="java -Xms512m -Xmx2048m -jar $GASVPATH/BAMToGASV.jar ${strBAM} -LIBRARY_SEPARATED all -OUTPUT_PREFIX ${strCurSampleDir}/${strSampleID}.bam"
   echo ${CMD}
   eval ${CMD}

   #2: Call structure variation
   echo "Call structure variation"
   CMD="java -jar $GASVPATH/GASV.jar --batch ${strCurSampleDir}/${strSampleID}.bam.gasv.in --outputdir ${strCurSampleDir}"
   echo ${CMD}
   eval ${CMD}

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