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

   #1: Create BAM Config
   BAMCONFIG=${strCurSampleDir}/bam.config
   echo "${strBAM} 400 1000Genome" > ${BAMCONFIG}

   #1: Call Pindel
   echo "2: Call Pindel"
   CHROM=ALL
   CMD="pindel -f ${strRef} -i ${BAMCONFIG} -c ${CHROM} -o ${strCurSampleDir}/${strSampleID} -T 24"
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