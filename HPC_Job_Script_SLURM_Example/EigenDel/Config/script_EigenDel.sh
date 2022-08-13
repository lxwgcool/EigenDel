#!/bin/bash

strRef="/data/COVID_ADHOC/Sequencing/Debug/lxwg/Data/1000Genome/Ref/hs37d5.fa"
strRefDir="/data/COVID_ADHOC/Sequencing/Debug/lxwg/Data/1000Genome/Ref"

#Run CNVnator by using different 1000 genome samples
EigenDel="/data/COVID_ADHOC/Sequencing/Debug/lxwg/SV/Callers/EigenDel/software/EigenDel/SourceCode/EigenDel/MakeFile/EigenDel"
dirOutput="/data/COVID_ADHOC/Sequencing/Debug/lxwg/SV/Callers/EigenDel/Output"
dirConfig="/data/COVID_ADHOC/Sequencing/Debug/lxwg/SV/Callers/EigenDel/Config"

function func_callSV()
{
   echo "Sample ID    : $1"
   echo "CurSample Dir: $2"
   echo "Bam File     : $3"
   echo

   strSampleID=$1
   strCurSampleDir=$2
   strBAM=$3

   # Step 1: Split BAM by chromosome
   echo "1: Split BAM by chromosome --> ${strSampleID}"
   #bamtools split -in ${strBAM} -reference

   # Step 2: Build Bai for each splited BAM
   # Get Main BAM Dir
   echo "2: Build Bai for each splited BAM"
   strBAMDir=$(dirname ${strBAM})
   echo "strBAMDir: ${strBAMDir}"
   for subBamFile in $(find ${strBAMDir} -iname '*REF_*.bam')
   do
     echo "Build Index for ${subBamFile}"
     #samtools index -@ 24 ${subBamFile} ${subBamFile}.bai
   done

   # Step 2: Run EigenDel
   echo "3: Run EigenDel --> ${strSampleID}"
   strConfig="${dirConfig}/Config_${strSampleID}.ini"
   ${EigenDel} ${strConfig}

   echo "============="
   echo "4: All Set Run EigenDel --> ${strSampleID}"
   echo
}

echo "==============="
date
echo "==============="
echo

strRoot="/data/COVID_ADHOC/Sequencing/Debug/lxwg/Data/1000Genome"
strBaseDir=$(dirname $(realpath $0))/../Output

#Build index  fastq
#samtools faidx ${strRef}
#bwa-mem2 index ${strRef}

for bamFile in $(find ${strRoot} -iname '*.bam')
do
   sampleID=$(basename $(dirname ${bamFile}))
   strCurSampleDir=${dirOutput}/${sampleID}

   # Create Folder
   mkdir -p ${strCurSampleDir}

   # Call
   SECONDS=0

   cd ${strCurSampleDir} && func_callSV "${sampleID}" "${strCurSampleDir}" "${bamFile}"

   result=$?
   duration=$SECONDS
   echo "Running Time: $(($duration / 3600))hrs $((($duration / 60) % 60))min $(($duration % 60))sec"
   if [[ ${result} -ne 0 ]]; then
      echo "Error: $(date) ${sampleID} run was failed!"
      exit 1
   else
      echo "$(date) ${sampleID} run was finished successfully!"
   fi
done
