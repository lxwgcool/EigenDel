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
   echo "Step 1: Call Delly to process BAM File"
   bcfFile="${strCurSampleDir}/delly.sv.callset.bcf"
   CMD="delly call -o ${bcfFile} -g ${strRef} ${strBAM}"
   echo ${CMD}
   #eval ${CMD}
   echo

   #2: Use bcf view to check results
   vcfFile="${strCurSampleDir}/delly.sv.callset.vcf"
   CMD="bcftools view ${bcfFile} > ${vcfFile}"
   echo ${CMD}
   #eval ${CMD}
   echo

   #2: Get True Set (The quality equals to PASS)
   echo "Call structure variation"
   vcfPassFile="${strCurSampleDir}/delly.sv.callset.pass.vcf"
   awk '{if($7 == "PASS") print $0}' ${vcfFile} > ${vcfPassFile}

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