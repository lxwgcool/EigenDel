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

   #1: Extract Fastq from BAM
   echo "1: Extract Fastq from BAM"
   bam2fastq="/data/COVID_ADHOC/Sequencing/Debug/lxwg/SV/Software/bam2fastq/bam2fastq"
   dirBAM=$(dirname ${bamFile})
   readsPE1="${dirBAM}/${strSampleID}.reads.PE_1.fq"
   readsPE2="${dirBAM}/${strSampleID}.reads.PE_2.fq"
   #CMD="bedtools bamtofastq -i ${strBAM} -fq ${readsPE1} -fq2 ${readsPE2}"
   CMD="${bam2fastq} ${strBAM} --output ${dirBAM}/${strSampleID}.reads.PE#.fq"
   echo ${CMD}
   eval ${CMD}

   #2: generate bam File (mapped, split bam and discordant)
   echo "2: generate bam File (mapped, split bam and discordant)"
   cd ${strCurSampleDir} && speedseq align -R "@RG\tID:id\tSM:sample\tLB:lib" -t 24 ${strRef} ${readsPE1} ${readsPE2}

   #3:run lumpy
   echo "3: run lumpy"
   SPLITBAM=${strCurSampleDir}/${strSampleID}.reads.PE_1.fq.splitters.bam
   DISCORDBAM=${strCurSampleDir}/${strSampleID}.reads.PE_1.fq.discordants.bam
   vcfFile=${strCurSampleDir}/${strSampleID}.sv.callset.vcf
   lumpyexpress -B ${strBAM} -S $SPLITBAM -D $DISCORDBAM -o ${vcfFile}

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
   echo
done