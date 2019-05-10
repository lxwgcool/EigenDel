***********************************************************
This is the introduction of how to compile and use EigenDel
***********************************************************

1: EigenDel use two tools below to calculate three global values base on bam file:
(1) Picard is used to calculate average insert size and the standard deviation of insert size.  
(2) Samtools is used to calculate the average depth (without include the range without coverage)

2: Dependent libraries
(1) Samtools v1.3
(2) Picard v2.0.1
(3) gcc v5.4.0 or higher 
(4) python 3 with scikit-learn be installed 

Notice: 
(1)Picard v2.0.1 requires Java 1.8 (jdk8u66) and R.
(2) Use UConn HPC as an example, the modules below are loaded:
module load samtools/1.3 gcc/5.4.0-alt java/1.8.0_31 sqlite/3.18.0 tcl/8.6.6.8606 python/3.6.1 r/3.1.1


3: How to compile EigenDel 
$ git clone https://github.com/lxwgcool/EigenDel.git
$ cd ./EigenDel/SourceCode/EigenDel/MakeFile
$ make clean
$ make

4: How to use EigenDel
(1) Set configure file: "config.ini"
[Debug]
VCF=/scratch/xil14026/sv/data/1000Genome/phase3_sv/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf
BamFile=/scratch/xil14026/sv/data/1000Genome/phase3_sv/Samples/NA12763/whole_genome_bam/NA12763.mapped.ILLUMINA.bwa.CEU.low_coverage.20130502.bam
Ref=/scratch/xil14026/sv/data/1000Genome/phase3_sv/ref/hs37d5.fa
PythonCode=../../../../Learning/SvLearning/SvClutering.py
Picard=/home/hpc-xin/software/picard/picard/build/libs/picard.jar
SampleName=NA12763
ChromIndex=
MultiAlignSeq=
PindelPath=
LumpyPath=
DellyPath=
GASVproPath=
CNVnatorPath=

*The meaning of each configure option
a) VCF: 
The SV annotation file provided by benchmark dataset, such as 1000 genomes project. This file is only for comparsion and generate scatter plot.
b) BamFile:
Reads alignment result.
c) Ref: 
Reference file used in BAM
d) PythonCode:
Python code used by supervised learning
e) Picard
The parth of picard
f) SampleName
The sample name you want to analysis
g) ChromIndex
The chromosome you want to analysis. Leave it be blank for analyzing the whole genome.
h) MultiAlignSeq -> discarded
i) PindelPath, LumpyPath, DellyPath, GASVproPath, CNVnatorPath
The path of each tool for comparison.

(2) Run EigenDel
$ EigenDel ./config.ini

(3) Results:
The result is named as "Result_Del_Sum.txt", and it is Located in the same folder as config.ini.
<<<<<<< HEAD
=======




>>>>>>> 46e8103808c5bb4d3769ac47d39781d0e163054f
