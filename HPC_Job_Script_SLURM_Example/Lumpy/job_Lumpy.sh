#!/bin/bash

#SBATCH --ntasks=24
#SBATCH --nodes=1
#SBATCH --mail-type=END
#SBATCH --mail-user=lix33@nih.gov
#SBATCH --output=./log/std.out
#SBATCH --error=./log/std.err
#SBATCH --job-name=Lumpy
#SBATCH --time=4-00:00:00

module load speedseq/0.1.2-20180208-4e60002
module load lumpy/0.2.13
module load bedtools/2.30.0
module load samtools/1.14

#Run CNVnator by using different 1000 genome samples
bash ./script_Lumpy.sh