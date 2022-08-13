#!/bin/bash

#SBATCH --ntasks=24
#SBATCH --nodes=1
#SBATCH --mail-type=END
#SBATCH --mail-user=lix33@nih.gov
#SBATCH --output=./log/std.out
#SBATCH --error=./log/std.err
#SBATCH --job-name=SVABA
#SBATCH --time=4-00:00:00

module load samtools/1.14
module load bwa-mem2/2.2.1

#Run CNVnator by using different 1000 genome samples
bash ./script_SVABA.sh
