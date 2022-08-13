#!/bin/bash

#SBATCH --ntasks=24
#SBATCH --nodes=1
#SBATCH --mail-type=END
#SBATCH --mail-user=lix33@nih.gov
#SBATCH --output=./log/std.out
#SBATCH --error=./log/std.err
#SBATCH --job-name=GASVPro
#SBATCH --time=4-00:00:00

module load java/12.0.1

#Run CNVnator by using different 1000 genome samples
bash ./script_GASVPro.sh