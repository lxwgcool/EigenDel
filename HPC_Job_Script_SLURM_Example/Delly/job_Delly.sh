#!/bin/bash

#SBATCH --ntasks=24
#SBATCH --nodes=1
#SBATCH --mail-type=END
#SBATCH --mail-user=lix33@nih.gov
#SBATCH --output=./log/std.out
#SBATCH --error=./log/std.err
#SBATCH --job-name=Delly
#SBATCH --time=4-00:00:00

module load delly/0.8.7
module load bcftools/1.13

#Run CNVnator by using different 1000 genome samples
bash ./script_Delly.sh