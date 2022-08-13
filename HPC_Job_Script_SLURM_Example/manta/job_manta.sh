#!/bin/bash

#SBATCH --ntasks=24
#SBATCH --nodes=1
#SBATCH --mail-type=END
#SBATCH --mail-user=lix33@nih.gov
#SBATCH --output=./log/std.out
#SBATCH --error=./log/std.err
#SBATCH --job-name=manta
#SBATCH --time=4-00:00:00

module load manta/1.6.0

#Run CNVnator by using different 1000 genome samples
bash ./script_manta.sh