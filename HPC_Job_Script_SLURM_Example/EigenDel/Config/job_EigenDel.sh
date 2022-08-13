#!/bin/bash

#SBATCH --ntasks=24
#SBATCH --nodes=1
#SBATCH --mail-type=END
#SBATCH --mail-user=lix33@nih.gov
#SBATCH --output=../log/std.out
#SBATCH --error=../log/std.err
#SBATCH --job-name=EigenDel
#SBATCH --time=4-00:00:00

module load picard/2.26.9
module load python/3.7
module load gcc/7.3.0
module load R/3.4.3
module load samtools/1.14
module load bamtools/2.5.2

# Run EigenDel
bash ./script_EigenDel.sh