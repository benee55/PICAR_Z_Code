#!/bin/sh
#SBATCH --job-name=ch
#SBATCH --partition=normal
#SBATCH --output=/scratch/slee287/logFiles/%j.out  # Output file
#SBATCH --error=/scratch/slee287/logFiles/%j.err   # Error file
#SBATCH --time=0-05:00  # Total time needed for job: Days-Hours:Minutes
#SBATCH   --mem-per-cpu=10GB             # Request nGB RAM per core

module load r
cd /home/slee287/ZeroInflatedPICAR/code/realExample/bivalve/run
Rscript --no-restore --quiet --no-save  count_hurdle_PICAR.R

