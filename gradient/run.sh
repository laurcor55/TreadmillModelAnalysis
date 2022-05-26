#!/bin/bash
#SBATCH -p common
#SBATCH -J matlab
#SBATCH -o slurm_ec.out
#SBATCH --time=24:00:00
#SBATCH --mem=64G 
#SBATCH --ntasks=1	
module load Matlab/R2019a
matlab -nodisplay -singleCompThread -r gradient 