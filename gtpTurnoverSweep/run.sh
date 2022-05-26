#!/bin/bash
#SBATCH –J matlab
#SBATCH –o slurm.out
#SBATCH --mem=64G #	4	GB	RAM		
/opt/apps/matlabR2016a/bin/matlab -nojvm -nodisplay –singleCompThread -r gtpTurnoverSweep