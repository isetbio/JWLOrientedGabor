#!/bin/bash
#SBATCH -p general #partition
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 0-20:00 # time (D-HH:MM)
#SBATCH --mem 128 # memory pool for all cores
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ek99@nyu.edu
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err

module load matlab/2016b

cd /scratch/ek99/JWLOrientedGabor/

# If the files you are running are not in the same folder as this script,
# you can insert "addpath(genpath('/PATH/TO/FILES/'));" before the command
# you want to run.
matlab -nodisplay -r "addpath(genpath('/scratch/ek99/JWLOrientedGabor')); addpath(genpath('/scratch/ek99/isetbio')); f_ogRGC_Classify([],0,{'100'},0); exit()"

exit

