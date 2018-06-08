#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time 20:00:00 # time (D-HH:MM)
#SBATCH --mem=125GB # memory pool for all cores
#SBATCH --job-name=ogRGC_classify_HPC
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ek99@nyu.edu
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err

module load matlab/2016b

cd /scratch/ek99/JWLOrientedGabor/

# If the files you are running are not in the same folder as this script,
# you can insert "addpath(genpath('/PATH/TO/FILES/'));" before the command
# you want to run.
matlab -nodisplay -r "addpath(genpath('/scratch/ek99/JWLOrientedGabor')); addpath(genpath('/scratch/ek99/isetbio')); expName = 'eccbasedcoverage'; subFolderName_toLoad = 'paddedStim'; subFolderName_toSave = '100trials_paddedStim'; ogRGC_Classify_HPC(expName, subFolderName_toLoad, subFolderName_toSave); exit()"

exit

