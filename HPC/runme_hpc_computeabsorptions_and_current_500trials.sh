#!/bin/bash
#SBATCH --nodes=1
#SBATCH --array=1-17
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time 24:00:00 # time (D-HH:MM)
#SBATCH --mem=64GB # memory pool for all cores
#SBATCH --job-name=ISETBIO_conedensity_LOnly_500trials
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ek99@nyu.edu
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err

module load matlab/2016b

cd /scratch/ek99/JWLOrientedGabor/

# If the files you are running are not in the same folder as this script,
# you can insert "addpath(genpath('/PATH/TO/FILES/'));" before the command
# you want to run.
matlab -nodisplay -r "addpath(genpath('/scratch/ek99/JWLOrientedGabor')); addpath(genpath('/scratch/ek99/isetbio')); runComputationalObserverModel('conedensitynophaseshiftlonly500', 'saveFolder', 'run2', 'seed', 1); exit()"

exit

