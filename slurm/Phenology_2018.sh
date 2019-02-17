#!/bin/bash
#SBATCH --job-name=PhenologyModel-%A-%a  #
#SBATCH --output=RunTimeOut/%A-%a.out  # Standard out goes to this file
#SBATCH --error=RunTimeOut/%A-%a.err   # Standard err goes to this file
#SBATCH -n 1 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 600 # Runtime in minutes
#SBATCH -p wolkovich # Partition to submit to (mine!)
#SBATCH --mem=1000 # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH --mail-type=END
#SBATCH --mail-user=donahuem@hawaii.edu
mkdir -p /scratch/wolkovich_lab/temporalvar/$SLURM_JOB_ID
module load R
export R_LIBS_USER=$HOME/apps/R:$R_LIBS_USER
Rscript /n/wolkovich_lab/temporalvar/R/PhenologyModel.R --quiet --no-restore --no-save
rm -rf /scratch/wolkovich_lab/temporalvar/$SLURM_JOB_ID
###run this with the following command from the temporalvar folder
###sbatch --array=1-10 slurm/Phenology_20160331.sh
###where 10 is the number of jobs in the array
###where Phenology_20130331.sh is the name of the batch script
