#!/bin/bash                                                                                               

#SBATCH --job-name=Test  #

#SBATCH --output=Test-%A-%a.out  # Standard out goes to this file

#SBATCH --error=Test-%A-%a.err		  # Standard err goes to this file                                          

#SBATCH -n 1 # Number of cores requested                                                                 

#SBATCH -N 1 # Ensure that all cores are on one machine                                                  

#SBATCH -t 60 # Runtime in minutes                                                                       

#SBATCH -p wolkovich # Partition to submit to (mine!)                                                    

#SBATCH --mem=1000 # Memory per cpu in MB (see also --mem-per-cpu)                                       

#SBATCH --mail-type=END
 
#SBATCH --mail-user=donahuem@hawaii.edu<
 
source new-modules.sh
module load R_packages

R CMD BATCH --quiet --no-restore --no-save R/PhenologyModel.r

###run this with the following command from the temporalvar folder

###sbatch --array=1-5 slurm/testbatch.sh

###where 5 is the number of jobs in the array
###where testbatch.sh is the name of the batch script
