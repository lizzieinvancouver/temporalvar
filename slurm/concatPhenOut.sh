#!/bin/bash                                                                                               

#SBATCH --job-name=concatPhenOut  #

#SBATCH --output=concatPhenOut-%A-%a.out  # Standard out goes to this file

#SBATCH --error=concatPhenOut-%A-%a.err		  # Standard err goes to this file                                          

#SBATCH -n 1 # Number of cores requested                                                                 

#SBATCH -N 1 # Ensure that all cores are on one machine                                                  

#SBATCH -t 10 # Runtime in minutes                                                                       

#SBATCH -p wolkovich # Partition to submit to (mine!)                                                    


#SBATCH --mem=10000 # Memory per cpu in MB (see also --mem-per-cpu)                                       

                                


#SBATCH --mail-type=END
 
#SBATCH --mail-user=donahuem@hawaii.edu<
 
source new-modules.sh
module load R_packages


R CMD BATCH --quiet --no-restore --no-save R/sourcefiles/concatPhenologyOut.R




###run this with the following command from the temporalvar folder



###sbatch slurm/concatPhenOut.sh

###where concatPhenOut.sh is the name of the batch script
