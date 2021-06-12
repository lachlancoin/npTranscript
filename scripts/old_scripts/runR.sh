#!/bin/bash

#SBATCH --job-name=corona2_R
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=9000 # mb
#SBATCH --time=100:00:00
#SBATCH --output=corona1.stdout
#SBATCH --error=corona1.stderr
#SBATCH --cpus-per-task=8
#SBATCH --partition=physical,snowy

# USAGE:
# sbatch --job-name=R  runR.sh

## Load required modules
module load gcc/8.3.0 
module load openmpi/3.1.4
module load r-bundle-bioconductor/3.9-r-3.6.0

if [ ! $npTranscript ] ; then
#export npTranscript=${HOME}/github/npTranscript
export npTranscript=/data/gpfs/projects/punim1068/npTranscript
fi 

libdir=/data/gpfs/projects/punim1068/Rlib
datasrc=${npTranscript}/data/SARS-Cov2/VIC01
Rscript ${npTranscript}/R/npTranscript.R  np.libdir=${libdir} np.source=${npTranscript}/R np.datasource=${datasrc}  $@ 

