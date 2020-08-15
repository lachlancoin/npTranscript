#!/bin/bash

#SBATCH --job-name=DE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=9000 # mb
#SBATCH --time=100:00:00
#SBATCH --output=DE.stdout
#SBATCH --error=DE.stderr
#SBATCH --cpus-per-task=8


##USAGE bash runDE.sh control infected

#pwd=$(pwd)
#cd ~/github/npTranscript
#git pull
#cd $pwd
export R_LIBS_USER=/usr/local/lib/R/site-library/
Rscript --vanilla   ~/github/npTranscript/R/npDE.R  $1 $2 betabinom 

