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


pwd=$(pwd)
cd ~/github/npTranscript
git pull
cd $pwd
Rscript ~/github/npTranscript/R/npTranscript.R 

