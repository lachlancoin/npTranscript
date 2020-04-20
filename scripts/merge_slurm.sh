#!/bin/bash

#SBATCH --job-name=corona2_analysis
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=9000 # mb
#SBATCH --time=100:00:00
#SBATCH --output=corona1.stdout
#SBATCH --error=corona1.stderr
#SBATCH --cpus-per-task=8


##script for running 
##tip - use symbolic link to put this in the directory with bam files
export JSA_MEM=8000m
npTranscript=${HOME}/github/npTranscript

##extra option
#--tocombine  results_host_infected1:results_host_control1"
#--filter \t5_3\t
#--ignore xxx

dat=$(date +%Y%m%d%H%M%S)
resdir="merged_${dat}"
opts="--inDir ."

##NOTE in dir contains the directories to be merged

bash ${npTranscript}/scripts/merge.sh --resdir ${resdir} ${opts} 


#can then run R scripts
#cd ${resdir}
#Rscript ~/github/npTranscript/R/npTranscript.R  ~/github/npTranscript/data/SARS-Cov2
