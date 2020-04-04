#!/bin/bash

#SBATCH --job-name=corona2_analysis
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32000 # mb
#SBATCH --time=100:00:00
#SBATCH --output=corona1.stdout
#SBATCH --error=corona1.stderr
#SBATCH --cpus-per-task=8


##script for running 
##tip - use symbolic link to put this in the directory with bam files

npTranscript=${HOME}/github/npTranscript

#bamfiles=bins/allreads_fq.bam:virion/sorted.virion_refmap.bam
bamdir="."
bamfiles=$(ls ${bamdir} | grep '.bam' | xargs -I {} echo ${bamdir}/{})



reference="${npTranscript}/data/SARS-Cov2/wuhan_coronavirus_australia.fasta.gz"
coord_file="${npTranscript}/data/SARS-Cov2/Coordinates.csv"

dat=$(date +%Y%m%d%H%M%S)
resdir="results_${dat}"
opts="--bin 10 --breakThresh 1000 --cluster_by_annotation true"
bash ${npTranscript}/run.sh --bamFile=${bamfiles}   --reference=${reference} --annotation ${coord_file} --resdir ${resdir} ${opts}

bash ${npTranscript}/run_extract_cluster.sh --inDir ${resdir} 
