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

#bamfiles=bins/allreads_fq.bam:virion/sorted.virion_refmap.bam
bamdir="."
bamfiles=$(ls ${bamdir} | grep '.bam$' | xargs -I {} echo ${bamdir}/{})
bamfiles1=$(echo $bamfiles | sed 's/ /:/g')
echo $bamfiles1
opts1="" 
a=$(ls  | grep '^reads_in.txt$'  | wc -l)
if [ $a -eq 1 ]; then 
  opts1="--readList reads_in.txt" 
fi

reference="${npTranscript}/data/SARS-Cov2/wuhan_coronavirus_australia.fasta.gz"
coord_file="${npTranscript}/data/SARS-Cov2/Coordinates.csv"

dat=$(date +%Y%m%d%H%M%S)
resdir="results_${dat}"
opts="--bin 10 --breakThresh 1000 --cluster_by_annotation true"
#opts="${opts} --maxReads 10000"
bash ${npTranscript}/scripts/run.sh --bamFile=${bamfiles1}   --reference=${reference} --annotation ${coord_file} --resdir ${resdir} ${opts} ${opts1}
bash ${npTranscript}/scripts/run_extract_cluster.sh --inDir ${resdir} 

cd ${resdir}
Rscript ~/github/npTranscript/R/npTranscript.R  ~/github/npTranscript/data/SARS-Cov2
