#!/bin/bash

#SBATCH --job-name=npTranscript
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

#following if you want to only analysis a subset of reads
a=$(ls  | grep '^reads_in.txt$'  | wc -l)
if [ $a -eq 1 ]; then 
  opts1="--readList reads_in.txt" 
fi

reference="../Chlorocebus_sabaeus.ChlSab1.1.dna.toplevel.fa.gz"
coord_file="../Chlorocebus_sabaeus.ChlSab1.1.99.gff3.gz"

dat=$(date +%Y%m%d%H%M%S)
resdir="results_${dat}"
opts="--bin=50 --breakThresh=100 --coronavirus=false --extra_threshold=500 --msaDepthThresh=50 --doMSA=false"
#opts="${opts} --maxReads 10000"
#for dRNA datasets
opts="${opts} --RNA=true"
bash ${npTranscript}/scripts/run.sh --bamFile=${bamfiles1}   --reference=${reference} --annotation ${coord_file} --resdir ${resdir} ${opts} ${opts1}


cd ${resdir}
Rscript ~/github/npTranscript/R/npDE.R  control infected betabinom
