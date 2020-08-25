#!/bin/bash

#SBATCH --job-name=corona2_analysis
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=15000 # mb
#SBATCH --time=200:00:00
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
bamfiles_="--bamFile=${bamfiles}"
bamfiles1=$(echo $bamfiles_ | sed 's/ /:/g')
echo $bamfiles1
opts1="" 
a=$(ls  | grep '^reads_in.txt$'  | wc -l)
if [ $a -eq 1 ]; then 
  opts1="--readList reads_in.txt" 
fi

#reference="${npTranscript}/data/SARS-Cov2/NC/NC_045512.fasta.gz" 
#coord_file="${npTranscript}/data/SARS-Cov2/NC/Coordinates.csv"

reference="${npTranscript}/data/SARS-Cov2/VIC01/wuhan_coronavirus_australia.fasta.gz"
coord_file="${npTranscript}/data/SARS-Cov2/VIC01/Coordinates.csv"
#alias abpoa='/sw/abpoa/v1.0.1/abpoa'
dat=$(date +%Y%m%d%H%M%S)
echo $dat
resdir="results_${dat}"
##note that to do separate msa for multiple input bams you type --doMsa=5_3:sep
opts="--bin=100 --RNA true --breakThresh=1000  --isoformDepthThresh=10000 --coverageDepthThresh=0 --extra_threshold=200 --msaDepthThresh=20 --doMSA=all:sep --reAlignExtra=true --bedChr=NC_045512v2"
opts2="--fail_thresh=7 --recordDepthByPosition=true" 
#opts="${opts} --maxReads 10000"
bash ${npTranscript}/scripts/run.sh ${bamfiles1}   --reference=${reference} --annotation ${coord_file} --resdir ${resdir} ${opts} ${opts1} ${opts2}

#cd ${resdir}
#bash ${npTranscript}/scripts/consensus.sh 0
#echo "now running R script"
#Rscript ~/github/npTranscript/R/npTranscript.R  ~/github/npTranscript/data/SARS-Cov2/VIC01
