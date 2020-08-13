#!/bin/bash

#SBATCH --job-name=npTranscript
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=31800 # mb
#SBATCH --time=100:00:00
#SBATCH --output=corona1.stdout
#SBATCH --error=corona1.stderr
#SBATCH --cpus-per-task=8


##script for running 
##tip - use symbolic link to put this in the directory with bam files
export JSA_MEM=30000m
npTranscript=${HOME}/github/npTranscript

#bamfiles=bins/allreads_fq.bam:virion/sorted.virion_refmap.bam
bamdir="."
bamfiles=$(ls ${bamdir} | grep '.bam$'  | grep 'primary' | grep 'sorted' | grep -v 'DNA'  |  xargs -I {} echo ${bamdir}/{})
bamfiles1=$(echo $bamfiles | sed 's/ /:/g')
echo $bamfiles1
opts1="" 

#following if you want to only analysis a subset of reads
a=$(ls  | grep '^reads_in.txt$'  | wc -l)
if [ $a -eq 1 ]; then 
  opts1="--readList reads_in.txt" 
fi

#reference="../Chlorocebus_sabaeus.ChlSab1.1.dna.toplevel.fa.gz"
#coord_file="../Chlorocebus_sabaeus.ChlSab1.1.99.gff3.gz"
export coord_file="gencode.v28.annotation.gff3"
export coord_file="/DataOnline/Data/Projects/Human_Direct_RNA_Sequencing/Analysis_Guppy_Data/Alignment_minimap2/gencode.v28.annotation.gff3" 
export coord_file=$(pwd)"/gencode.v28.annotation.gff3.gz"
export reference=$(pwd)"/GRCh38_reference_genome/GRCh38_chromosome_only.fa.gz" 
export reference=$(pwd)"/GRCh38_full_analysis_set_plus_decoy_hla.fa"
#reference="/DataOnline/Data/H.sapiens/Reference/GRCh38_reference_genome/GRCh38_chromosome_only.fa"

dat=$(date +%Y%m%d%H%M%S)
resdir="results_${dat}"
opts="--bin=100 --breakThresh=100 --coronavirus=false --maxThreads=8 --extra_threshold=200 --writePolyA=true --msaDepthThresh=50 --doMSA=false --numExonsMSA=1:2:3:4:5 --msa_source=RNA;cDNA  --GFF_features=ID:description:ID:gene_type:gene_name --useExons=true --span=protein_coding --includeStart=false --isoformDepthThresh 100:100:100:100:1:100"
#for dRNA datasets
opts="${opts} --RNA=true"

#opts="${opts} --maxReads 10000"
bash ${npTranscript}/scripts/run.sh --bamFile=${bamfiles1}   --reference=${reference} --annotation ${coord_file} --resdir ${resdir} ${opts} ${opts1}


cd ${resdir}
#Rscript ~/github/npTranscript/R/npDE.R  control infected betabinom
