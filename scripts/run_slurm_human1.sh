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
bamc=$(ls  | grep '.bam' | wc  -l | cut -f 1)
if [ $bamc -gt 0 ]; then 
	bamfiles=$(ls ${bamdir} | grep '.bam$' | xargs -I {} echo ${bamdir}/{})
	bamfiles_="--bamFile=${bamfiles}"
else
        bamfiles=$(ls ${bamdir} | grep '.fastq$' | xargs -I {} echo ${bamdir}/{})
        bamfiles_="--fastqFile=${bamfiles}"
fi
bamfiles1=$(echo $bamfiles_ | sed 's/ /:/g')
echo $bamfiles1
opts1="" 

#following if you want to only analysis a subset of reads
a=$(ls  | grep '^reads_in.txt$'  | wc -l)
if [ $a -eq 1 ]; then 
  opts1="--readList reads_in.txt" 
fi

#reference="../Chlorocebus_sabaeus.ChlSab1.1.dna.toplevel.fa.gz"
#coord_file="../Chlorocebus_sabaeus.ChlSab1.1.99.gff3.gz"
reference="./human_virus_sequin_ensembl_pri_merged_genome.fasta"
coord_file="./gencode.v28.annotation.gff3.gz"
dat=$(date +%Y%m%d%H%M%S)
resdir="results_${dat}"
opts="--bin=100 --breakThresh=100 --coronavirus=false --maxThreads=8 --extra_threshold=1000 --writePolyA=true --msaDepthThresh=1000 --doMSA=false --numExonsMSA=1:2:3:4:5 --msa_source=RNA --useExons=true --span=protein_coding --includeStart=false --isoformDepthThresh 1000"
#for dRNA datasets
opts="${opts} --RNA=true"
opts2="--fail_thresh=0  --chromsToRemap=MT007544.1  --mm2_memory=10g --recordDepthByPosition=true"
opts3="--GFF_features=ID:description:ID:gene_type:gene_name"

#opts="${opts} --maxReads 10000"
bash ${npTranscript}/scripts/run.sh ${bamfiles1}   --reference=${reference} --annotation ${coord_file} --resdir ${resdir} ${opts} ${opts1} ${opts2} ${opts3}


cd ${resdir}
#Rscript ~/github/npTranscript/R/npDE.R  control infected betabinom
