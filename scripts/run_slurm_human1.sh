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
dat=$(date +%Y%m%d%H%M%S)

bamdir="."
bamfiles=$(ls ${bamdir} | grep '.bam$' | xargs -I {} echo ${bamdir}/{})
bamfiles_="--bamFile=${bamfiles}"
bamfiles1=$(echo $bamfiles_ | sed 's/ /:/g')

##SPECIFY LOCATION OF COMBINED AND VIRUS ONLY DB
#reference="../Chlorocebus_sabaeus.ChlSab1.1.dna.toplevel.fa.gz"
#coord_file="../Chlorocebus_sabaeus.ChlSab1.1.99.gff3.gz"
reference="./human_virus_sequin_ensembl_pri_merged_genome.fasta"
coord_file="./gencode.v28.annotation.gff3.gz"
reference_virus="${npTranscript}/data/SARS-Cov2/VIC01/wuhan_coronavirus_australia.fasta.gz"
coord_file_virus="${npTranscript}/data/SARS-Cov2/VIC01/Coordinates.csv"
cov_chr=$(zcat ${coord_file_virus} | head -n 1 | cut -f 1 -d ' ' | sed 's/>//g')
echo "coronavirus chr id ${cov_chr}" 
resdir="results_${dat}"
opts="--bin=100 --breakThresh=100 --coronavirus=false --maxThreads=8 --extra_threshold=1000 --writePolyA=true --msaDepthThresh=1000 --doMSA=false --numExonsMSA=1:2:3:4:5 --msa_source=RNA --useExons=true --span=protein_coding --includeStart=false --isoformDepthThresh 1000"

#for dRNA datasets
opts="${opts} --RNA=true"
opts2="--fail_thresh=0  --chromsToRemap=${cov_chr}  --mm2_memory=10g --recordDepthByPosition=true"
opts3="--GFF_features=gene_id:description:gene_name:gene_type:gene_name"

bash ${npTranscript}/scripts/run.sh ${bamfiles1}   --reference=${reference} --annotation ${coord_file} --resdir ${resdir} ${opts} ${opts1} ${opts2} ${opts3}



cd ${resdir}
Rscript ${npTranscript}/R/npDE.R  control infected betabinom
