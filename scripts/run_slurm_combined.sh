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
mm2_path="/sw/minimap2/current/minimap2"

mode=$1
shift
if [ ! $mode ];then
 mode="combined"
fi

bamdir="."
bamfiles=$(ls ${bamdir} | grep '.bam$' | xargs -I {} echo ${bamdir}/{})
bamfiles_="--bamFile=${bamfiles}"
bamfiles1=$(echo $bamfiles_ | sed 's/ /:/g')

##SPECIFY LOCATION OF COMBINED AND VIRUS ONLY DB
#reference="../Chlorocebus_sabaeus.ChlSab1.1.dna.toplevel.fa.gz"
#coord_file="../Chlorocebus_sabaeus.ChlSab1.1.99.gff3.gz"
reference="/DataOnline/Data/Projects/corona_invitro/host_analysis/db/merged/human_virus_sequin_ensembl_pri_merged_genome.fasta"
coord_file="../gencode.v28.annotation.gff3.gz"
##FOLLOWING IS FEATURES IN GFF FILE.THIS MAY NEED TO BE CUSTOMISED
##IF uSING USE_eXONS=true, THEY SHOULD REFER TO GENE FEATURES, NOT EXON FEATURES
GFF_features="--GFF_features=gene_name:description:ID:gene_type:Parent"

reference_virus="${npTranscript}/data/SARS-Cov2/VIC01/wuhan_coronavirus_australia.fasta.gz"
coord_file_virus="${npTranscript}/data/SARS-Cov2/VIC01/Coordinates.csv"
cov_chr=$(zcat ${reference_virus} | head -n 1 | cut -f 1 -d ' ' | sed 's/>//g')
echo "coronavirus chr id ${cov_chr}" 
resdir="results_${dat}"
opts="--bin=100 --breakThresh=100 --coronavirus=false --maxThreads=8 --extra_threshold=100 --writePolyA=true --msaDepthThresh=1000 --doMSA=false --numExonsMSA=1:2:3:4:5 --msa_source=RNA --useExons=true --span=protein_coding --includeStart=false --isoformDepthThresh 50"

#for dRNA datasets
opts="${opts} --RNA=true"
opts2="--fail_thresh=0  --chromsToRemap=${cov_chr}  --mm2_memory=10g --recordDepthByPosition=true"
opts="${opts} $@"
echo $opts
bash ${npTranscript}/scripts/run.sh ${bamfiles1}   --reference=${reference} --annotation=${coord_file} --resdir=${resdir} ${opts} ${opts1} ${opts2} ${opts3} ${GFF_features}



cd ${resdir}


if [ $mode == "combined" ]; then
##VIRAL ANALYSIS ON VIRAL READS
bamdir="."
bamfiles=$(wc -l *fastq | tr -s ' ' | grep -v ' 0 '  | cut -f 3 -d ' ' | grep -v 'polyA'  | grep -v 'leftover' | grep -v 'total' |  xargs -I {} echo ${bamdir}/{})
bamfiles_="--fastqFile=${bamfiles}"
bamfiles_virus=$(echo $bamfiles_ | sed 's/ /:/g')


resdir_virus="results_virus"
resdir_leftover="results_leftover"

opts="--bin=100 --breakThresh=1000  --isoformDepthThresh=10000 --coverageDepthThresh=0 --extra_threshold=200 --msaDepthThresh=20 --doMSA=all:sep --reAlignExtra=true"
opts2="--fail_thresh=0 --recordDepthByPosition=true"
opts3="--mm2_path=${mm2_path}"
bash ${npTranscript}/scripts/run.sh ${bamfiles_virus}   --reference=${reference_virus} --annotation ${coord_file_virus} --resdir ${resdir_virus} ${opts} ${opts1} ${opts2} ${opts3}


bamdir="."
bamfiles=$(wc -l *fastq | tr -s ' ' | grep -v ' 0 '  | cut -f 3 -d ' ' | grep  'leftover' | grep -v 'total' |  xargs -I {} echo ${bamdir}/{})
bamfiles_="--fastqFile=${bamfiles}"
bamfiles_leftover=$(echo $bamfiles_ | sed 's/ /:/g')
bash ${npTranscript}/scripts/run.sh ${bamfiles_leftover}   --reference=${reference_virus} --annotation ${coord_file_virus} --resdir ${resdir_leftover} ${opts} ${opts1} ${opts2} ${opts3}
fi
