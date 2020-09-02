#!/bin/bash

#SBATCH --job-name=npTr_comb
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

if [ ! $npTranscript ] ; then
export  npTranscript=${HOME}/github/npTranscript
fi
if [ ! $reference_virus ]; then
  export reference_virus="${npTranscript}/data/SARS-Cov2/VIC01/wuhan_coronavirus_australia.fasta.gz"
  export coord_file_virus="${npTranscript}/data/SARS-Cov2/VIC01/Coordinates.csv"
fi

species=$1
shift;
if [ ! $species ]; then
        echo "first option needs to be 'monkey'  or ' human' "
        exit;
fi
if [ $species == "human" ]; then
        reference="/DataOnline/Data/Projects/corona_invitro/host_analysis/db/merged/human_virus_sequin_ensembl_pri_merged_genome.fasta"
        coord_file="/DataOnline/Data/Projects/corona_invitro/host_analysis/db/human/ensembl/Homo_sapiens.GRCh38.100.gtf.gz"
#coord_file="/home/lcoin/Homo_sapiens.GRCh38.100.gtf.gz"
        GFF_features="--GFF_features=gene_name:description:gene_id:gene_biotype:gene_id"
elif [ $species == "monkey" ]; then
        reference="/DataOnline/Data/Projects/corona_invitro/host_analysis/db/merged/monkey_virus_sequin_genome.fasta"
        coord_file="/DataOnline/Data/raw_external/Coronavirus/monkey/newdb/Chlorocebus_sabaeus.ChlSab1.1.99.gff3.gz"
        GFF_features="-GFF_features=Name:description:ID:biotype:Parent"
else
        echo "first option needs to be 'monkey'  or ' human' "
        exit
fi
echo $reference


MM2='/sw/minimap2/minimap2-2.11_x64-linux/minimap2'
npT=$npTranscript
HOST="/DataOnline/Data/Projects/corona_invitro/host_analysis/db/merged/human_virus_sequin_ensembl_pri_merged_genome.mmi"
tag=".leftover.fastq"

grep '^@\w\w\w\w\w\w\w\w-' *.leftover.fastq  | cut -d ':' -f 2  | cut -f 1,3 -d ' '  | sed 's/ /        /g' | sed 's/@//g' > leftover_readIds.txt
du -h *fastq | grep  '^0 '  | cut -f 2 > deleted_empty_fastq.txt
du -h *fastq | grep  '^0 '  | cut -f 2 | xargs -I {} rm -f {}

READS=$(find . -maxdepth 1 -type f,l -size +0b | grep "${tag}$")
#map to host
for i in $READS; do
	out=$i
	echo  $i $out
$MM2 -ax splice -uf $HOST $i | samtools view -hF4 -b | samtools sort - > $out.sorted.bam
done

tag=".leftover.fastq.sorted.bam"
bamfiles=$(find . -maxdepth 1 -type f,l -size +0b | grep "${tag}$" )
bamfiles1="--bamFile=$(echo $bamfiles | sed 's/ /:/g')"
resdir="results_${tag}"
opts="--bin=100 --breakThresh=100 --coronavirus=false --maxThreads=8 --extra_threshold=500 --writePolyA=true --msaDepthThresh=1000 --doMSA=false --numExonsMSA=1:2:3:4:5 --msa_source=RNA --useExons=true --span=protein_coding --includeStart=false --isoformDepthThresh 50"

#for dRNA datasets
opts="${opts} --RNA=true --readList=leftover_readIds.txt"
opts2="--fail_thresh=0  --recordDepthByPosition=true"
opts1=""
echo $opts $opts1 $opts2 $@
bash ${npTranscript}/scripts/run.sh ${bamfiles1}   --reference=${reference} --annotation=${coord_file} --resdir=${resdir} ${opts} ${opts1} ${opts2} ${GFF_features} $@

