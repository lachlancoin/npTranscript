#!/bin/bash




#SBATCH --job-name=npTr_comb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=31800 # mb
#SBATCH --time=100:00:00
#SBATCH --output=corona1.stdout
#SBATCH --error=corona1.stderr
#SBATCH --cpus-per-task=16
#SBATCH --partition=physical,snowy


## Load required modules
module load samtools/1.9
module load minimap2/2.17
module load java

##script for running 
##tip - use symbolic link to put this in the directory with bam files
#run as sbatch run_slurm_combined.sh todo.txt human combined --RNA=true
#  sbatch run_slurm_combined.sh human combined --RNA=false
#NOTE: need to have an opts_host.txt in the run directory
#to find bams for todo.txt
#find /data/gpfs/projects/punim1068/host_analysis/direct_RNA_sequin_genome/Vero/ -name *bam > todo.txt
export JSA_MEM=30000m

if [ ! $npTranscript ] ; then
#export	npTranscript=${HOME}/github/npTranscript
export npTranscript=/data/gpfs/projects/punim1068/npTranscript
fi 
if [ ! $reference_virus ]; then
  export reference_virus="${npTranscript}/data/SARS-Cov2/VIC01/wuhan_coronavirus_australia.fasta.gz"
  export coord_file_virus="${npTranscript}/data/SARS-Cov2/VIC01/Coordinates.csv"
fi

todo=$1
shift;
species=$1
shift;
if [ ! $species ]; then 
	echo "first option needs to be 'monkey'  or ' human' "
	exit;
fi

db_path=/data/gpfs/projects/punim1068/db

	
if [ $species == "human" ]; then
	reference="${db_path}/merged/sequin_genome/VIC_virus_human_sequin_genome_ensembl_pri_merged.fasta"
	coord_file="${db_path}/merged/human_ensembl_sequin_merged.gtf"
elif [ $species == "monkey" ]; then
	reference="${db_path}/merged/sequin_genome/VIC_virus_monkey_sequin_genome_merged.fasta"
	coord_file="${db_path}/merged/monkey_sequin_merged.gtf"
elif [ $species == "human_only" ]; then
	reference="./GRCh38_full_analysis_set_plus_decoy_hla.fa.gz"
	coord_file="./Homo_sapiens.GRCh38.100.gtf.gz"
else 
	echo "first option needs to be 'monkey'  or ' human' "
	exit
fi
echo $reference



dat=$(date +%Y%m%d%H%M%S)
mm2_path="/sw/minimap2/current/minimap2"
mm2_path="minimap2"
mode=$1
shift
if [ ! $mode ];then
 mode="combined"
fi

##SPECIFY LOCATION OF COMBINED AND VIRUS ONLY DB
cov_chr=$(zcat ${reference_virus} | head -n 1 | cut -f 1 -d ' ' | sed 's/>//g')
echo "coronavirus chr id ${cov_chr}" 
#chroms_to_include="chrIS:${cov_chr}" 
chroms_to_include="all" ## this includes all chromosomes 
chroms_to_ignore="none"   ##
resdir="results_${dat}"
#opts="--bin=100 --breakThresh=100 --coronavirus=false --maxThreads=13 --extra_threshold=2000 --writePolyA=false --msaDepthThresh=1000 --doMSA=false --msa_source=RNA --useExons=true --span=protein_coding --includeStart=false --isoformDepthThresh=50 --chroms_to_ignore=${chroms_to_ignore} --chroms_to_include=${chroms_to_include}"

opts_host="opts_host.txt"
if [ ! -f $opts_host ]; then
  opts_host="${npTranscript}/opts_host.txt"
fi

opts=$(grep -v '^#' ${opts_host})

#for dRNA datasets
#opts="${opts} --RNA=true"
#opts2="--chromsToRemap=${cov_chr}  --mm2_memory=10g --writeIsoforms=true --writeGFF=true --recordDepthByPosition=false --coverageDepthThresh=100 --qualThresh=7 --trainStrand=true"
#echo $opts
if [ -f $todo ]; then
	bamfiles=$(cat $todo | sed "s/[[:space:]]\+/\n/g" | grep ".bam$")
else
	tag=".bam"
	bamfiles=$(find . -maxdepth 1  -size +0b | grep "${tag}$" )
fi
bamfiles1="--bamFile=$(echo $bamfiles | sed 's/^M//g' | sed 's/^ //g' | sed 's/ /:/g')"
#bamfiles1=$(bash ${npTranscript}/scripts/getInputFiles.sh ${tag})

bash ${npTranscript}/scripts/run.sh ${bamfiles1} --chromsToRemap=${cov_chr}  --reference=${reference} --annotation=${coord_file} --resdir=${resdir} ${opts}

cd ${resdir}

##following reassigns any leftover portions of original reads
#bash ${npTranscript}/scripts/run_slurm_leftover.sh $species


if [ $mode == "combined" ]; then
	bash ${npTranscript}/scripts/run_slurm_virus_fastq.sh $species $@
	
fi




