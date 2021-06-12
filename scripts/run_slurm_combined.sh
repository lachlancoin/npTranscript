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
#run as sbatch run_slurm_combined.sh human combined --RNA=true
#  sbatch run_slurm_combined.sh human combined --RNA=false
#NOTE: need to have an opts_run.txt in the run directory
#to find bams for todo.txt
#find /data/gpfs/projects/punim1068/host_analysis/direct_RNA_sequin_genome/Vero/ -name *bam > todo.txt
export JSA_MEM=30000m

if [ ! $npTranscript ] ; then
#export	npTranscript=${HOME}/github/npTranscript
export npTranscript=/data/gpfs/projects/punim1068/npTranscript_streaming
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
elif [ $species == "SARSCov2" ]; then
	reference="${npTranscript}/data/SARS-Cov2/wuhan/wuhan_coronavirus.fasta.gz"
	coord_file="${npTranscript}/data/SARS-Cov2/wuhuan/Coordinates.csv"
elif [ $species == "U13369" ]; then
	#reference="${npTranscript}/data/U13369/Human_ribosomal_DNA_complete_repeating_unit.fasta.gz"
	reference="${npTranscript}/data/U13369/18S_28S.fa.gz"
elif [ $species == "combined" ]; then
	reference="${npTranscript}/data/U13369/combined.fa.gz"
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
#cov_chr=$(zcat ${reference_virus} | head -n 1 | cut -f 1 -d ' ' | sed 's/>//g')
#echo "coronavirus chr id ${cov_chr}" 


opts_host="opts_run.txt"
if [ ! -f $opts_host ]; then
  echo "need to have a file opts_run.txt"
  exit ;
fi

todo="todo.txt"
if [ ! -f $todo ]; then
  echo "need to have a file todo.txt which lists input files"
  exit ;
fi

dat=$(date +%Y%m%d%H%M%S)
resdir="res_${species}_${dat}"

bash ${npTranscript}/scripts/run.sh   --reference=${reference} --resdir=${resdir} --optsFile=${opts_host} --todo=${todo}



opts_annot="opts_annot.txt"

if [ ! -f $opts_annot ]; then
  echo "need to have a file opts_annot.txt for annotation step"
  exit ;
fi

h5file=${resdir}"/0.isoforms.h5"

if [ ! -f $h5file ]; then
  echo "h5 file not created"
  exit ;
fi

if [ ! $coord_file ]; then
  echo "no coordinate file"
  exit ;
fi


bash ${npTranscript}/scripts/run_annot.sh  --reference=${reference} --annotation=${coord_file} --resdir=${resdir} --optsFile=${opts_annot} --h5file=${h5file}




