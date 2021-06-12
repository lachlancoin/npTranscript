#!/bin/bash

#SBATCH --job-name=VERO
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=31800 # mb
#SBATCH --time=100:00:00
#SBATCH --output=corona1.stdout
#SBATCH --error=corona1.stderr
#SBATCH --cpus-per-task=8
#SBATCH --partition=physical,snowy


## Load required modules
#module load samtools/1.9
module load minimap2/2.17
module load java


##script for running
##tip - use symbolic link to put this in the directory with bam files
export JSA_MEM=30000m
if [ ! $npTranscript ] ; then
	export npTranscript=/data/gpfs/projects/punim1068/npTranscript
fi
if [ ! $reference_virus ]; then
  export reference_virus="${npTranscript}/data/SARS-Cov2/VIC01/wuhan_coronavirus_australia.fasta.gz"
  export coord_file_virus="${npTranscript}/data/SARS-Cov2/VIC01/Coordinates.csv"
fi

#species only used for the leftover remapping
species=$1
shift;

if [ ! $species ]; then
	species="human"
fi

dat=$(date +%Y%m%d%H%M%S)
mm2_path="/sw/minimap2/current/minimap2"
mm2_path="minimap2"
opts_virus="../opts_virus.txt"
if [ ! -f $opts_virus ]; then
  opts_virus="${npTranscript}/opts_virus.txt"
fi

opts=$(grep -v '^#' ${opts_virus})


##VIRAL ANALYSIS ON VIRAL READS
#opts="--bin=100 --breakThresh=1000 --gffThresh=100:100  --isoformDepthThresh=10000 --coverageDepthThresh=0 --extra_threshold=200 --msaDepthThresh=20 --doMSA=false --reAlignExtra=true"
#opts2="--recordDepthByPosition=false --maxThreads=8 --qualThresh=7 --trainStrand=true"
# opts3="--mm2_path=${mm2_path}"

tag=".primary.fastq"
tag=".fastq" 
resdir_virus="results_${tag}"
#bamfiles_virus=$(bash ${npTranscript}/scripts/getInputFiles.sh ${tag})
bamfiles=$(find . -maxdepth 1  -size +0b | grep "${tag}$" )
bamfiles_virus="--fastqFile=$(echo $bamfiles | sed 's/ /:/g')"

bash ${npTranscript}/scripts/run.sh ${bamfiles_virus}   --reference=${reference_virus} --annotation ${coord_file_virus} --resdir ${resdir_virus} 
#cd ${resdir_virus}
#bash ${npTranscript}/scripts/run_slurm_leftover.sh $species $@
#cd ..

#tag=".secondary.fastq"
#resdir_virus="results_${tag}"
#bamfiles_virus=$(bash ${npTranscript}/scripts/getInputFiles.sh ${tag})
#bamfiles=$(find . -maxdepth 1 -type f,l -size +0b | grep "${tag}$" )
#bamfiles_virus="--fastqFile=$(echo $bamfiles | sed 's/ /:/g')"
#bash ${npTranscript}/scripts/run.sh ${bamfiles_virus}   --reference=${reference_virus} --annotation ${coord_file_virus} --resdir ${resdir_virus} ${opts} ${opts1} ${opts2} ${opts3}

tag=".supplementary.fastq"
resdir_virus="results_${tag}"
#bamfiles_virus=$(bash ${npTranscript}/scripts/getInputFiles.sh ${tag})
bamfiles=$(find . -maxdepth 1 -size +0b | grep "${tag}$" )
bamfiles_virus="--fastqFile=$(echo $bamfiles | sed 's/ /:/g')"
#bash ${npTranscript}/scripts/run.sh ${bamfiles_virus}   --reference=${reference_virus} --annotation ${coord_file_virus} --resdir ${resdir_virus} ${opts} ${opts1} ${opts2} ${opts3} $@
#cd ${resdir_virus}
#bash ${npTranscript}/scripts/run_slurm_leftover.sh $species $@
#cd ..
