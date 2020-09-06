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


##script for running
##tip - use symbolic link to put this in the directory with bam files
export JSA_MEM=30000m
if [ ! $npTranscript ] ; then
        npTranscript=${HOME}/github/npTranscript
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


##VIRAL ANALYSIS ON VIRAL READS
opts="--bin=100 --breakThresh=1000 --gffThresh=100  --isoformDepthThresh=10000 --coverageDepthThresh=0 --extra_threshold=200 --msaDepthThresh=20 --doMSA=all:sep --reAlignExtra=true"
opts2="--fail_thresh=0 --recordDepthByPosition=true"
opts3="--mm2_path=${mm2_path}"

tag=".primary.fastq"
resdir_virus="results_${tag}"
#bamfiles_virus=$(bash ${npTranscript}/scripts/getInputFiles.sh ${tag})
bamfiles=$(find . -maxdepth 1 -type f,l -size +0b | grep "${tag}$" )
bamfiles_virus="--fastqFile=$(echo $bamfiles | sed 's/ /:/g')"

bash ${npTranscript}/scripts/run.sh ${bamfiles_virus}   --reference=${reference_virus} --annotation ${coord_file_virus} --resdir ${resdir_virus} ${opts} ${opts1} ${opts2} ${opts3}
cd ${resdir_virus}
bash ${npTranscript}/scripts/run_slurm_leftover.sh $species
cd ..

#tag=".secondary.fastq"
#resdir_virus="results_${tag}"
#bamfiles_virus=$(bash ${npTranscript}/scripts/getInputFiles.sh ${tag})
#bamfiles=$(find . -maxdepth 1 -type f,l -size +0b | grep "${tag}$" )
#bamfiles_virus="--fastqFile=$(echo $bamfiles | sed 's/ /:/g')"
#bash ${npTranscript}/scripts/run.sh ${bamfiles_virus}   --reference=${reference_virus} --annotation ${coord_file_virus} --resdir ${resdir_virus} ${opts} ${opts1} ${opts2} ${opts3}

tag=".supplementary.fastq"
resdir_virus="results_${tag}"
#bamfiles_virus=$(bash ${npTranscript}/scripts/getInputFiles.sh ${tag})
bamfiles=$(find . -maxdepth 1 -type f,l -size +0b | grep "${tag}$" )
bamfiles_virus="--fastqFile=$(echo $bamfiles | sed 's/ /:/g')"
bash ${npTranscript}/scripts/run.sh ${bamfiles_virus}   --reference=${reference_virus} --annotation ${coord_file_virus} --resdir ${resdir_virus} ${opts} ${opts1} ${opts2} ${opts3}
cd ${resdir_virus}
bash ${npTranscript}/scripts/run_slurm_leftover.sh $species
cd ..
