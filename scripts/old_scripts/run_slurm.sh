#!/bin/sh

#SBATCH --job-name=npTr_sars
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=122800 # mb
#SBATCH --time=100:00:00
#SBATCH --output=corona1.stdout
#SBATCH --error=corona1.stderr
#SBATCH --cpus-per-task=16
#SBATCH --partition=physical,snowy


## Load required modules
module load gcc/8.3.0
module load samtools/1.9
module load minimap2/2.17
module load java

##script for running 
##tip - use symbolic link to put this in the directory with bam files
#run as sbatch run_slurm_combined.sh todo.txt human combined --RNA=true
#  sbatch run_slurm_combined.sh human combined --RNA=false
#to find bams for todo.txt
#find /data/gpfs/projects/punim1068/host_analysis/direct_RNA_sequin_genome/Vero/ -name *bam > todo.txt
export JSA_MEM=100000m

refdir=$1
shift

loc=`scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2}'`
echo $loc
export npTranscript=`echo $loc | sed 's/\(npTranscript\).*/\1/'  `
echo ${npTranscript}

if [ $1 ]; then
todo=$1
shift
else
todo="fdagdfdgzgfdd"
fi

if [ -f $todo ]; then
        bamfiles=$(cat $todo | sed "s/[[:space:]]\+/\n/g" | grep ".bam$")
else
        tag=".bam"
        bamfiles=$(find . -maxdepth 1  -size +0b | grep "${tag}$" )
fi
bamfiles1="--bamFile=$(echo $bamfiles | sed 's/^M//g' | sed 's/^ //g' | sed 's/ /:/g')"
echo $bamfiles1
opts1="" 
a=$(ls  | grep '^reads_in.txt$'  | wc -l)
if [ $a -eq 1 ]; then 
  opts1="--readList=reads_in.txt" 
fi

#reference="${npTranscript}/data/SARS-Cov2/NC/NC_045512.fasta.gz" 
#coord_file="${npTranscript}/data/SARS-Cov2/NC/Coordinates.csv"
#refdir="SARS-Cov2/VIC01"
#refdir="229E_CoV"
reference=$(find "${npTranscript}/data/${refdir}/" | grep  'fasta.gz')
coord_file="${npTranscript}/data/${refdir}/Coordinates.csv"

if [ ! -f $reference ]; then
echo "usage requires refdir eg SARS-Cov2/VIC01"
fi
#alias abpoa='/sw/abpoa/v1.0.1/abpoa'
dat=$(date +%Y%m%d%H%M%S)
echo $dat
resdir="results_${dat}"
##note that to do separate msa for multiple input bams you type --doMsa=5_3:sep
#opts="--bin=100 --RNA true --breakThresh=1000  --isoformDepthThresh=10000 --coverageDepthThresh=0 --extra_threshold=200 --msaDepthThresh=20 --doMSA=all:sep --reAlignExtra=true"
opts="--bin=100 --library=./ --writePolyA=true --RNA=true --breakThresh=500  --isoformDepthThresh=10000 --coverageDepthThresh=1 --extra_threshold=200 --msaDepthThresh=20 --doMSA=false --reAlignExtra=true"
opts2="-gffThresh=10:10 --writeGFF=true --attempt5rescue=false  --fail_thresh=7 --recordDepthByPosition=true --maxThreads=14"
#opts2="${opts2} --readList=readToCluster.txt.gz" 
#opts="${opts} --maxReads 10000"
bash ${npTranscript}/scripts/run.sh ${bamfiles1}   --reference=${reference} --annotation=${coord_file} --resdir=${resdir} ${opts} ${opts1} ${opts2} $@

cd ${resdir}
#bash ${npTranscript}/scripts/consensus.sh 0
#echo "now running R script"
#Rscript ${npTranscript}/R/npTranscript.R  ${npTranscript}/data/SARS-Cov2/VIC01
