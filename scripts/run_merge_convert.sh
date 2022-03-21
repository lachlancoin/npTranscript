#!/bin/bash
#SBATCH --job-name=npTr_conv
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32000 # mb
#SBATCH --time=3:59:00
#SBATCH --output=%x.%j.stdout
#SBATCH --error=%x.%j.stderr
#SBATCH --cpus-per-task=8
#SBATCH --partition=physical,snowy


##CONVERT TO MATRIZ
module load gcc/8.3.0
module load openmpi/3.1.4
module load r-bundle-bioconductor/3.9-r-3.6.0
#cd $resdir
#Rscript ${npTranscript}/R/convertHDF5Matrix.R
if [ ! $npTranscript ]; then
 npTranscript=/home/lcoin/github/npTranscript
fi

for i in merged_all_lib*; do zcat $i/*transcripts.mod.txt.gz | cut -f 1 ; done  | sort | uniq > all_transcripts.txt

for i in merged_all_lib*; do cd $i; Rscript  ${npTranscript}/R/convertHDFToMatrix.R  transcript_file=../all_transcripts.txt;  cd ..; done

Rscript ${npTranscript}/R/merge_sc.R --prefix=merge_all_lib


