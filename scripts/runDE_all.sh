#!/bin/bash

#SBATCH --job-name=DE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=9000 # mb
#SBATCH --time=100:00:00
#SBATCH --output=DE.stdout
#SBATCH --error=DE.stderr
#SBATCH --cpus-per-task=8
#SBATCH --partition=physical,snowy

# USAGE:
# sbatch --job-name=DE  runDE.sh

## Load required modules
module load gcc/8.3.0 
module load openmpi/3.1.4
module load r-bundle-bioconductor/3.9-r-3.6.0

##USAGE sbatch runDE.sh np.control=control np.case=infected np.virus=TRUE

if [ ! $npTranscript ] ; then
#export npTranscript=${HOME}/github/npTranscript
export npTranscript=/data/gpfs/projects/punim1068/npTranscript
fi 


#pwd=$(pwd)
#cd ~/github/npTranscript
#git pull
#cd $pwd
params="np.analysis=betabinom np.depth_thresh=100 np.isoformDepth=100 np.isoform.test=chisq.test np.dm.test=chisq.test np.maxIsoformGroups=5 np.adjustMethod=BH"
#params1="np.control=control np.case=infected"
#"np.prefix_keep=ENSC
#"np.prefix_remove"="^[0-9]{1,}\\."
#"np.prefix_sequins"="^R[0-9]_"
#"np.install=TRUE"
#"np.libs_to_install"

libdir=/data/gpfs/projects/punim1068/Rlib
pwd=$(pwd)

dirs=$(find . -name results_.primary.fastq*)
for i in $dirs; do
cd $i
Rscript --vanilla   ${npTranscript}/R/npDE.R  ${params}  np.libdir=${libdir} np.source=${npTranscript}/R $@
cd $pwd
done
