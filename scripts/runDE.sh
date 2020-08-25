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


##USAGE bash runDE.sh control infected

#pwd=$(pwd)
#cd ~/github/npTranscript
#git pull
#cd $pwd
params="np.control=control np.case=infected np.exclude=none np.analysis=betabinom np.depth_thresh=100 np.isoformDepth=100 np.isoform.test=chisq.test np.dm.test=chisq.test np.maxIsoformGroups=5 np.adjustMethod=BH"
#"np.prefix_keep=ENSC
#"np.prefix_remove"="^[0-9]{1,}\\."
#"np.prefix_sequins"="^R[0-9]_"
params1=""
if [ $1 ]; then
 if [ $1 == "install" ]; then
	params1="np.install=TRUE"
 fi
fi
#"np.libs_to_install"

Rscript --vanilla   ~/github/npTranscript/R/npDE.R  ${params} ${params1}

