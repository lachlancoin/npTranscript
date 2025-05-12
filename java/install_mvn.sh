#!/bin/bash
module load Java
#module load java/1.8.0_241
#module load java
module load Maven/3.9.3
#module load maven/3.6.3
reinstall=$1
#module load anaconda2/2019.10
#conda activate build
git pull

if [ $reinstall ]; then

 mvn install:install-file -Dfile=./libs/japsacov-1.9.6a.jar -DgroupId=coin -DartifactId=japsacov         -Dversion=1.9.6a        -Dpackaging=jar
 mvn install:install-file -Dfile=./libs/sis-base-18.09.0.jar  -DgroupId=ch.ethz.wiki-bsse -DartifactId=sis-base     -Dversion=18.09.0      -Dpackaging=jar
 mvn install:install-file -Dfile=./libs/sis-jhdf5-19.04.0.jar  -DgroupId=ch.ethz.wiki-bsse -DartifactId=sis-jhdf5         -Dversion=19.04.0      -Dpackaging=jar
fi


mvn clean package install


#export mainclass=npTranscript.run.Barcodes
#bash ../scripts/run.sh
#export mainclass=npTranscript.run.ViralTranscriptAnalysisCmd2
