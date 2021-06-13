#!/bin/bash


module load java
module load maven/3.6.3
mvn install:install-file -Dfile=./libs/japsacov-1.9.5e.jar -DgroupId=coin -DartifactId=japsacov         -Dversion=1.9.5e        -Dpackaging=jar
mvn install:install-file -Dfile=./libs/sis-base-18.09.0.jar  -DgroupId=ch.ethz.wiki-bsse -DartifactId=sis-base     -Dversion=18.09.0      -Dpackaging=jar
mvn install:install-file -Dfile=./libs/sis-jhdf5-19.04.0.jar  -DgroupId=ch.ethz.wiki-bsse -DartifactId=sis-jhdf5         -Dversion=19.04.0      -Dpackaging=jar


mvn clean package install
