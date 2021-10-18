#!/bin/bash


module load java
module load maven/3.6.3
mvn install:install-file -Dfile=./libs/japsacov-1.9.5e.jar -DgroupId=coin -DartifactId=japsacov         -Dversion=1.9.5e        -Dpackaging=jar
mvn install:install-file -Dfile=./libs/sis-base-18.09.0.jar  -DgroupId=ch.ethz.wiki-bsse -DartifactId=sis-base     -Dversion=18.09.0      -Dpackaging=jar
mvn install:install-file -Dfile=./libs/sis-jhdf5-19.04.0.jar  -DgroupId=ch.ethz.wiki-bsse -DartifactId=sis-jhdf5         -Dversion=19.04.0      -Dpackaging=jar


jardir=~/github/JNAerator/jnaerator/target
java -jar ${jardir}/jnaerator-0.13-SNAPSHOT-shaded.jar -mode Jar -runtime JNA ~/github/edlib/edlib/include/edlib.h ~/github/edlib/meson-build/libedlib.so
mvn install:install-file -Dfile=out.jar  -DgroupId=edlib -DartifactId=edlib         -Dversion=1.2.7      -Dpackaging=jar


mvn clean package install

jar xf out.jar lib/linux_x64/libedlib.so
mkdir target/classes/linux-x86-64
mv lib/linux_x64/libedlib.so target/classes/linux-x86-64
rm -rf lib
rm out.jar
