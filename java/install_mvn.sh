mvn install:install-file -Dfile=./libs/japsa-1.9.4e.jar  -DgroupId=coin -DartifactId=japsa         -Dversion=1.9.4e        -Dpackaging=jar
mvn install:install-file -Dfile=./libs/sis-base-18.09.0.jar  -DgroupId=ch.ethz.wiki-bsse -DartifactId=sis-base     -Dversion=18.09.0      -Dpackaging=jar
mvn install:install-file -Dfile=./libs/sis-jhdf5-19.04.0.jar  -DgroupId=ch.ethz.wiki-bsse -DartifactId=sis-jhdf5         -Dversion=19.04.0      -Dpackaging=jar


mvn clean package install
