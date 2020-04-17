#!/bin/bash -e

declare -a READS
for read_file in ${@}; do READS+=($(readlink -f $read_file)); done

echo $READS

MM2='/sw/minimap2/minimap2-2.11_x64-linux/minimap2'
ORIG_DIR=$(pwd)

cd $(dirname $0)


#map to whole genome
$MM2 -ax splice -un ../data/SARS-Cov2/wuhan_coronavirus_australia.fasta.gz $READS | samtools view -hF4 -b | samtools sort - > $ORIG_DIR/sorted.whole_genome_mapped.bam

#map to leaderr
$MM2 -ax map-ont -un ../data/SARS-Cov2/leader.wuhan_coronavirus_australia.fasta.gz $READS | samtools view -hF4 -b | samtools sort - > $ORIG_DIR/sorted.leader_mapped.bam


cd $ORIG_DIR
