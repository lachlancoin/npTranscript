#!/bin/bash -e

declare -a READS=${@}
echo "Mapping reads from: ${READS[@]}"
MM2='/sw/minimap2/minimap2-2.11_x64-linux/minimap2'
SCRIPT_DIR=$(dirname $0)

#map to whole genome
$MM2 -ax splice -un $SCRIPT_DIR/../data/SARS-Cov2/wuhan_coronavirus_australia.fasta.gz $READS | samtools view -hF4 -b | samtools sort - > sorted.whole_genome_mapped.bam

#map to leaderr
$MM2 -ax map-ont -un $SCRIPT_DIR/../data/SARS-Cov2/leader.wuhan_coronavirus_australia.fasta.gz $READS | samtools view -hF4 -b | samtools sort - > sorted.leader_mapped.bam

