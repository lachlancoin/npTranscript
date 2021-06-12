#!/bin/bash -e

declare -a READS=${@}
echo "Mapping reads from: ${READS[@]}"
MM2='/sw/minimap2/minimap2-2.11_x64-linux/minimap2'
npT=$(dirname $0)
HOST=''


#map to whole genome
$MM2 -ax splice -un $npT/../data/SARS-Cov2/VIC01/wuhan_coronavirus_australia.fasta.gz $READS | samtools view -hF4 -b | samtools sort - > whole_genome_mapped.bam

#map to host
$MM@ -ax splice -uf $HOST $READS | samtools view -hF4 -b | samtools sort - > host_mapped.bam
