#!/bin/bash -e

MM2='/sw/minimap2/minimap2-2.11_x64-linux/minimap2'

READS="0.consensus.fasta"
N_READS=$(grep ">" 0.consensus.fasta | wc -l)
echo "Remapping "${N_READS}" consensus sequences"
MM2='/sw/minimap2/minimap2-2.11_x64-linux/minimap2'
npT=$(dirname $0)
HOST='/DataOnline/Data/raw_external/Coronavirus/monkey/newdb/Chlorocebus_sabaeus.ChlSab1.1.dna.toplevel.fa.gz'
VIRUS="$npT/../data/SARS-Cov2/wuhan_coronavirus_australia.fasta.gz"

#map to whole genome. Include all mappings
$MM2 -ax splice -un $npT/../data/SARS-Cov2/wuhan_coronavirus_australia.fasta.gz $READS | samtools view -h -b | samtools sort - > consensus_remapped.bam


