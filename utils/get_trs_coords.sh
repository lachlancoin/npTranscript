#!/bin/bash -e

#script location
npT=`dirname $0`

#MEME SUITE location
MS='/home/lcoin/meme-5.1.1/'

VIRUS='/home/daniel/wuhan_coronavirus_australia.fasta'

BAM='/DataOnline/Data/raw_external/Coronavirus/Direct_RNA_Sequence_Cellular/doherty_analysis/whole_genome_map/sorted.allreads.bam'

#check meme suite is accessible
command -v fimo >/dev/null 2>&1 || { command -v $MS/src/fimo >/dev/null 2>&1; } || { echo >&2 "FIMO not found.  Aborting."; exit 1; }
command -v jaspar2meme >/dev/null 2>&1 || { command -v $MS/scripts/jaspar2meme >/dev/null 2>&1; } || { echo >&2 "jaspar2meme not found.  Aborting."; exit 1; }

#Find TRS. Shift this to other script
#TRS-L coords. Exact string match. Ouputs co-ord file "TRS_short.coords.csv"
python3 $npT/../python/trs-coords.py --virus_genome $VIRUS

#build short motif
#coords of first match (ie. TRS-L)
declare -a coord_array=( `sed '2q;d' TRS_short.coords.csv | sed -e 's/,[A-Z]*$//g' -e 's/,/ /g'` )
echo "Building PWM from coordiates ${coord_array[@]} of reference"

python3 $npT/../python/motif_from_bam.py --jaspar TRS_short.cm --coords "${coord_array[@]}" $BAM


#build long motif
(( coord_array[0]-=5 )) #extend out to include preceeding conserved hexamer
echo "Building PWM from coordiates ${coord_array[@]} of reference"

python3 $npT/../python/motif_from_bam.py --jaspar TRS_long.cm --coords "${coord_array[@]}" $BAM



echo 'Searching genome with FIMO'
$MS/scripts/jaspar2meme -strands 1 -cm . > TRS.meme
$MS/src/fimo --norc --thresh 1e-3 TRS.meme $VIRUS


echo "Finished!"
