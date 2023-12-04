#!/bin/bash

#Jamie Pike
#Command for filtering conatminant blobtools hits. This groups all of the standard BlobTools commands I used to remove a contaminat contigs analysis. 

python -c "print('=' * 75)"
echo "Blobtools Contaminant Filter"
echo "----------------------------"
echo $(date)
echo "Usage: ContaminantFilter.sh <FASTA file> <blobtools json file> <outfile prefix>"
python -c "print('=' * 75)"

infile=${1?Please provide the assembly input fasta.} #Input Assembly.
json=${2?Please provide a BlobTools json file.} #Input BAM.
prefix=${3?Please provide a prefix for you outputs.} 

########################
#Escape the script if there are any errors. 
set -e 

echo "Creating species table..."
#Use BlobTools view to generate a blobtools table.txt file for filtering.
blobtools view -i ${json} -r species -o ${prefix}

echo "Filetering for Fusarium and no-hit only contigs..."
#Use grep to extact only the Fusarium and no-hit lines from the hit table.
grep -E 'Fusarium|no-hit' ${prefix}*.table.txt > ${prefix}.FusariumHits.table.txt

echo "Generating list file..."
#Create a list of the Fusarium and no-hit only contigs. 
awk '{print $1}' ${prefix}.FusariumHits.table.txt >  ${prefix}.FusariumHits.list.txt

echo "Generating contaminant filtered fasta..."
#Create a list of the Fusarium and no-hit only contigs. 
blobtools seqfilter -i ${infile} -l ${prefix}.FusariumHits.list.txt -o ${prefix}.ContaminantFiltered

echo "Generating contaminant sequences fasta..."
#Create a list of the Fusarium and no-hit only contigs. 
blobtools seqfilter -v -i ${infile} -l ${prefix}.FusariumHits.list.txt -o ${prefix}.ContaminantSequences

echo "ContaminantFilter.sh finished."
echo $(date)
echo "+++++
It is advisable that you check the number of contigs in the ${prefix}.ContaminantFiltered.fasta matches the number of lines in the ${prefix}.FusariumHits.table.txt file."
python -c "print('=' * 75)" 



