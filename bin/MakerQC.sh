#!/bin/bash

#Script to map fastq files to a reference. 
#Builds the index, performs mapping and notifies when complete. 
#Miniconda AssemblyEnv should be activated. 

#To exit if any of the pipeline fails. 
set -e

#Set the variables
gff=${1? Please eneter a gff3 File.} #Reference Assembley for mapping

echo "Processing gff3..."

# Count the number of gene models:
echo "Gene models:"
cat ${gff} | awk '{ if ($3 == "gene") print $0 }' | awk '{ sum += ($5 - $4) } END { print NR, sum / NR }'

# Count the number of gene models:
echo "Gene models:"
cat ${gff} | awk '{ if ($3 == "gene") print $0 }' | awk '{ sum += ($5 - $4) } END { print NR, sum / NR }'


echo "Done."




