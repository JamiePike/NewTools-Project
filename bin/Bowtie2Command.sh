#!/bin/bash

#Script to map fastq files to a reference. 
#Builds the index, performs mapping and notifies when complete. 
#Miniconda AssemblyEnv should be activated. 

#To exit if any of the pipeline fails. 
set -e

#Set the variables
Ref=${1? Please eneter a reference Assembly.} #Reference Assembley for mapping
One=${2?Please provide the fastq 1.} #Where One is the first fastq.
Two=${3?Please provide the fastq 2.} #Where two is the second fastq.
Prefix=${4?Please provide an output prefix} #Prefix for outputs


if [ "$CONDA_DEFAULT_ENV" = "AssemblyEnv" ]; then #Check the conda environment is correct.
    python -c "print('=' * 78)"
    echo "Bowtie2 Mapping"
    echo "---------------"
    echo $(date)
    echo "Usage: Bowtie2Command.sh <reference fasta> <fastq_1> <fastq_2> <output prefix>"
    python -c "print('=' * 78)"
    echo "Building Bowtie2 Index..."
    #Make Ouput Directory for index
    mkdir ${Prefix}_Bowtie2Index
    #Build the index 
    bowtie2-build ${Ref} ${Prefix}_Bowtie2Index/${Prefix}.Index >${Prefix}_Bowtie2-Build.log

    echo "Performing Bowtie2 Alignment..."
    #Perform Bowtie2 Mapping using S. maltophilia as the reference.
    bowtie2 --local -x ./${Prefix}_Bowtie2Index/${Prefix}.Index  -1 $One -2 $Two | samtools view -Shu | samtools sort -o ${Prefix}.sorted.bam 

    python -c "print('=' * 78)"
    #Send email notifciation of job completion.
    SendEmail.py "Bowtie2 Complete" "Bowtie2 Mapping complete for ${Prefix}. Check for results."
    python -c "print('=' * 78)"      
else
    echo "You have not activated the conda environment."
    exit 1
fi