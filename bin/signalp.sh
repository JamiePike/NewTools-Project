#!/bin/bash

for i in *.fasta ; do 
    echo "Processing ${i}..." ;  
    signalp -fasta ./${i} 1>./${i}-signalp.out 2>./${i}-signalp.log #Identify sequences with signal peptide.
    echo "Generating index files..." ;
    awk '$2 ~ /SP\(Sec\/SPI\)/ {print $1}' ./${i}/${i}-Translated_summary.signalp5 > ./${i}-signalPepPositives.list;
    echo "Generating fasta...";
    for j in $(cat ./${i}-signalPepPositives.list); do 
        samtools faidx ./${i} ${j}; 
        done > ./${i}-signalp.filtered.fasta;
    echo "done." ;
done