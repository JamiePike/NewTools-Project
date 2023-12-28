#!/bin/bash

for i in *.fasta ; do 
    echo "processing ${i}..." ;
    EffectorP.py -i ./${i}-signalp.filtered.fasta -E ../effectorp/${i}-EffectorP.filtered.fasta > ../effectorp/${i}-EffectorP.filtered.log ; 
    echo "done." ;
done 