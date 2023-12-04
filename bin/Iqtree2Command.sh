#!/bin/bash

#Command for running iqtree2 phylogeny with 1000 bootstrap and generating iqtree2 log. 

Aln=${1?No alinmnent file provided.}
Prefix=${2?No Prefix Provided.}

python -c "print('=' * 75)"
echo "Iqtree2 phylogeny"
echo "-----------------"
echo $(date)
python -c "print('=' * 75)"
echo "iqtree2 version:"
iqtree2 --version
echo "Command:"
echo "iqtree2 -B 1000 -s Fo_TEF.aln"

python -c "print('=' * 75)"
echo "Iqtree2 log"
echo "----------"

#Run the iqtree2 command:
iqtree2 -B 1000 --prefix ${Prefix}  -s ${Aln}

#Signify job has finished and send email notification
python -c "print('=' * 75)"
SendEmail.py "iqtree2 Phylogeny Complete"  "Iqtree2 Phylogeny complete. Check directory for output files" 
python -c "print('=' * 75)"
echo $(date)
python -c "print('=' * 75)"
