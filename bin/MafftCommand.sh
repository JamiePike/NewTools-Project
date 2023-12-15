#!/bin/bash

#Command for running MAFFT alignment and generating MAFFT log. 

infile=${1?Please ensure you have provided an input fasta.}

python -c "print('=' * 75)"
echo "MAFFT alignment"
echo "---------------"
echo $(date)
python -c "print('=' * 75)"
echo "MAFFT version:"
mafft --version 
echo "MAFFT Command:"
echo "mafft --adjustdirectionaccurately --reorder " ${infile} ">"  ${infile}".aln"

python -c "print('=' * 75)"
echo "MAFFT log"
echo "---------"

#Run the MAFFT command:
mafft --adjustdirectionaccurately --reorder ${infile} > ${infile}.aln

#Signify job has finished and send email notification
python -c "print('=' * 75)"
SendEmail.py "MAFFT Alignment Complete"  "Mafft Alignment complete. Check directory for alignment file" 
python -c "print('=' * 75)"
echo $(date)
python -c "print('=' * 75)"
