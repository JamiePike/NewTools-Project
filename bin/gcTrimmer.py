#!/bin/Biopython

#Extracts sequences from a FASTA file and order in to two different FASTAs depending on GC content

###############################################
#Set up and import modules.
import re, sys, colorama
from datetime import date
from tqdm import tqdm
from colorama import Fore
from Bio import SeqIO, SeqUtils
from Bio.SeqUtils import GC
from Bio.SeqRecord import SeqRecord
###############################################
#Establish input file (FASTA) and gc Threshold.
infile = sys.argv[1]
gcContentThreshold = sys.argv[2]
###############################################
#Establish colour reset.
colorama.init(autoreset=True)
###############################################
# Establish Count for sequences above and below threshold.
SeqAboveCount = 0
SeqBelowCount = 0

###############################################
#Progress:
#Record date
today = date.today()
startTime = today.strftime("%B %d, %Y")
#Create progress bar
num = len([1 for line in open(infile) if line.startswith(">")])
print("="*80)
print(f'Running gcTrimmer.py.\n{startTime}.\n\nInput File:\t\t{infile}\nNumber of Contigs:\t{num}\nGC threshold:\t\t{gcContentThreshold}%\n')

#Start progress bar
with tqdm(total=num) as pbar:
###############################################
    #Open two .fasta files using the name of the infile.
    with open(f'{infile}_GCcontentBelow{gcContentThreshold}perc.fasta'.format(), "a") as LessThanFile, open(f'{infile}_GCcontentAbove{gcContentThreshold}perc.fasta'.format(), "a") as GreaterThanFile, open(f'gcContentReport.txt'.format(), "a") as report  :

    #Parse and loop through the FASTA input file. Search each scaffold/contig for sequences with length => than the desired GC content and deposit them in the FASTA
    #containing the less than the desired gc content. If the sequence gc content is greater than the gc content threshold, then it is added to the "greater than" FASTA file.
        for sequence in SeqIO.parse(infile, "fasta"):
            print("="*10, file=report, sep="\n")
            gccontent = (sequence.seq.count('G') + sequence.seq.count('C')) / len(sequence) * 100
            pbar.update(1)
            if gccontent < int(gcContentThreshold):
                print(f'{sequence.description} is below the threshold GC content. GC content is {Fore.RED}{gccontent}{Fore.BLACK}.\n', file=report, sep="\n")
                print(sequence.format("fasta"), file=LessThanFile, sep="\n")
                SeqBelowCount = SeqBelowCount + 1
            else:
                print(f'{sequence.description} is above the threshold GC content. GC content is {Fore.GREEN}{gccontent}{Fore.BLACK}.\n', file=report, sep="\n")
                print(sequence.format("fasta"), file=GreaterThanFile, sep="\n")
                SeqAboveCount = SeqAboveCount + 1


#Print the total number of sequences above and below the GC threshold.

print(f'\nThe number of sequences above the GC threshold: {Fore.GREEN}{SeqAboveCount}{Fore.BLACK}.\nThe number of sequences below the GC threshold: {Fore.RED}{SeqBelowCount}{Fore.BLACK}.')
print("="*80)
