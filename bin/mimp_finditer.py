#!/bin/python

#### MIMP FINDER - SEARCHES FOR ALL MIMP IN SEQUENCE ####

#Search through a FASTA file and generate a bed file containing the location of mimps.

#Set up and import correct modules.

from datetime import datetime
startTime = datetime.now()
import re, sys

from Bio import SeqIO
# from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

############################################
#Establish input file (FASTA) and start time.

infile = sys.argv[1]

print(f'\n\nRunning Mimp_finder.py on {infile} at {startTime}.')

#Open a .bed file using the name of the infile. The loctation of the mimp finder will be desposited within.

with open(f'{infile}_mimp_regex_hits.bed'.format(), "a") as file_object:

#Parse and loop through the FASTA input file. Search through each scaffold/contig of the FASTA looking for sequences which match the regex/mimp TIR "CAGTGGG..GCAA[TA]AA" (re.IGNORECASE remove case sensitivity issues).
#mimp_length = Create a variable for mimp length - Mimps cannot exceed 400 nucleotides.
#mimp_count = Create a variable to count the number of mimps found on each contid/scaffold.

    for sequence in SeqIO.parse(infile, "fasta"):
        print("="*80 + "\n")
        print(f'Searching {sequence.description} for mimps.\n')
        hit = re.finditer(r"CAGTGGG..GCAA[TA]AA", str(sequence.seq), re.IGNORECASE)
        mimp_length = 400
        mimp_reverse_length = -400
        mimp_count = 0

#For every instance where the mimp TIR "CAGTGGG..GCAA[TA]AA" is found, print the location of the hit and search that contig for the other mimp TIR.
#h_start = record the start location of the mimp

        for h in hit:
            if h:
                h_start = h.start()
                print(f'--Mimp TIR found: "{h.group()}" at position {h.start()} to {h.end()}--')
                hit_rc = re.finditer(r"TT[TA]TTGC..CCCACTG", str(sequence.seq), re.IGNORECASE)

#Ensure that the distance between the start of the first mimp TIR and end of the second mimp TIR is not greater than mimp_length (400 nucleotides).
#If TRUE, add 1 to the mimp_count variable and print the contig/scaffold description along with the start location of the mimp and the end location of the mimp in the .bed file opened (ensure fields are \t seperated to maintain .bed file convention).
#If the distance between the start of the first mimp TIR and the end of the second mimp TIR is greater than 400 nucleotides, it is not a mimp.
#If the second mimp TIR comes before the first mimp TIR found (second mimp TIR end - first mimp TIR start < 0), it is not a mimp.

                for h_rc in hit_rc:
                    print('Looking for reverse complement sequence(s)...')
                    h_rc_end = h_rc.end()
                    if h_rc:
                        print(f'Mimp reverse complement TIR found: "{h_rc.group()}" at position {h_rc.start()} to {h_rc.end()}.  Is this complete mimp?')
                        length = h_rc_end - h_start
                        if length > mimp_reverse_length:
                            if length < mimp_length:
                                print(f'Yes! Full mimp found at position {h.start()} to {h_rc.end()}.')
                                print("The mimp length is:")
                                print(f'{length}bp\n')
                                mimp_count = mimp_count+1
                                print(sequence.description, h.start(), h_rc.end(), "+", file=file_object, sep="\t")
                            else:
                                print("No, TIRs are at a greater distance than 400bp. \nOnly a partial mimp.\n")
                        else:
                            print("No, TIRs are at a greater distance than 400bp")

#Print the number of mimps found on that contig/scaffold.
#When whole FASTA has been looped through, print the file searched and time finished.


                print(f'\n{sequence.description} contains {mimp_count} mimp(s)')


    print("="*80 + "\n")

    print(f'\nMimp_finder.py complete. \n\nFile searched: {infile}\n\nFinished at: {startTime}')
	
