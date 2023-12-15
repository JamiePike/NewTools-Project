#!usr/bin/envs/python3

#Jamie Pike - 23/5/2022
#Covert Maei output to matrix in csv format

###########################################################
#Set up and import correct modules.
import re, sys, os, glob
import numpy as np
import pandas as pd

from Bio import SeqIO
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
# from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
###########################################################
#Parse command line arguments
parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("-p", "--path", help="Path the the EffectorBlast.sh output directory., e.g., ~/EffectorSearch/EffectorBlast_13-05-2022/", required = True)
parser.add_argument("-i", "--input_fasta", help="Input the query sequence FASTA used for EffectorBlast.sh", required = True)
parser.add_argument("-o", "--output", help="output file name", required = True)
args = vars(parser.parse_args())


###########################################################
#Determine input files
path1 = args["path"]
path2 = path1 + "*/Occurrences_of*.txt"

files_names = []
for file in glob.glob(path2):
    files_names.append(os.path.relpath(file))

infile1 = args["input_fasta"]
output = args["output"]
###########################################################
#Define the list of sequences

#Create array to store sequences in
seqList=[]
#Loop through FASTA file to get compelete list of effectors
for sequence in SeqIO.parse(infile1, "fasta"):
    seqList.append(str(sequence.id))

###########################################################
#Create dataframes

#Create an empty dataframe for merging at the end
df = pd.DataFrame(columns=["sequence"])

#Loop through all of the input files.
for infile in files_names:
    #Create the dataframe from the sequence hits in the genomes (identified in the Occurrences file input file (infile)).
    Occurrences=pd.read_csv(infile, delimiter='\t', index_col=False)    #Read the input file as a tab separated dataframe.
    pd.set_option("display.max_colwidth", None) #Ensure that the sequence names are not cut off.
    Occurrences.rename(columns = {list(Occurrences)[0]: 'sequence'}, inplace = True) #Name the sequences column
    #Trim the file path so that the column names are just th genomes.
    InName = os.path.basename(infile)   #Trim to just the file name
    Trimmed = InName.replace('Occurrences_of_SIX_Lycopersici_sequences_in_','') #Trim the file name preffix
    ColName = Trimmed.replace('.txt','')    #Trim the file name suffix
    #Create the data frame that contains the hits and the sequences
    Occurrences.rename(columns = {list(Occurrences)[1]: ColName}, inplace = True) #Name the Occurrences column
    Occurrences['sequence'] = Occurrences['sequence'].str.strip()   #Remove any white space from the sequence column
    Occurrences = Occurrences.set_index('sequence').reindex(seqList, fill_value="0").reset_index()  #Reset indexes
    df = df.merge(Occurrences, on='sequence', how='outer') #merge the current data frame with the existing data frame.

#Create output csv.
df.to_csv(output+'.csv', index = False)
