#!/home/u1983390/miniconda3/envs/biopythonEnv/bin/python

import pandas as pd
import numpy as np
import sys

infile = sys.argv[1] #Infile argument for the command line. 

column_names=["Cluster","header"] #Name the two columns pandas will detect. 

clstr= pd.read_csv(infile, sep='\t', names=column_names) #Read in the .clstr file. 

clstr.loc[clstr['Cluster'].str.isnumeric(), 'Cluster'] = np.nan # Fill each section in the cluster column that does not have the ">Cluster *" string (i.e., the number in the cluster), with a null vlaue.
clstr['Cluster'].fillna(method = 'ffill', inplace=True) #Replace all null values with the string in the pervious cell. 
clstr.dropna(inplace=True) #Drop the rows with remaining NaN values (the orgignal rows with just the >Cluster * name)
clstr['Cluster'] = clstr['Cluster'].map(lambda x: x.lstrip('>')) #Remove ">" from the Cluster column 
clstr[['Size (aa)','Genome', 'Name', 'Identity']] = clstr['header'].str.split(', >|-candidateEffectors_| ...', expand=True) #Split the header column in to size and sequence Name
clstr.drop(columns='header', inplace=True) #Drop the header column, as we have now split this column 
clstr['Name'] = clstr['Name'].str.split('\.\.\.').str[:-1].str.join('\.\.\.') #Remove the additional characters added by cd-hit.
clstr['Size (aa)'] = clstr['Size (aa)'].str.split('aa').str[:-1].str.join('aa') #Remove "aa," added to the aa size by cdhit.


occur = clstr.groupby(['Cluster', 'Genome'], as_index=False)["Name"].count() #Count the number of candNameate effectors from each genome in each cluster. 
pivot = occur.pivot(index='Cluster', columns='Genome', values='Name').fillna("0") #Reformat the dataframe as a matrix which can be used to build a heatmap.
clstr['Cluster size'] = clstr.groupby('Cluster')['Cluster'].transform('count') #Calculate the total number of sequences in each cluster
clstr['Total isolate seqs in cluster'] = clstr.groupby(['Cluster', 'Genome'], as_index=False)["Name"].transform('count') #Calculate the number of sequences from each genome in each cluster
clstr['Isolate percentage in Cluster'] = clstr['Total isolate seqs in cluster']/clstr['Cluster size'] * 100 #Calculate the percentage of seqs in the cluster from each genome. 

clstr = clstr[['Cluster', 'Cluster size', 'Genome', 'Total isolate seqs in cluster', 'Name', 'Identity', 'Isolate percentage in Cluster']]


clstr.to_csv(infile+"_table.tsv", header=True, index=False, sep="\t") #Write out the table of clusters, genomes, pNameetns and candNameate effectors 
pivot.to_csv(infile+"_heatmapdata.tsv", header=True, index=True, sep="\t") #Write out the matrix of occurances per genome. 