# Phylogeny of TNAU isolates

TEF1-a and RBPii, common genetic barcodes in *Fusarium*, were used to generate phylogenies for the TNAU genomes.

Databases for TEF1-a and RBPii were compiled from NCBI reference sequences as part of a Undergraduate project with Dr John Clarkson ([see TEF1-a fasta](../../data/phylogenies/Tef1a_db.fasta)).

The Fusarium assemblies were symlinked from the [data directory](../../data/genomes/All_Fusarium_Genomes/)

## BLASTN Search

Using BLASTN (v2.9.0+; `1e-6` cut-off), homologs were identified in our Fusarium assembly database and the TNAU genomes using the `BigBlast.sh` script ([see bin](../../bin/BigBlast.sh)). Hits with >70% identity and >90% coverage were manually extracted via Samtools (v1.15.1).

```bash

# tef
BigBlast.sh -q ../../data/phylogenies/Tef1a_db.fasta -g FastaList.txt -i 70 -c 90 -blastn -b 

# rbp2
BigBlast.sh -q ../../data/phylogenies/Ena-RPB2.fasta  -g FastaList.txt -i 70 -c 90 -blastn -b 
```

TEF1-a and RBPii regions were extracted from each genome and a single FASTA created for each sequence. These FASTAs were concatenated to generate a TEF1-a and RBPii multiFASTA files.

```bash
# I looped through each directory, as they all started with an "F" and the file ending for each extracted barcode was the same. 

# tef
for i in F* ; do cat ${i}/*._TEF.fna >> Fusarium_TEF.fna ; done 

#rbp2
for i in F* ; do cat ${i}/*._RBP2.fna >> Fusarium_RBP2.fna ; done 
```

## MAFFT alignment

MAFFT (v7.505) was used align the sequences in the concatenated FASTA, with the command in the shell script `MafftCommand.sh`

```bash
# tef
MafftCommand.sh Fusarium_TEF.fna

# rbp2
MafftCommand.sh Fusarium_RBP2.fna
```

Overhanging sequences were trimmed manually, and the trimmed output saved for each barcode into a new FASTA files. Some of the isolate names were altered to shorten them, as MAFFT can sometimes trim off longer names.

## IQtree2 Analysis

IQ-TREE was used to (v2.2.0.3) infer maximum-likelihood phylogenies with 1000 ultrafast bootstrap replicates. The `Iqtree2Command.sh` shell script was used for this also.

```bash
#tef
Iqtree2Command.sh Fusarium_TEF_Additional_Genomes_MAFFT.msa Fusarium_Additional_Species_TEF

#rbp2
Iqtree2Command.sh  Fusarium_RPB2-MAFFTaln.Phylogeny2.msa Fusarium_RPB2-MAFFTaln.Phylogeny2.msa
```

Phylogenies were then visualised and annotated using [iTOL](https://itol.embl.de)
