# Phylogeny of TNAU isolates

TEF1-a and RPBii, common genetic barcodes in *Fusarium*, were used to generate phylogenies for the TNAU genomes.

Databases for TEF1-a and RPBii were compiled from NCBI reference sequences as part of a Undergraduate project with Dr John Clarkson ([see TEF1-a fasta](../../data/phylogenies/Tef1a_db.fasta)).

The Fusarium assemblies were symlinked from the [data directory](../../data/genomes/All_Fusarium_Genomes/)

## BLASTN Search

Using BLASTN (v2.9.0+; `1e-6` cut-off), homologs were identified in our Fusarium assembly database and the TNAU genomes using the `BigBlast.sh` script ([see bin](../../bin/BigBlast.sh)). Hits with >70% identity and >90% coverage were manually extracted via Samtools (v1.15.1).

```bash
# tef
mkdir TEF # create dir
cd TEF # move into it
cp ../../data/* ./ #copy data
BigBlast.sh -q Tef1a_db.fasta -g FastaList.txt -i 70 -c 90 -t blastn -b #run search
cd .. # move out of dir

# rbp2
mkdir RPBii # create dir
cd RPBii # move into it 
cp ../../data/* ./ #copy data 
BigBlast.sh -q Ena-RPB2.fasta  -g FastaList.txt -i 70 -c 90 -t blastn -b #run search
```

To increase space, I removed the copied fasta's and created a symlink to the originals.

```bash
#remove copies
for i in $(cat FastaList.txt); do rm ${i} ; done
#create link
ln -s ../../data/* ./
```

Top scoring TEF1-a and RPBii regions were extracted from each genome manually and a single FASTA created for each sequence. These FASTAs were concatenated to generate a TEF1-a and RPBii multiFASTA files.

```bash
#example for manual extraction
grep -v "partial" Ena-RPB2_vs_F._oxysporum_Fo47.blastn.outfmt6.evalue_used_1e-6.out | column -t #check the top scoring hits for full length sequences and ensure they are in a similar location and full length query.

samtools faidx F._oxysporum_Fo47.fna # index the fasta

samtools faidx F._oxysporum_Fo47.fna CP052042.1_Fusarium_oxysporum_Fo47_chromosome_V:3933350-3935771 >F._oxysporum_Fo47_RPB2.fna # extract the top scoring full length hit and save to new file. 
```

```bash
# I looped through each directory, as they all started with an "F" and the file ending for each extracted barcode was the same. 

# tef
cd ../../TEF
cd BigBlast*
for i in F* ; do awk '/>/{sub(">","&"FILENAME"_");sub(/\.fna/,x)}1' ${i}/*_TEF.fna >> Fusarium_TEF.fna ; done 

#rbp2
cd ../../RPBii
cd BigBlast*
for i in F* ; do awk '/>/{sub(">","&"FILENAME"_");sub(/\.fna/,x)}1' ${i}/*_RPB2.fna >> Fusarium_RPB2.fna ; done
``B

## MAFFT alignment

Before alignment, the barcodes sequences from the novel species [F. mindanaoense](https://www.mdpi.com/2309-608X/9/4/443) were added.

MAFFT (v7.505) was used align the sequences in the concatenated FASTA, with the command in the shell script `MafftCommand.sh`

```bash
# tef
MafftCommand.sh Fusarium_TEF.fna

# rbp2
MafftCommand.sh Fusarium_RPB2.fna
```

Overhanging sequences were trimmed manually, and the trimmed output saved for each barcode into a new FASTA files. Some of the isolate names were altered to shorten them, as MAFFT can sometimes trim off longer names.

## IQtree2 Analysis

IQ-TREE was used to (v2.2.0.3) infer maximum-likelihood phylogenies with 1000 ultrafast bootstrap replicates. The `Iqtree2Command.sh` shell script was used for this also.

```bash
#tef
Iqtree2Command.sh Fusarium_TEF_Additional_Genomes_MAFFT.msa Fusarium_Additional_Species_TEF

#rbp2
Iqtree2Command.sh  RPBii-Trimmed.fasta RPBii-iqtree2Phylogeny
```

Phylogenies were then visualised and annotated using [iTOL](https://itol.embl.de).
