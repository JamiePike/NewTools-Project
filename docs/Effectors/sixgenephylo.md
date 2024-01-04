# SIX gene prediction and phylogenies

## SIX gene identification

*SIX* genes were identified in the *Fusarium* assemblies using tBLASTn. A reference set of *SIX* genes using in [this paper] were downloaded from NCBI.

| **_SIX gene_** | **Genbank Reference** | **Isolate** |
|---|---|---|
| _SIX1_ | AJ608702.1 | Foly |
| _SIX2_ | GQ268949.1 | Foly_BFOL-51 |
| _SIX3_ | AM234063.1 | Foly |
| _SIX4_ | GQ268951.1 | Foly_BFOL-51 |
| _SIX5_ | GQ268952.1 | Foly_BFOL-51_5 |
| _SIX6_ | GQ268953.1 | Foly_BFOL-51 |
| _SIX7_ | GQ268954.1 | Foly_BFOL-51_7 |
| _SIX8_ | FJ755837.1 | Foly_007 |
| _SIX9_ | KC701447.1 | Foly_007 |
| _SIX10_ | KC701448.1 | Foly_007 |
| _SIX11_ | KC701449.1 | Foly_007 |
| _SIX12_ | KC701450.1 | Foly_007 |
| _SIX13_ | KC701451.1 | Foly_007 |
| _SIX14_ | KC701452.1 | Foly_007 |
| _SIX15_ | KY073750.1 | Foly_1943 |

```bash
BigBlast.sh -q SIX_lycopersici.fasta -g FastaList.txt -i 1 -c 1 -t tblastn -b 
```

Top scoring SIX gene regions were extracted from each genome manually and a FASTA created for each gene sequence. These FASTAs were concatenated to generate a multiFASTA files for each *SIX* gene.

```bash
# example of manual extraction
# move into dir.
cd F._oxysporum_f._sp._cepae_FoC_Fus2
# examine the output, including locations
cat SIX_lycopersici_vs_F._oxysporum_f._sp._cepae_FoC_Fus2.hits_within_threshold.bed 
# create fasta index.
samtools faidx F._oxysporum_f._sp._cepae_FoC_Fus2.fna 
# extract and view the sequence.
samtools faidx F._oxysporum_f._sp._cepae_FoC_Fus2.fna MRCU01000010.1_Fusarium_oxysporum_f._sp._cepae_strain_FoC_Fus2_contig_10,_whole_genome_shotgun_sequence:809491-809836 
# save the sequence into a fasta file.
samtools faidx F._oxysporum_f._sp._cepae_FoC_Fus2.fna MRCU01000010.1_Fusarium_oxysporum_f._sp._cepae_strain_FoC_Fus2_contig_10,_whole_genome_shotgun_sequence:809491-809836 > F._oxysporum_f._sp._cepae_FoC_Fus2_SIX9.fasta 
```

I then created a multifasta file containing all copies of each *SIX* gene from the assemblies. 

```bash 
# generate SIX1 fasta
for i in F* ; do awk '/>/{sub(">","&"FILENAME"_");sub(/\.fna/,x)}1' ${i}/*SIX1.fasta >> FusSIX1.fasta; done

# generate six2 fasta
for i in F* ; do awk '/>/{sub(">","&"FILENAME"_");sub(/\.fna/,x)}1' ${i}/*SIX2.fasta >> FusSIX2.fasta; done

# generate six4 fasta
for i in F* ; do awk '/>/{sub(">","&"FILENAME"_");sub(/\.fna/,x)}1' ${i}/*SIX4.fasta >> FusSIX4.fasta; done

# generate six6 fasta
for i in F* ; do awk '/>/{sub(">","&"FILENAME"_");sub(/\.fna/,x)}1' ${i}/*SIX6.fasta >> FusSIX6.fasta; done

# generate six9 fasta
for i in F* ; do awk '/>/{sub(">","&"FILENAME"_");sub(/\.fna/,x)}1' ${i}/*SIX9.fasta >> FusSIX9.fasta; done
```

I then aligned each fasta file using the `MafftCommand.sh` shell script.

```bash
# loop through the fasta and align each. 
for i in FusSIX* ; do echo Processing $i ; ../../../../../bin/MafftCommand.sh ${i} ; echo done. ; done
```

I then trimmed the alignments manually before the iqtree2 analysis. Further, these sequence headers are too long for Iqtree2, so i trimmed the fasta headers manually.

```bash
# loop through files for iqtree2
for i in FusSIX*.aln ; do echo Processing $i ; ../../../../../bin/Iqtree2Command.sh ${i} ${i}.phylo ; echo done. ; done
```

I tidied up the directory next by making a directory for each *SIX* gene. 

```bash 
# make dirs
mkdir -p phylos SIX1Phylo SIX2Phylo SIX4Phylo SIX6Phylo SIX9Phylo

# move files
mv FusSIX1.* ./SIX1Phylo/
mv FusSIX2.* ./SIX2Phylo/
mv FusSIX4.* ./SIX4Phylo/
mv FusSIX6.* ./SIX6Phylo/
mv FusSIX9.* ./SIX9Phylo/
mv *Phylo ./phylos/
mv phylos/ ../
```

I realised is trimmed the SIX4 file, and ran the iqtree2 loop on the untrimmed, so repeated in the new dir. 

```bash
#clear out dir of iqtree2 output 
rm FusSIX4.fasta.aln.phylo*

#rerun iqtree2
../../../../../../bin/Iqtree2Command.sh FusSIX4-trimmed.fasta.aln.fa FusSIX4-trimmed.phylo
```

I then visualised the trees using [iTol](https://itol.embl.de).

