# Annotations

For annotations, the MAKER pipeline (version 3.01.04) was used. The Repeat Masking step was skipped, as we have already modeled and masked repeats in these assemblies.

We do not have transcriptome data for these isolates and, as we are not clear on their species, have used a *Fusarium* wide reference proteome for annotation.

## TOC

- [DirStructure](#directory-structure--set-up)
  - [Symlinks](#set-up-a-symbolic-link-symlink)
- [Maker](#maker)
  - [Maker Inputs](#maker-inputs)
  - [Running Maker](#running-maker)
  - [Processing Outputs](#processing-outputs)

## Directory Structure & Set up

First, I generated directories for each isolate in the `GenomeAnnotations` directory.

As no RNA sequencing data were available for these isolates, a reference proteome set was generated using the NCBI `RefSeq nr` database, using the search term ``Fusarium AND srcdb_refseq[PROP]`` with the protein option selected. The output was saved as `Fusarium_RefSeq-nr-fullProtdb.fasta` in the `AnnotationsRef` directory.

### Set up a symbolic link (symlink)

I then created a [symlink](https://www.futurelearn.com/info/courses/linux-for-bioinformatics/0/steps/201767) to the TNAU masked FASTA files, to save disk space.

``` bash
ln -s ../../GenomeAssemblies/S32/FASTA_versions/S32_V5-FS66-Contiglabelled.FullMask.fasta ./
```

## MAKER

I followed [this tutorial](https://weatherby.genetics.utah.edu/MAKER/wiki/index.php/MAKER_Tutorial_for_GMOD_Online_Training_2014) for the majority of the process.

### Maker inputs

Maker requires at least three input `.ctl` files. These files were generated using the `maker -CTL` command.

I edited the `maker_opts.ctl` file as required. The default output for the other files was left as is. Below I show the lines in the `maker_opts.ctl` file for the S16 assembly I edited.

```bash
#-----Genome (these are always required)
genome=/home/u1983390/Projects/NewTools-Proj/IndianGenomeAssemblies/GenomeAnnotations/S16/S16_V4-Contiglabelled.FullMask.fasta #genome sequence (fasta file or fasta embeded in GFF3 file)

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=/home/u1983390/Projects/NewTools-Proj/IndianGenomeAssemblies/GenomeAnnotations/AnnotationsRef/Fusarium_RefSeq-nr-fullProtdb.fasta  #protein sequence file in fasta format (i.e. from mutiple organisms)

#-----Gene Prediction
augustus_species=fusarium #Augustus gene prediction species model
protein2genome=1 #infer predictions from protein homology, 1 = yes, 0 = no
```

### Running Maker

Once the `.ctl` files had been edited, I ran maker using the `MakerCommand.sh`. The original script is saved in the `bin` directory, but the file was symlinked in each isolate directory. The `MakerCommand.sh` script contains the following;

```bash
 #!/bin/bash 

#Maker Command | Jamie Pike

#Run like this: nohup MakerCommand.sh 1>MakerRun_X.log &

export LIBDIR=/home/u1983390/apps/RepeatMaskerSetup/RepeatMasker/Libraries/ 
export REPEATMASKER_LIB_DIR=/home/u1983390/apps/RepeatMaskerSetup/RepeatMasker/Libraries/
export REPEATMASKER_MATRICES_DIR=/home/u1983390/apps/RepeatMaskerSetup/RepeatMasker/Libraries/Matrices

set -e #Set escape in case the maker fails.

python -c "print('=' * 75)"
echo "Maker Run" $(date)
python -c "print('=' * 75)"

/home/u1983390/apps/maker/bin/maker -RM_off maker_bopts.ctl maker_exe.ctl maker_opts.ctl 

python -c "print('=' * 75)"


SendEmail.py "Maker Annotation Complete" "Hi Jamie, the Maker annotation is now complete. :)"
```

---

For the S16 isolate, I ran `MakerCommand.sh` as follows:

```bash

#Ensure that i am using the AnnotationsEnv
nohup MakerCommand.sh 1>S16_Maker_01_141123.log &
```

There were some issues with the `augustus` executable stored on Vettel, I therefore also had to edit the `Maker_exe.ctl`, pointing to anther version of `augustus`.

```bash
augustus=/home/u1983390/miniconda3/envs/MaeiEnv/bin/augustus #location of augustus executable
```

#### Time logs

S16:

I set the maker annotation for S16 running at 17:06 on the 14th of November 2023.

When I checked at 08:55 on 20th November, MAKER had finished. I had checked on 15th, 16th, and 17th and MAKER was still running.

---

S32:

```bash
nohup MakerCommand.sh 1>S32_Maker_01_211123.log &
```

Here: `/home/u1983390/Projects/NewTools-Proj/IndianGenomeAssemblies/GenomeAnnotations/S32`

Start: 16:54, 21/11/2023
Finish: 22:21, 24/11/2023

---

S6:

```bash
nohup MakerCommand.sh 1>S6_Maker_01_271123.log &
```

Start: 10:37, 27/11/2023
Finish: 09:42, 01/12/2023

### Processing Outputs

Summary of outputs:

- S6: all contigs were successfully annotated.
- S16: all contigs were successfully annotated.
- S32: all contigs were successfully annotated.
- SY-2:
- S32_2:
- S56:

I therefore merged the output FASTA and GFF files, using the maker inbuilt `fasta_merge` and `gff_merge`.

```bash
# example for S16

fasta_merge -d S16_V4-Contiglabelled.FullMask_master_datastore_index.log # -g ensure we take only the Maker outputs, and not BLAST alignments etc.

gff3_merge -g -d S16_V4-Contiglabelled.FullMask_master_datastore_index.log
```

I then tried to gauge the quality of the annotations using the following steps.

First, I checked the number of gene models for, and compared it to similarly sized Fusarium assemblies.

```bash
#check the number of predicted genes
awk '{ if ($3 == "gene") print $0 }' S16_V4-Contiglabelled.FullMask.all.gff | awk '{ sum += ($5 - $4) } END { print NR, sum / NR }'
```

Results

- S6: 17891 1418.74
- S16: 15727 1521.13
- S32: 15824 1491.34
- S32_2
- SY-2
- S56

I then checked the AED score output using the `AED_cdf_generator.pl` script supplied with MAKER, and used [Sadik's script](/Volumes/Jamie_EXT/Projects/NewToolsProject/exp/GenomeAnnotations/bin/SadikShared) to generate a figure.

```bash
# use AED generator.
perl AED_cdf_generator.pl -b 0.010 S16_V4-Contiglabelled.FullMask.all.gff > S16_V4-Contiglabelled.FullMask.AED-0.010.txt
```

I also ran BUSCO on the `S16_V4-Contiglabelled.FullMask.all.maker.transcripts.fasta`, using the following command:

```bash
# BUSCO command for transcripts fasta

nohup busco -i  S16_V4-Contiglabelled.FullMask.all.maker.transcripts.fasta -o busco_Hypocreales -l hypocreales -m transcriptome -c 1 -f 1>BuscoOfTranscripts.log &
```

### BUSCO results for isolate predicted transcripts.

```bash
#  S16 transcripts:
 --------------------------------------------------
 |Results from dataset hypocreales_odb10           |
 --------------------------------------------------
 |C:98.7%[S:98.5%,D:0.2%],F:0.6%,M:0.7%,n:4494     |
 |4434 Complete BUSCOs (C)                       |
 |4426 Complete and single-copy BUSCOs (S)       |
 |8 Complete and duplicated BUSCOs (D)        |
 |26 Fragmented BUSCOs (F)                     |
 |34 Missing BUSCOs (M)                        |
 |4494 Total BUSCO groups searched               |
 --------------------------------------------------

# S32 transcripts:

 --------------------------------------------------
 |Results from dataset hypocreales_odb10           |
 --------------------------------------------------
 |C:98.7%[S:98.5%,D:0.2%],F:0.5%,M:0.8%,n:4494     |
 |4437 Complete BUSCOs (C)                       |
 |4427 Complete and single-copy BUSCOs (S)       |
 |10 Complete and duplicated BUSCOs (D)        |
 |23 Fragmented BUSCOs (F)                     |
 |34 Missing BUSCOs (M)                        |
 |4494 Total BUSCO groups searched               |
 --------------------------------------------------

# S6 transcripts:

 --------------------------------------------------
 |Results from dataset hypocreales_odb10           |
 --------------------------------------------------
 |C:99.1%[S:98.6%,D:0.5%],F:0.3%,M:0.6%,n:4494     |
 |4454 Complete BUSCOs (C)                       |
 |4432 Complete and single-copy BUSCOs (S)       |
 |22 Complete and duplicated BUSCOs (D)        |
 |13 Fragmented BUSCOs (F)                     |
 |27 Missing BUSCOs (M)                        |
 |4494 Total BUSCO groups searched               |
 --------------------------------------------------

```

I also visualised each MAKER2 gff3 on the genome using the genome browser IGV, to ensure the output's look sensible.
