# Prediction of candidate effectors

The gene models produced from the MAKER annotations were used to identify candidate effectors in the TNAU datasets. Predicted genes were filtered based on size (>30aa and <450aa), then submitted to SignalP (v4.1), sequences which were predicted to contain a signal peptide were parsed to EffectorP (v2.0.1).

## Directory set up

I created the `Effectors` directory in the `NewTools-project/exp/` directory. A subdirectory was then generated to store symlinks the annotation files, and avoid repeatedly typing lengthy paths.

```bash
# make the effectors directories 
mkdir -p Effectors Effectors/data

# create symlink for annotations. 

ln -s ../../GenomeAnnotations/S6/S6_V4.2-FocR1.Contiglabelled.FullMask.maker.output/S6_V4.2-FocR1.Contiglabelled.FullMask.all.maker.proteins.fasta ./
ln -s ../../GenomeAnnotations/S16/S16_V4-Contiglabelled.FullMask.maker.output/S16_V4-Contiglabelled.FullMask.all.maker.proteins.fasta ./
ln -s ../../GenomeAnnotations/S32/S32_V5-FS66-Contiglabelled.FullMask.maker.output/S32_V5-FS66-Contiglabelled.FullMask.all.maker.proteins.fasta ./
ln -s ../../GenomeAnnotations/SY-2/SY-2_V5-RepeatMasked.maker.output/SY-2_V5-RepeatMasked.all.maker.proteins.fasta ./
```

## Size filtering

Once the annotation files had been linked, i reduced the total number of candidates from each by filtering the effectors based on size.

```bash
# generate the processing directory.
mkdir sizeFilter

# generate index file for size filtering.
cd data 
for i in * ; do echo $i ; samtools faidx ${i} ; done 

# filter all sequences by size and create sorted bed file of short sequences . <450aa and >30aa accepted. 
for i in *.fai ; do echo "${i}..." ; awk '{if($2 < 450 && $2 > 30) print $1 "\t0\t" $2 "\t"}' ${i} > ../sizeFilter/${i}.sizeFilter.bed ; bedtools sort -i ../sizeFilter/${i}.sizeFilter.bed > ../sizeFilter/${i}.sizeFilter.sorted.bed ; echo "done." ; done

# generate fasta of filtered sequences
#---
# S6
bedtools getfasta -s -fi S6_V4.2-FocR1.Contiglabelled.FullMask.all.maker.proteins.fasta -bed ../sizeFilter/S6_V4.2-FocR1.Contiglabelled.FullMask.all.maker.proteins.fasta.fai.sizeFilter.sorted.bed -fo ../sizeFilter/S6_V4.2-FocR1.Contiglabelled.FullMask.all.maker.proteins.sizeFilter.sorted.fasta

# S16 
bedtools getfasta -s -fi S16_V4-Contiglabelled.FullMask.all.maker.proteins.fasta -bed ../sizeFilter/S16_V4-Contiglabelled.FullMask.all.maker.proteins.fasta.fai.sizeFilter.sorted.bed -fo ../sizeFilter/S16_V4-Contiglabelled.FullMask.all.maker.proteins.sizeFilter.fasta

# S32
bedtools getfasta -s -fi S32_V5-FS66-Contiglabelled.FullMask.all.maker.proteins.fasta -bed ../sizeFilter/S32_V5-FS66-Contiglabelled.FullMask.all.maker.proteins.fasta.fai.sizeFilter.sorted.bed -fo ../sizeFilter/S32_V5-FS66-Contiglabelled.FullMask.all.maker.proteins.sizeFilter.fasta

# SY-2
bedtools getfasta -s -fi  SY-2_V5-RepeatMasked.all.maker.proteins.fasta -bed ../sizeFilter/SY-2_V5-RepeatMasked.all.maker.proteins.fasta.fai.sizeFilter.sorted.bed -fo ../sizeFilter/SY-2_V5-RepeatMasked.all.maker.proteins.sizeFilter.fasta

cd ../sizeFilter/

# check the number of short sequences for each isolate
grep -c ">" *.fasta
S16_V4-Contiglabelled.FullMask.all.maker.proteins.sizeFilter.fasta:9148
S32_V5-FS66-Contiglabelled.FullMask.all.maker.proteins.sizeFilter.fasta:9253
S6_V4.2-FocR1.Contiglabelled.FullMask.all.maker.proteins.sizeFilter.sorted.fasta:10946
SY-2_V5-RepeatMasked.all.maker.proteins.sizeFilter.fasta:9115
```

## Signalp

I created a directory for the analysis, `mkdir signalp`. The size filtered fasta files were then symlinked in the `signalp` directory.

I  filtered the short sequences through signalp (v5) using the following shell script `signalp.sh`.

```bash
#!/bin/bash

for i in *.fasta ; do 
    echo "Processing ${i}..." ;  
    signalp -fasta ./${i} 1>./${i}-signalp.out 2>./${i}-signalp.log #Identify sequences with signal peptide.
    echo "Generating index files..." ;
    awk '$2 ~ /SP\(Sec\/SPI\)/ {print $1}' ./${i}/${i}-Translated_summary.signalp5 > ./${i}-signalPepPositives.list;
    echo "Generating fasta...";
    for j in $(cat ./${i}-signalPepPositives.list); do 
        samtools faidx ./${i} ${j}; 
        done > ./${i}-signalp.filtered.fasta;
    echo "done." ;
done
```

The shell script was run from the command line as follows,

```bash
# run using this command
../../bin/signalp.sh 2>1& | tee signalp.log
```

## EffectorP

Likely effectors were then predicted from this small secreted protein set using EffectorP (v2.0.1), run in the `signalp` directory using the following command:

```bash
# make directory for effectorp output
mkdir ../effectorp

# run effectorp
for i in *-signalp.filtered.fasta ; do \
    echo "processing ${i}..." ; \
    EffectorP.py -i ./${i} -E ../effectorp/${i}-EffectorP.filtered.fasta > ../effectorp/${i}-EffectorP.filtered.log ; \
    echo "done." ; \
done 
```

I checked the total number of candidates per assembly using grep.

```bash
cd ../effectorp/
grep -c  ">" *.fasta
S16_V4-Contiglabelled.FullMask.all.maker.proteins.sizeFilter.fasta-signalp.filtered.fasta-EffectorP.filtered.fasta:289
S32_V5-FS66-Contiglabelled.FullMask.all.maker.proteins.sizeFilter.fasta-signalp.filtered.fasta-EffectorP.filtered.fasta:314
S6_V4.2-FocR1.Contiglabelled.FullMask.all.maker.proteins.sizeFilter.sorted.fasta-signalp.filtered.fasta-EffectorP.filtered.fasta:333
SY-2_V5-RepeatMasked.all.maker.proteins.sizeFilter.fasta-signalp.filtered.fasta-EffectorP.filtered.fasta:289
```

As the file names have become very large, I renamed the fastas. 

```bash
# S6
mv S6_V4.2-FocR1.Contiglabelled.FullMask.all.maker.proteins.sizeFilter.sorted.fasta-signalp.filtered.fasta-EffectorP.filtered.fasta S6_V4.2-candidateEffectors.fasta

# S16
mv S16_V4-Contiglabelled.FullMask.all.maker.proteins.sizeFilter.fasta-signalp.filtered.fasta-EffectorP.filtered.fasta S16_V4-candidateEffectors.fasta

# S32
mv S32_V5-FS66-Contiglabelled.FullMask.all.maker.proteins.sizeFilter.fasta-signalp.filtered.fasta-EffectorP.filtered.fasta S32_V5-FS66-candidateEffectors.fasta

# SY-2
mv SY-2_V5-RepeatMasked.all.maker.proteins.sizeFilter.fasta-signalp.filtered.fasta-EffectorP.filtered.fasta SY-2_V5-candidateEffectors.fasta
```

## Cd-hit clustering

I then clustered the candidate effectors at 80% identity to identify shared candidates. Initially, i had to create a single fasta.

```bash
# create an empty file to hold all of the candidate effectors.
touch ./AllCandidateEffectorSets.fasta 
echo "Clustering final effector sets..."
for i in *.fasta; 
# Add genome filename to start of Fasta Headers so we know which isolate this came from and  combine all of the individual candidate effector sets by adding them to the empty candidate effector fasta.
  do
  awk '/>/{sub(">","&"FILENAME"_");sub(/\.fasta/,x)}1' ${i} >> AllCandidateEffectorSets.fasta ; 
done 
```

Then cluster the effectors by sequence identity (80%).

```bash
# make a directory for the cdhit output
mkdir cdhit

# symlink the cdhit input file (from the effectorp dir).
ln -s AllCandidateEffectorSets.fasta ../cdhit/

cdhit -i ../effectorp/AllCandidateEffectorSets.fasta -d 0 -o ./AllCandidateEffectorSets -c 0.80 -n 5  -G 1 -g 0 -b 20 -l 10 -s 0.0 -aL 0.0 -aS 0.0 1> cd-hit.log
```

## Orthofinder MCL

Finding orthoglous groups of the candidate effectors from each TNAU assembly.

### Build Venn Diagram

I processed the `cd-hit` data using a custom python script, [ProcessingCdhit.py](https://github.com/JamiePike/NewTools-Project/blob/master/bin/ProcessingCdhit.py).

I renamed the cdhit output fasta too `mv AllCandidateEffectorSets AllCandidateEffectorSets_cdhit.fasta`

I then tried to identify some of the interesting groups or candidates by plotting as a Venn diagram, using the [`CandidateEffectorsVenn.R`](https://github.com/JamiePike/NewTools-Project/blob/master/bin/CandidateEffectorsVenn.R) script.

---

Going forward...

- BLAST against focub genomes?
- BLAST focub set against these assemblies?
