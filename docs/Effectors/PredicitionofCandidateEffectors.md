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
