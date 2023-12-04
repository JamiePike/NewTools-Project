#!/bin/bash

#Generate a fastq of mapped reads and unmapped reads. 

bam=${1?Please provide a BAM file input.}
prefix=${2?Please provide an output file prefix}

python -c "print('=' * 75)"
echo "Fastq Generator"
echo "---------------"
echo $(date)
python -c "print('=' * 75)"

echo "Extracting..."
#Extract the mapped reads
samtools view -h -b -f 1 -F 12 ${bam} > ${prefix}_MappedPairs.bam 
# R1 unmapped, R2 mapped
samtools view -h -b -f 4 -F 264 ${bam}> ${prefix}_MappedRight.tmp.bam
# R1 mapped, R2 unmapped
samtools view -h -b -f 8 -F 260 ${bam}> ${prefix}_MappedLeft.tmp.bam
# R1 & R2 unmapped
samtools view -h -b -f 12 -F 256 ${bam}> ${prefix}_UnmappedPairs.tmp.bam
echo "Merging..."
#Now merge all of the unmapped reads into one bam 
samtools merge -u ${prefix}_AllUnmappedReads.bam ${prefix}_MappedRight.tmp.bam ${prefix}_MappedLeft.tmp.bam ${prefix}_UnmappedPairs.tmp.bam
echo "Sorting..."
#Next, these BAM files must be resorted so that they are ordered by read ID instead of location in the reference.
samtools sort -n ${prefix}_MappedPairs.bam ${prefix}_MappedPairs.sorted
samtools sort -n ${prefix}_AllUnmappedReads.bam ${prefix}_AllUnmappedReads.sorted

#Validate that the correct total number of reads have been extracted and sorted. 
python -c "print('-' * 50)"
echo "It is a good idea to check that you have the correct number of reads and no redundancy. Therefore, the orginal BAM file has been summarised and the read numbers from the mapped bam and unmapped bam provided."
echo ${bam} "Stats:"
samtools flagstat ${bam}
echo ${bam} "Stats:" 
samtools view -c ${prefix}_MappedPairs.sorted.bam
echo ${bam} "Stats:"
samtools view -c ${prefix}_AllUnmappedReads.sorted.bam 
python -c "print('-' * 50)"

#Continue to generate fastq files. 
echo "Generating fastqs..."
bamToFastq -i ${prefix}_MappedPairs.sorted.bam -fq ${prefix}_Mapped.1.fastq -fq2 ${prefix}_Mapped.2.fastq 
bamToFastq -i ${prefix}_AllUnmappedReads.sorted.bam -fq ${prefix}_Unmapped.1.fastq -fq2 ${prefix}_Unmapped.2.fastq 
echo "Tidying up..."
rm *.tmp.bam

python -c "print('=' * 75)"

#Send email notifciation of job completion.
SendEmail.py "Fastqs Generated" "Your fastq files for ${prefix} have now been generated."

python -c "print('=' * 75)" 
