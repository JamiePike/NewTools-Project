 #!/bin/bash

Assemb=${1?Please provide the assembly input fasta.} #Input Assembly.
Bam=${2?Please provide a bam file of raw reads mapped back to the assembly} #Input BAM.
prefix=${3?Please provide a prefix for you outputs.}
python -c "print('=' * 75)"
echo "Blobtools BLAST"
echo "---------------"
echo $(date)
python -c "print('=' * 75)"
echo "Blasting..."
#Perform BLASTN for Blobtools.
blastn -task megablast -query ${Assemb} -db /shared/reference/NCBI_NT/BLASTDB/nt  -outfmt '6 qseqid staxids bitscore std' -max_target_seqs 1 -max_hsps 1 -num_threads 16 -evalue 1e-25 -out ${prefix}.vs.nt.mts1.hsp1.1e25.megablast.out 

python -c "print('=' * 75)"
echo "Running Blobtools..."

#Perform the Blobtools Analysis
/home/u1983390/apps/blobtools/blobtools create -i ${Assemb} -b ${Bam}  -t ${prefix}.vs.nt.mts1.hsp1.1e25.megablast.out  -o ${prefix}.blobtools

echo "Generating Blobtools figures..."

#To view results.
/home/u1983390/apps/blobtools/blobtools view -i ${prefix}.blobtools.blobDB.json -r species
/home/u1983390/apps/blobtools/blobtools view -i ${prefix}.blobtools.blobDB.json -r phylum

#To generate graphs.
/home/u1983390/apps/blobtools/blobtools plot -i ${prefix}.blobtools.blobDB.json -r species
/home/u1983390/apps/blobtools/blobtools plot -i ${prefix}.blobtools.blobDB.json -r phylum

#To generate PDF versions.
/home/u1983390/apps/blobtools/blobtools plot -i ${prefix}.blobtools.blobDB.json -r species --format pdf
/home/u1983390/apps/blobtools/blobtools plot -i ${prefix}.blobtools.blobDB.json -r phylum --format pdf

python -c "print('=' * 75)"
#Send email notifciation of job completion.
SendEmail.py "Blobtools Complete" "Blobtools analysis now complete. Check for results."
echo $(date)
python -c "print('=' * 75)"      
