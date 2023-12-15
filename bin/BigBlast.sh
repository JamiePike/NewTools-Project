#!/bin/bash

#Jamie Pike - Last developed 07/06/2022

#This script can be used to search a query FASTA file against a list of genome/subject files.
#The user must determine the type of BLAST search to be performed and can list the subject FASTAS
#in a .txt file. The aim of this script is to automate any repetative BLAST searches that must be
#undertaken for a number of subject sequences, producing not only the BLAST results in the standard
#output, as outfmt 6 and a bed file for each input FASTA, but also but can also filter results based
#on percentage identity and coverage, producing a bed file for hits within the threshold as well as
#a tsv indicating the number of times a query sequence has been idenfited (within that threshold).

#The output tsvs can be used, in tandem with the FASTA list .txt file, by my BuildTable.py script to
#create a data matrix indicating the sequences that have been identified within the threshold and
#within which subject FASTAs. This can then be used to generate a heatmap with various tools in R. 

function main(){
  check $@

  #Set Variables

  date=$(date '+%d-%m-%Y')
  echo "
  BigBlast_v2.sh
  Date: $date"
  query_filename="${query##*/}"
  query_extension="${query_filename##*.}"
  query_filename="${query_filename%.*}"


  red='\033[01;31m'
  none='\033[0m'

  #Make blastdbs
  #-----------

  mkdir BigBlast_${BLAST}_${date} #Make directories for the output.


  #Indicate the bedfile option was selected.
    if [ $b ]
    then
      echo "Bed file option selected." ;
      #Log the flags used and the files searched
      echo "
    Date of Search: ${date}

    Query Input File: ${query}
    Genome List File: ${genomes}
    BLAST Search Type: ${BLAST}
    Percentage Identity Threshold:  ${PercentageIdent}%
    Percentage Coverage Threshold:  ${PercentageCov}%

    Bed files generated.

    Subject FASTAs  searched:" > BigBlast_${BLAST}_${date}/BigBlastSummary.txt ;
    for i in $(cat ${genomes}) ; do echo "    "${i}; done >> BigBlast_${BLAST}_${date}/BigBlastSummary.txt ;
    else
      echo "No bed files being generated. Please use the -b flag if you would like a bed file for each genome search." ;
      #Log the flags used and the files searched
      echo "
      Date of Search: ${date}

      Query Input File: ${query}
      Genome List File: ${genomes}
      BLAST Search Type: ${BLAST}
      Percentage Identity Threshold:  ${PercentageIdent}%
      Percentage Coverage Threshold:  ${PercentageCov}%

      Bed files not generated.

      Subject FASTAs  searched:" > BigBlast_${BLAST}_${date}/BigBlastSummary.txt
      for i in $(cat ${genomes}) ; do echo "  "${i}; done >> BigBlast_${BLAST}_${date}/BigBlastSummary.txt
    fi

  #Indicate the genomes being searched
    echo "Databases searched:"
    for i in $(cat ${genomes}) ;
      do echo "   "${i} ;
    done
    python -c "print('=' * 80)" ;


  #Make Blast Databases.
  #--------------------
    if [ $BLAST = "tblastx" ]
    then
      #For the query sequences...
      echo "Making ${BLAST} database for ${query}..." ;
      makeblastdb -dbtype nucl -in ${query} ;
      python -c "print('=' * 80)" ;
      #For the subject sequences
      for i in $(cat ${genomes}) ;
        do
        echo "Making ${BLAST} database for ${i}..." ;
        makeblastdb -dbtype nucl -in ${i} ;
        python -c "print('=' * 80)" ;
        done
    elif [ $BLAST = "tblastn" ]
    then
      #For the query sequences...
      echo "Making tblastn database for ${query}..." ;
      makeblastdb -dbtype prot -in ${query} ;
      python -c "print('=' * 80)" ;
      #For the subject sequences
      for i in $(cat ${genomes}) ;
        do
        echo "Making tblastn database for ${i}..." ;
        makeblastdb -dbtype nucl -in ${i} ;
        python -c "print('=' * 80)" ;
        done
    elif [ $BLAST = "blastn" ]
    then
      #For the query sequences...
      echo "Making blastn database for ${query}..." ;
      makeblastdb -dbtype nucl -in ${query} ;
      python -c "print('=' * 80)" ;
      #For the subject sequences
      for i in $(cat ${genomes}) ;
        do
        echo "Making blastn database for ${i}..." ;
        makeblastdb -dbtype nucl -in ${i} ;
        python -c "print('=' * 80)" ;
        done
    elif [ $BLAST = "blastp" ]
    then
      #For the query sequences...
      echo "Making blastp database for ${query}..." ;
      makeblastdb -dbtype prot -in ${query} ;
      python -c "print('=' * 80)" ;
      #For the subject sequences
      for i in $(cat ${genomes}) ;
        do
        echo "Making blastp database for ${i}..." ;
        makeblastdb -dbtype prot -in ${i} ;
        python -c "print('=' * 80)" ;
        done
    fi

  #Run ${BLAST}
  #-----------

  echo "//Performing ${BLAST} search..."
  #Indicate the bedfile option was selected.
  if [ $b ]
  then
      for i in $(cat ${genomes}) ; #Loop through the genomes provided in the input txt file (-g flag).
        do
        python -c "print('=' * 80)" ;
        #Create output file names
        subject_filename="${i##*/}" ;
        subject_extension="${subject_filename##*.}" ;
        subject_filename="${subject_filename%.*}" ;
        #Create output directories
        mkdir ./BigBlast_${BLAST}_${date}/${subject_filename} ;
        mkdir ./BigBlast_${BLAST}_${date}/${subject_filename}/BlastDB ;
        #Print statements
        echo "Running ${BLAST} for ${i}..." ;
        echo "out: ./BigBlast_${BLAST}_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.${BLAST}.outfmt6.evalue_used_1e-6.out"
        #Generate generic BLAST output.
        ${BLAST} -query ${query} -db ${i} -evalue 1e-6  -out ./BigBlast_${BLAST}_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.${BLAST}.evalue_1e-6.summarised.out ;
        #Perform ${BLAST}; output in outfmt 6 including all usual colummns plus the query coverage column.
        ${BLAST} -query ${query} -db ${i}  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" -evalue 1e-6  -out ./BigBlast_${BLAST}_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.${BLAST}.outfmt6.evalue_used_1e-6.out ;
        #Check whether any hits were found. If found, proceed with the next step.
        if [[ -s ./BigBlast_${BLAST}_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.${BLAST}.outfmt6.evalue_used_1e-6.out ]] ;
        then
          #Extract the hits which meet the percentage identity and percentage coverage threshold and prepare a txt file that has both contigs (effector and contig in genome) and a bed file that contains all hits.
          awk -v p="$PercentageIdent" -v c="$PercentageCov" '{ if ($3 >= p && $13>= c) print $0}' ./BigBlast_${BLAST}_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.${BLAST}.outfmt6.evalue_used_1e-6.out | awk '{ if ($11 < 1e-6) print $1,"\t",$2}' > ./BigBlast_${BLAST}_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.query_and_subject_hits.out ;
          cut -f 1,2,9,10 ./BigBlast_${BLAST}_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.${BLAST}.outfmt6.evalue_used_1e-6.out | awk '{print $2,$3,$4,$1}' OFS='\t' | awk '{if ($2>$3)print $1,$3,$2,$4,".","-";else print $1,$2,$3,$4,".","+";}' OFS='\t' |awk '{a=$2-1;print $1,a,$3,$4,$5,$6;}' OFS='\t'| bedtools sort > ./BigBlast_${BLAST}_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.all_hits.bed ;
          #Prepare a bedfile that has all the hits within the threshold set.
          awk -v p="$PercentageIdent" -v c="$PercentageCov" '{ if ($3 >= p && $13>= c) print $0}' OFS='\t' ./BigBlast_${BLAST}_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.${BLAST}.outfmt6.evalue_used_1e-6.out | awk '{ if ($11 < 1e-6) print $0}' OFS='\t' | cut -f 1,2,9,10 | awk '{print $2,$3,$4,$1}' OFS='\t' | awk '{if ($2>$3)print $1,$3,$2,$4,".","-";else print $1,$2,$3,$4,".","+";}' OFS='\t' |awk '{a=$2-1;print $1,a,$3,$4,$5,$6;}' OFS='\t'| bedtools sort > ./BigBlast_${BLAST}_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.hits_within_threshold.bed
          #Create a tab separated table which contains the candidate effector hit and the number of hits for that effector.
          awk '{a[$1]++}END{for(i in a){print i,"\t",a[i]}}' ./BigBlast_${BLAST}_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.query_and_subject_hits.out > ./BigBlast_${BLAST}_${date}/${subject_filename}/Occurrences_of_${query_filename}_sequences_in_${subject_filename}.txt
          #Print the number of hits for that genome to the screen.
          echo -e "Total number of unique hits in ${i}:" ;
          #Use wc -l to count the lines in the tab separated table created using awk to determine the total number of indvidual effectors found.
          wc -l < ./BigBlast_${BLAST}_${date}/${subject_filename}/Occurrences_of_${query_filename}_sequences_in_${subject_filename}.txt
          echo "Bed file generated: ./BigBlast_${BLAST}_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.bed"
          cp ${i} ./BigBlast_${BLAST}_${date}/${subject_filename}/${i} ;
          mv ${i}.* ./BigBlast_${BLAST}_${date}/${subject_filename}/BlastDB ;
        else
          #If no hits are found, print that no hits are found to the screen and end.
          echo "
          No hits found.
          "
          echo "##No hits found." >> ./BigBlast_${BLAST}_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.${BLAST}.outfmt6.evalue_used_1e-6.out
        fi
        done
    else
      for i in $(cat ${genomes}) ;
        do
        python -c "print('=' * 80)" ;
        #Create output file names
        subject_filename="${i##*/}" ;
        subject_extension="${subject_filename##*.}" ;
        subject_filename="${subject_filename%.*}" ;
        #Create output directories.
        mkdir ./BigBlast_${BLAST}_${date}/${subject_filename} ;
        mkdir ./BigBlast_${BLAST}_${date}/${subject_filename}/BlastDB ;
        #Print that ${BLAST} is being run and provide run information.
        echo "Running ${BLAST} for ${i}..." ;
        echo "out: ./BigBlast_${BLAST}_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.${BLAST}.outfmt6.evalue_used_1e-6.out"
        #Generate generic BLAST output.
        ${BLAST} -query ${query} -db ${i} -evalue 1e-6  -out ./BigBlast_${BLAST}_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.${BLAST}.evalue_1e-6.summarised.out ;
        #Perform ${BLAST}; output in outfmt 6 including all usual colummns plus the query coverage column.
        ${BLAST} -query ${query} -db ${i}  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" -evalue 1e-6  -out ./BigBlast_${BLAST}_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.${BLAST}.outfmt6.evalue_used_1e-6.out ;
        #Check whether any hits were found. If found, proceed with the next step.
        if [[ -s ./BigBlast_${BLAST}_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.${BLAST}.outfmt6.evalue_used_1e-6.out ]] ;
        then
          #Extract the hits which meet the percentage identity and percentage coverage threshold and prepare a txt file that has both contigs (effector and contig in genome) and a bed file that contains all hits.
          awk -v p="$PercentageIdent" -v c="$PercentageCov" '{ if ($3 >= p && $13>= c) print $0}' ./BigBlast_${BLAST}_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.${BLAST}.outfmt6.evalue_used_1e-6.out | awk '{ if ($11 < 1e-6) print $1,"\t",$2}' > ./BigBlast_${BLAST}_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.query_and_subject_hits.out ;
          #Create a tab separated table which contains the candidate effector hit and the number of hits for that effector.
          awk '{a[$1]++}END{for(i in a){print i,"\t",a[i]}}' ./BigBlast_${BLAST}_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.query_and_subject_hits.out > ./BigBlast_${BLAST}_${date}/${subject_filename}/Occurrences_of_${query_filename}_sequences_in_${subject_filename}.txt
          #Print the number of hits for that genome to the screen.
          echo -e "Total number of unique hits in ${i}:" ;
          #Use wc -l to count the lines in the tab separated table created using awk to determine the total number of indvidual effectors found.
          wc -l < ./BigBlast_${BLAST}_${date}/${subject_filename}/Occurrences_of_${query_filename}_sequences_in_${subject_filename}.txt
          cp ${i} ./BigBlast_${BLAST}_${date}/${subject_filename}/${i} ;
          mv ${i}.* ./BigBlast_${BLAST}_${date}/${subject_filename}/BlastDB ;
        else
          #If no hits are found, print that no hits are found to the screen and end.
          echo "
          No hits found.
          "
          echo "##No hits found." >> ./BigBlast_${BLAST}_${date}/${subject_filename}/${query_filename}_vs_${subject_filename}.${BLAST}.outfmt6.evalue_used_1e-6.out
        fi
        done
      fi

  mkdir ./BigBlast_${BLAST}_${date}/$query_filename
  cp $query ./BigBlast_${BLAST}_${date}/$query_filename
  mv $query.* ./BigBlast_${BLAST}_${date}/$query_filename

}

  function help(){
    echo "
    ###################################################################################################
    BigBlast.sh

    Usage: [PATH]/BigBlast.sh -q <FASTA> -g <TEXT FILE> -i <VALUE> -c <VALUE> [OPTIONAL ARGUMENTS]

    -q  FASTA file as input which contains a set of query sequences (effectors you're searching for).
    -g  TEXT file listing the subject FASTA files to search.
    -i  Percentage indentity threshold.
    -c  Percentage coverage threshold.
    -t  Type of BLAST search to be performed.
        You can choose from the following BLAST search types:
          * blastn  - nucleotide to nucleotide search
          * blastp  - protein to protein search
          * tblastn - translated protein to nuclotide
          * tblastx - translate nuclotide to protein to nucleotide. This is used for a nucleotide to
                      nucleotide search and will take into account frame shift.


    OPTIONAL ARGUMENT:

    -b  Generate a bedfile from each subject FASTA containing hits within the percenatge identity and
        coverage threshold with an evalue of < 1e-6.

    Note the subject FASTA files should be in the directory you wish to run the BigBlast.sh or
    the text files should contain the full path to the subject FASTA files.

    DEPENDENCIES
    BLAST: version 2.9.0+
    BEDTOOLS: version 2.25.0
    ##################################################################################################
    "
  }

  function help2(){
    echo "
      You can choose from the following BLAST search types:
        * blastn  - nucleotide to nucleotide search
        * blastp  - protein to protein search
        * tblastn - translated protein to nuclotide
        * tblastx  - translate nuclotide to protein to nucleotide. This is used for a nucleotide to
                    nucleotide search and will take into account frame shift.
    "
  }

  function check(){
  #  local OPTIND opt i
    while getopts ":uq:g:i:c:t:b" opt; do
      case $opt in
        u) usage=true ;;
        q) query="$OPTARG" ;;
        g) genomes="$OPTARG" ;;
        i) PercentageIdent="$OPTARG" ;;
        c) PercentageCov="$OPTARG" ;;
        t) BLAST="$OPTARG" ;;
        b) b=true ;;
        \?) help; exit 1
        ;;
      esac
    done
  #  shift $((OPTINF -1))

    if [ $usage ]
    then
      echo "
      Usage: [PATH]/BigBlast.sh -q <FASTA> -g <TEXT FILE> -i <VALUE> -c <VALUE> -t <BLASTSEARCHTYPE> [OPTIONAL ARGUMENTS]
      " ;
      exit
    fi

    if [ "$query" = "" ]
    then
      echo "
      Error: Please provide a FASTA file containing query sequences.
      " ;
      exit 2
    fi

    if [ "$genomes" = "" ]
    then
      echo "
      Error: Please provide a genome to be searched.
      " ;
      exit 3
    fi

    if [ "$PercentageIdent" = "" ]
    then
      echo "
      Error: Please provide a percentage identity threshold.
      " ;
      exit 3
    fi

    if [ "$PercentageCov" = "" ]
    then
      echo "
      Error: Please provide a percentage coverage threshold.
      " ;
      exit 3
    fi

    if [ "$BLAST" = "" ]
    then
      echo "
      Error: Please indicate the type of BLAST search." ;
      help2 ;
      exit 3
    fi
  }

  main $@


  mkdir BigBlast_${BLAST}_${date} #Make directories for the output.
