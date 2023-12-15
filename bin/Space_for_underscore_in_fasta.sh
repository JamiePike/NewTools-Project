#!/bin/bash

#Jamie Pike 16/12/2020
#Script for replacing spaced and tabs in FASTA heading for underscores

function main(){
  check "$@"

  if [ $a ]
  then
    if [[ $added = *[[:space:]]* ]]
    then
        to_add=$(sed $'s/ /_/g; s/\t/_/g;' <(printf '%s\n' "$added"))
    else
        to_add=$added
    fi

    python -c "print('=' * 80)"
    echo //Generating backup file: $infile.bak.
    echo //Converting $infile.
    echo //Adding text to $infile

    sed -i.bak $"s/>.*/&_"${to_add}"/; s/ /_/g; s/$(printf '\t')/_/g;" $infile  #Take the input fasta, append the extra text, then replace spaces and tabs with an underscore.

    python -c "print('=' * 80)"

  else

    python -c "print('=' * 80)"
    echo //Generating backup file: $infile.bak.
    echo //Converting $infile.

    sed -i.bak $'s/ /_/g; s/\t/_/g;'  $infile        #Take the input fasta and replace spaces and tabs with an underscore.

    python -c "print('=' * 80)"
  fi

 echo "It is recomended that $infile.bak is saved in a new directory.
Should Space for Underscore be run again, $infile.bak will be replaced and the original file lost."
 python -c "print('=' * 80)"
}

function help(){
  python -c "print('=' * 80)"
  echo "
  //Error: Argument not recognised.

  ##############################################################################
  Space for Underscore
  --------------------

  Usage: cmd -i <FASTA> [OPTIONS]

  -i Provide a FASTA file as input.

  OPTIONAL:

  -a Text which you wish to add to the end of the FASTA header, e.g. -a ""APPENED TEXT""

  ##############################################################################
  "
  python -c "print('=' * 80)"
}

function help2(){
  echo "
  ##############################################################################
  Space for Underscore
  --------------------

  Usage: cmd -i <FASTA> [OPTIONS]

  -i Provide a FASTA file as input.

  OPTIONAL:

  -a Text which you wish to add to the end of the FASTA header, e.g. -a ""APPENED TEXT""

  ##############################################################################
  "
  python -c "print('=' * 80)"
}

function check(){
#  local OPTIND opt i
  while getopts ":i:a:" opt; do
    case $opt in
      i) infile="$OPTARG" ;;
      a) a=true; added="$OPTARG";;
      \?) help; exit 1
      ;;
    esac
  done
#  shift $((OPTINF -1))

  if [ "$infile" = "" ]
  then
    python -c "print('=' * 80)"
    echo "
    //Please provide a input FASTA file." ;
    help2 ;
    exit 2
  fi
}

main "$@"
