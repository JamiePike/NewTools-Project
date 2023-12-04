#!/bin/bash 

#Command used to generate the S16 SPAdes Assembly 

one=${1?Enter first raw reads file.}
two=${2?Enter second raw reads file.} 
Iso=${3?Please enter the isolate code.}


python -c "print('=' * 75)"
echo "Spades Assembly for "${Iso}
echo "----------------------"
echo $(date)
python -c "print('=' * 75)"

#Run the Spades  command.
spades.py -1 ${one} -2 ${two} --isolate --cov-cutoff auto -o ./ 

python -c "print('=' * 75)"
#Send an email notification of job completion.
SendEmail.py "${Iso} SPAdes Assembly Complete" "Spades assembly for ${Iso} is now complete"

python -c "print('=' * 75)"
