#!/bin/bash 

#Quast Command used to analyse SPAdes Assembly

Assemb=${1?No Assembly file given.}

python -c "print('=' * 75)"
echo "Quast Analysis: " ${Assemb}
echo "-------------------------------------------"
echo $(date)
python -c "print('=' * 75)"

echo "Quast Analysis"
echo "--------------"
quast.py --fungus --circos -o QuastResult ${Assemb} 

python -c "print('=' * 75)"
echo "BUSCO Analysis"
echo "--------------"
busco -i ${Assemb} -o busco_Hypocreales -l hypocreales -m geno -c 1 -f 

python -c "print('=' * 75)"
SendEmail.py "Assembly Quality Analysis Complete"  "Quast and BUSCO assessment complete for current spades Assembly" 
python -c "print('=' * 75)"
echo $(date)
python -c "print('=' * 75)"
