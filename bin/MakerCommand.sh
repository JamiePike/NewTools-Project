 #!/bin/bash 

#Maker Command | Jamie Pike

#Run like this: nohup MakerCommand.sh 1>MakerRun_X.log

export LIBDIR=/home/u1983390/apps/RepeatMaskerSetup/RepeatMasker/Libraries/ 
export REPEATMASKER_LIB_DIR=/home/u1983390/apps/RepeatMaskerSetup/RepeatMasker/Libraries/
export REPEATMASKER_MATRICES_DIR=/home/u1983390/apps/RepeatMaskerSetup/RepeatMasker/Libraries/Matrices

set -e

python -c "print('=' * 75)"
echo "Maker Run" $(date)
python -c "print('=' * 75)"

/home/u1983390/apps/maker/bin/maker -RM_off maker_bopts.ctl maker_exe.ctl maker_opts.ctl 

python -c "print('=' * 75)"


SendEmail.py "Maker Annotation Complete" "Hi Jamie, the Maker annotation is now complete. :)"
