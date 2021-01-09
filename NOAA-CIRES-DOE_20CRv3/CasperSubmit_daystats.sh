#!/bin/bash -l
#SBATCH -J RF
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH -t 24:00:00
#SBATCH -A P66770001
#SBATCH -p dav


# to start this scrip run  sbatch CasperSubmit_daystats.sh variable"
# variable can be [PRES PSL PW Q2 Q500 T2 T2min T2max T500 U10 U200 U500 V10 V200 V500 V850 ZG500]

# to chack the status run "squeue -u $USER"

VAR="$1"

PythonName=3-hourly_to_daily_files.py

module load python/3.7.5
ncar_pylib
ml ncl nco
echo $VAR
srun ./$PythonName $VAR

exit
