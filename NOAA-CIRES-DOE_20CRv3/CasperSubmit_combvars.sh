#!/bin/bash -l
#SBATCH -J RF
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH -t 24:00:00
#SBATCH -A P66770001
#SBATCH -p dav


# to start this scrip run  sbatch CasperSubmit_daystats.sh program"
# variable can be [PRES PSL PW Q2 Q500 T2 T2min T2max T500 U10 U200 U500 V10 V200 V500 V850]

# to chack the status run "squeue -u $USER"

PythonName="$1"

module load python/3.7.5
ncar_pylib
ml ncl nco
echo $PythonName
srun ./$PythonName

exit
