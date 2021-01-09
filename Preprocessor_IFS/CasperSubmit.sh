#!/bin/bash -l
#SBATCH -J RF
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH -t 24:00:00
#SBATCH -A P66770001
#SBATCH -p dav


# to start this scrip run:     for yy in {1..12}; do sbatch CasperSubmit.sh $yy; done

# to chack the status run "squeue -u $USER"

VAR="$1"

PythonName=Preprocessor_IFS.py

module load python/3.7.5
ncar_pylib
ml ncl nco
echo $VAR
srun ./$PythonName $VAR

exit
