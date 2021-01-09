#!/bin/bash -l
#SBATCH -J RF
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50G
#SBATCH -t 24:00:00
#SBATCH -A P66770001
#SBATCH -p dav
#SBATCH --output=R-%x.%j.out
#SBATCH --error=R-%x.%j.err


# to start this scrip run -- sbatch CasperSubmit_regions.sh REG"
# to run all regions run -- for ii in {0..7}; do sbatch CasperSubmit_regions.sh $ii; done

# to chack the status run "squeue -u $USER"

REG="$1"

PythonName=NOAA20C_XWTs_in_Model.py

module load python/3.7.5
ncar_pylib
ml ncl nco
echo $VAR
srun ./$PythonName $REG

exit
