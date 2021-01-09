#!/bin/bash -l
#SBATCH -J RF
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100G
#SBATCH -t 24:00:00
#SBATCH -A P66770001
#SBATCH -p dav


# to start this scrip run "for ii in {0..7}; do echo $ii; sbatch CasperSubmit.sh $ii; done"
# to chack the status run "squeue -u $USER"

HUCnr="$1"

echo $HUCnr

# module load python/3.7.5
ncar_pylib
ml ncl nco
ncar_pylib /glade/work/$USER/PYTHON_CASPER_clone
srun ./SearchOptimum_XWT-Restarts.py $HUCnr

exit
