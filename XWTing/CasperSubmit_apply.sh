#!/bin/bash -l
#SBATCH -J RF
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50G
#SBATCH -t 0:25:00
#SBATCH -A P66770001
#SBATCH -p dav
#SBATCH --constraint=x11


# to start this scrip run "for ii in {0..18}; do echo $ii; sbatch CasperSubmit_apply.sh $ii; done"
# to chack the status run "squeue -u $USER"

HUCnr="$1"

echo $HUCnr

module load python/2.7.14
ncar_pylib
ml ncl nco
srun ./Centroids-and-Scatterplot.py $HUCnr

exit
