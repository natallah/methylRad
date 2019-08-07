#!/bin/sh -l
#PBS -N createSites
#PBS -q pccr
#PBS -l nodes=1:ppn=4
#PBS -l walltime=10:00:00
# Print the hostname of the compute node on which this job is running.
/bin/hostname

cd $PBS_O_WORKDIR

module load anaconda/5.1.0-py36

python3 4_CompareSites.py -i F9_UD_readsCatalogue &
python3 4_CompareSites.py -i F9_D4_readsCatalogue &
python3 4_CompareSites.py -i F9_D4_PG_readsCatalogue &
python3 4_CompareSites.py -i F9_D4_TCP_readsCatalogue &

wait