#!/bin/sh -l
#PBS -N countsites
#PBS -q pccr
#PBS -l nodes=1:ppn=1
#PBS -l walltime=100:00:00
# Print the hostname of the compute node on which this job is running.
/bin/hostname

cd $PBS_O_WORKDIR

module load anaconda/5.1.0-py36

python3 4b_VizCounts.py > runVizCounts.out
