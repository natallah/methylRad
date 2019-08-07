#!/bin/sh -l
#PBS -N countsites
#PBS -q pccr
#PBS -l nodes=1:ppn=1
#PBS -l walltime=20:00:00
# Print the hostname of the compute node on which this job is running.
/bin/hostname

cd $PBS_O_WORKDIR

module load anaconda/5.1.0-py36

# python3 3_conv2CountMat.py -i tabFiles/F9_UD_MAPQ10_coord_sorted.tsv -o F9_UD_readsCatalogue &
# python3 3_conv2CountMat.py -i tabFiles/F9_D4_MAPQ10_coord_sorted.tsv -o F9_D4_readsCatalogue &
python3 3_conv2CountMat.py -i tabFiles/F9_D4_PG_MAPQ10_coord_sorted.tsv -o F9_D4_PG_readsCatalogue &
# python3 3_conv2CountMat.py -i tabFiles/F9_D4_TCP_MAPQ10_coord_sorted.tsv -o F9_D4_TCP_readsCatalogue &

wait
