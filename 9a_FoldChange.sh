#!/bin/bash
#PBS -q pccr
#PBS -l walltime=336:00:00
#PBS -l nodes=1:ppn=10
#PBS -l naccesspolicy=shared
#PBS -M psudyant@purdue.edu

starts=$(date +"%s")
start=$(date +"%r, %m-%d-%Y")

module load anaconda/5.1.0_py36

cd $PBS_0_WORKDIR

python 9_FoldChange.py > foldChange.out

ends=$(date +"%s")
end=$(date +"%r, %m-%d-%Y")
diff=$(($ends-$starts))
hours=$(($diff / 3600))
dif=$(($diff % 3600))
minutes=$(($dif / 60))
seconds=$(($dif % 60))
printf "\n\t===========Time Stamp===========\n"
printf "\tStart\t:$start\n\tEnd\t:$end\n\tTime\t:%02d:%02d:%02d\n" "$hours" "$minutes" "$seconds"
printf "\t================================\n\n"

qstat -f1 $PBS_JOBID | egrep 'Job Id|Job_Name|resources_used|exec_host'