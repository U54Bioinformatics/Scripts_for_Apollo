#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time=00:10:00
#SBATCH --workdir=./
#SBATCH -p abild,fast,all
#SBATCH -J export_ID
#SBATCH --output=ssGSEA_%j.log

N=$SLURM_ARRAY_TASK_ID
if [ ! $N ]; then
    N=1
fi

echo "N: $N"

sample=`cat samples.txt | head -n $N | tail -n 1`

module R 

Rscript Export_ID.R $sample