#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=InferCNV.sh.%A_%a.stdout
#SBATCH -p abild,fast,all
#SBATCH --workdir=./
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=anath@coh.org

start=`date +%s`

CPU=$SLURM_NTASKS
if [ ! $CPU ]; then
   CPU=2
fi

N=$SLURM_ARRAY_TASK_ID
if [ ! $N ]; then
    N=1
fi

echo "CPU: $CPU"
echo "N: $N"

sample=`cat samples.ic.txt | head -n $N | tail -n 1`

/home/anath/miniconda3/envs/INFERCNV/bin/Rscript InferCNV.R $sample 8

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"
