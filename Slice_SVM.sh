#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --workdir=./
#SBATCH -p abild,fast,all
#SBATCH -J slice_SVM
#SBATCH --output=slice_SVM_%j.log

N=$SLURM_ARRAY_TASK_ID
if [ ! $N ]; then
    N=1
fi

echo "N: $N"

sample=`cat samples.txt | head -n $N | tail -n 1`

export PYTHONPATH=/home/ibishara/miniconda3/envs/betsy/bin/python2.7

python slice_svm.py --keep_samples_from_file_lowmem $PWD/$sample/${sample}_cells.txt $PWD/$sample/phylogenetics.svm.txt.gz > $PWD/$sample/${sample}_clean.svm.txt
gzip $PWD/$sample/${sample}_clean.svm.txt