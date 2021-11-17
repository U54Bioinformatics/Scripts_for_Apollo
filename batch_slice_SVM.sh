#!/bin/bash
echo "Run batch slice SVM"
sbatch --array 1-62 Slice_SVM.sh