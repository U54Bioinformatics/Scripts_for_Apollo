#!/bin/bash
echo "Run batch InferCNV"
sbatch --array 1-61 InferCNV.batch.sh
