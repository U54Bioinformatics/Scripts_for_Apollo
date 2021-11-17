#!/bin/bash
echo "Run batch phylo filter"
sbatch --array 1-62 phylo.filter.sh