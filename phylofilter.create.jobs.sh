#!/bin/bash
echo "Run batch ssGSEA"
sbatch --array 1-31 phylo.tree.sh