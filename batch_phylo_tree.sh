#!/bin/bash
echo "Run batch phylo tree"
sbatch --array 1-62 phylo.tree.sh