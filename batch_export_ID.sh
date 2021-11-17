#!/bin/bash
echo "Run batch export ID"
sbatch --array 1-62 Export_ID.sh
