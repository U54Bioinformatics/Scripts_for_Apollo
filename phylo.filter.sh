#!/bin/bash
#SBATCH --ntasks=20
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH -p fast, abild
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ibishara@coh.org
#SBATCH -J phylo.filter
#SBATCH --workdir=./
#SBATCH --output=phylo_%j_filter.log

N=$SLURM_ARRAY_TASK_ID
if [ ! $N ]; then
    N=1
fi

echo "N: $N"

sample=`cat samples.txt | head -n $N | tail -n 1`


module load singularity

echo $(singularity --version)

#export PYTHONPATH=/home/ibishara/miniconda3/envs/betsy/bin/python2.7
export PATH="${PATH}:/home/ibishara/miniconda3/envs/betsy/lib64/python2.7/site-packages/genomicode/bin"

echo $(python --version)

betsy_run.py --network_plot $sample/parameter_screen/${sample}.pdf --num_cores 20 --receipt $sample/parameter_screen/${sample}.txt \
    --input PhyloSimpleVariantMatrix --input_file $sample/${sample}_clean.svm.txt \
    --dattr PhyloSimpleVariantMatrix.is_pseudobulk=yes \
    --input PhyloMetadataFile --input_file $sample/${sample}_full_phylo_meta.txt \
    --output PhyloFilterAnalysis --output_file $sample/parameter_screen/${sample}_30_15_10 \
    --dattr PhyloFilterAnalysis.is_pseudobulk=yes \
    --mattr keep_cells_with_perc_vars=30 \
    --mattr keep_vars_seen_in_perc_cells=15 \
    --mattr keep_vars_with_mixed_calls=10 \
    --mattr discard_close_vars=5 \
    --mattr discard_multiallelic_vars=yes \
    --run


betsy_run.py --network_plot $sample/parameter_screen/${sample}.pdf --num_cores 20 --receipt $sample/parameter_screen/${sample}.txt \
    --input PhyloSimpleVariantMatrix --input_file $sample/${sample}_clean.svm.txt \
    --dattr PhyloSimpleVariantMatrix.is_pseudobulk=yes \
    --input PhyloMetadataFile --input_file $sample/${sample}_full_phylo_meta.txt \
    --output PhyloFilterAnalysis --output_file $sample/parameter_screen/${sample}_20_15_10 \
    --dattr PhyloFilterAnalysis.is_pseudobulk=yes \
    --mattr keep_cells_with_perc_vars=20 \
    --mattr keep_vars_seen_in_perc_cells=15 \
    --mattr keep_vars_with_mixed_calls=10 \
    --mattr discard_close_vars=5 \
    --mattr discard_multiallelic_vars=yes \
    --run

betsy_run.py --network_plot $sample/parameter_screen/${sample}.pdf --num_cores 20 --receipt $sample/parameter_screen/${sample}.txt \
    --input PhyloSimpleVariantMatrix --input_file $sample/${sample}_clean.svm.txt \
    --dattr PhyloSimpleVariantMatrix.is_pseudobulk=yes \
    --input PhyloMetadataFile --input_file $sample/${sample}_full_phylo_meta.txt \
    --output PhyloFilterAnalysis --output_file $sample/parameter_screen/${sample}_15_15_10 \
    --dattr PhyloFilterAnalysis.is_pseudobulk=yes \
    --mattr keep_cells_with_perc_vars=15 \
    --mattr keep_vars_seen_in_perc_cells=15 \
    --mattr keep_vars_with_mixed_calls=10 \
    --mattr discard_close_vars=5 \
    --mattr discard_multiallelic_vars=yes \
    --run

betsy_run.py --network_plot $sample/parameter_screen/${sample}.pdf --num_cores 20 --receipt $sample/parameter_screen/${sample}.txt \
    --input PhyloSimpleVariantMatrix --input_file $sample/${sample}_clean.svm.txt \
    --dattr PhyloSimpleVariantMatrix.is_pseudobulk=yes \
    --input PhyloMetadataFile --input_file $sample/${sample}_full_phylo_meta.txt \
    --output PhyloFilterAnalysis --output_file $sample/parameter_screen/${sample}_30_10_10 \
    --dattr PhyloFilterAnalysis.is_pseudobulk=yes \
    --mattr keep_cells_with_perc_vars=30 \
    --mattr keep_vars_seen_in_perc_cells=10 \
    --mattr keep_vars_with_mixed_calls=10 \
    --mattr discard_close_vars=5 \
    --mattr discard_multiallelic_vars=yes \
    --run

betsy_run.py --network_plot $sample/parameter_screen/${sample}.pdf --num_cores 20 --receipt $sample/parameter_screen/${sample}.txt \
    --input PhyloSimpleVariantMatrix --input_file $sample/${sample}_clean.svm.txt \
    --dattr PhyloSimpleVariantMatrix.is_pseudobulk=yes \
    --input PhyloMetadataFile --input_file $sample/${sample}_full_phylo_meta.txt \
    --output PhyloFilterAnalysis --output_file $sample/parameter_screen/${sample}_20_10_10 \
    --dattr PhyloFilterAnalysis.is_pseudobulk=yes \
    --mattr keep_cells_with_perc_vars=20 \
    --mattr keep_vars_seen_in_perc_cells=10 \
    --mattr keep_vars_with_mixed_calls=10 \
    --mattr discard_close_vars=5 \
    --mattr discard_multiallelic_vars=yes \
    --run

betsy_run.py --network_plot $sample/parameter_screen/${sample}.pdf --num_cores 20 --receipt $sample/parameter_screen/${sample}.txt \
    --input PhyloSimpleVariantMatrix --input_file $sample/${sample}_clean.svm.txt \
    --dattr PhyloSimpleVariantMatrix.is_pseudobulk=yes \
    --input PhyloMetadataFile --input_file $sample/${sample}_full_phylo_meta.txt \
    --output PhyloFilterAnalysis --output_file $sample/parameter_screen/${sample}_20_5_5 \
    --dattr PhyloFilterAnalysis.is_pseudobulk=yes \
    --mattr keep_cells_with_perc_vars=20 \
    --mattr keep_vars_seen_in_perc_cells=5 \
    --mattr keep_vars_with_mixed_calls=5 \
    --mattr discard_close_vars=5 \
    --mattr discard_multiallelic_vars=yes \
    --run


betsy_run.py --network_plot $sample/parameter_screen/${sample}.pdf --num_cores 20 --receipt $sample/parameter_screen/${sample}.txt \
    --input PhyloSimpleVariantMatrix --input_file $sample/${sample}_clean.svm.txt \
    --dattr PhyloSimpleVariantMatrix.is_pseudobulk=yes \
    --input PhyloMetadataFile --input_file $sample/${sample}_full_phylo_meta.txt \
    --output PhyloFilterAnalysis --output_file $sample/parameter_screen/${sample}_10_5_5 \
    --dattr PhyloFilterAnalysis.is_pseudobulk=yes \
    --mattr keep_cells_with_perc_vars=10 \
    --mattr keep_vars_seen_in_perc_cells=5 \
    --mattr keep_vars_with_mixed_calls=5 \
    --mattr discard_close_vars=5 \
    --mattr discard_multiallelic_vars=yes \
    --run


betsy_run.py --network_plot $sample/parameter_screen/${sample}.pdf --num_cores 20 --receipt $sample/parameter_screen/${sample}.txt \
    --input PhyloSimpleVariantMatrix --input_file $sample/${sample}_clean.svm.txt \
    --dattr PhyloSimpleVariantMatrix.is_pseudobulk=yes \
    --input PhyloMetadataFile --input_file $sample/${sample}_full_phylo_meta.txt \
    --output PhyloFilterAnalysis --output_file $sample/parameter_screen/${sample}_15_5_5 \
    --dattr PhyloFilterAnalysis.is_pseudobulk=yes \
    --mattr keep_cells_with_perc_vars=15 \
    --mattr keep_vars_seen_in_perc_cells=5 \
    --mattr keep_vars_with_mixed_calls=5 \
    --mattr discard_close_vars=5 \
    --mattr discard_multiallelic_vars=yes \
    --run

betsy_run.py --network_plot $sample/parameter_screen/${sample}.pdf --num_cores 20 --receipt $sample/parameter_screen/${sample}.txt \
    --input PhyloSimpleVariantMatrix --input_file $sample/${sample}_clean.svm.txt \
    --dattr PhyloSimpleVariantMatrix.is_pseudobulk=yes \
    --input PhyloMetadataFile --input_file $sample/${sample}_full_phylo_meta.txt \
    --output PhyloFilterAnalysis --output_file $sample/parameter_screen/${sample}_15_10_10 \
    --dattr PhyloFilterAnalysis.is_pseudobulk=yes \
    --mattr keep_cells_with_perc_vars=15 \
    --mattr keep_vars_seen_in_perc_cells=10 \
    --mattr keep_vars_with_mixed_calls=10 \
    --mattr discard_close_vars=5 \
    --mattr discard_multiallelic_vars=yes \
    --run

betsy_run.py --network_plot $sample/parameter_screen/${sample}.pdf --num_cores 20 --receipt $sample/parameter_screen/${sample}.txt \
    --input PhyloSimpleVariantMatrix --input_file $sample/${sample}_clean.svm.txt \
    --dattr PhyloSimpleVariantMatrix.is_pseudobulk=yes \
    --input PhyloMetadataFile --input_file $sample/${sample}_full_phylo_meta.txt \
    --output PhyloFilterAnalysis --output_file $sample/parameter_screen/${sample}_30_5_5 \
    --dattr PhyloFilterAnalysis.is_pseudobulk=yes \
    --mattr keep_cells_with_perc_vars=30 \
    --mattr keep_vars_seen_in_perc_cells=5 \
    --mattr keep_vars_with_mixed_calls=5 \
    --mattr discard_close_vars=5 \
    --mattr discard_multiallelic_vars=yes \
    --run

