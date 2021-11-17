#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time=1-12:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ibishara@coh.org
#SBATCH -J phylo.filter
#SBATCH --workdir=./
#SBATCH --output=phylo_%j_filter.log
#SBATCH -p abild, fast



N=$SLURM_ARRAY_TASK_ID
if [ ! $N ]; then
    N=1
fi

echo "N: $N"

sample=`cat samples.txt | head -n $N | tail -n 1`




module load singularity

echo $(singularity --version)

export PYTHONPATH=/home/ibishara/miniconda3/envs/betsy/bin/python2.7
export PATH="${PATH}:/home/ibishara/miniconda3/envs/betsy/lib64/python2.7/site-packages/genomicode/bin"

echo $(python --version)

iter=100000000
burnin=95000000

betsy_run.py --network_plot $sample/BEAST2Analysis.pdf --num_cores 1 \
	--input PhyloSimpleVariantMatrix --input_file $sample/${sample}_filtered_svm.txt \
  	--dattr PhyloSimpleVariantMatrix.is_pseudobulk=yes \
	--input PhyloMetadataFile --input_file $sample/${sample}_full_phylo_meta.txt \
	--output BEAST2Analysis --output_file $sample/BEAST2Analysis \
	--dattr BEAST2Analysis.is_pseudobulk=yes \
	--dattr BEAST2Analysis.variant_annotations=cancer_cosmic \
	--mattr keep_cells_with_perc_vars=20 \
	--mattr keep_vars_seen_in_perc_cells=10 \
	--mattr keep_vars_with_mixed_calls=10 \
	--mattr discard_close_vars=5 \
	--mattr discard_multiallelic_vars=yes \
	--mattr beast2_clock_model=rln \
	--mattr beast2_site_model=gtr \
	--mattr beast2_tree_prior=yule \
	--mattr beast2_iterations=$iter \
	--mattr beast2_burnin=$burnin \
	--mattr tree_plot_height=6 \
	--mattr tree_plot_width=8 \
	--mattr tree_layout=roundrect \
	--mattr color_tree_by_category=yes \
	--mattr annovar_buildver=hg19 \
	--mattr cancer_genes_file=/temporarynfsmount/labs/abild/aritro/Phylogeny/cancer_genes.txt \
	--mattr cosmic_variants_file=/temporarynfsmount/labs/abild/aritro/Phylogeny/cosmic.v79.grch37.mutation_data.txt.gz \
	--run
