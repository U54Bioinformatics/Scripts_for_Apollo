# Scripts_for_Apollo

## 1. Running inferCNV in parallel on multiple samples using SBATCH array 


```bash
anath@ppxhpcacc01$ bash InferCNV.create.jobs.sh
```

Requires these files: 
  1. Raw (filtered) counts file: *SAMPLE*.ic.counts.txt 
  2. InferCNV annotations file (label anticipated cancer cells as "malignant": *SAMPLE*.ic.annot.txt 
  3. sample.ic.txt: contains a list of *SAMPLE*
  4. hg19.RefSeq.NM_pos_unique_sort.txt: gene order file (static for each run)

  
`InferCNV.create.jobs.sh` creates an array of **N** jobs and invokes `InferCNV.batch.sh`. Modify **N** based on number of samples you will be submitting. 


`InferCNV.batch.sh` is the actual job submission script. Modify the SBATCH parameters as required. 


`InferCNV.R` is the R-script calling InferCNV code. Modify the path to the Conda enviroment and parmeters for InferCNV::run() based on requirements. The script will rename all non-malignant cells as "normal" in the annotations file. 


`InferCNV_Custom_Heatmap.R` can be used to display additional gene-level annotations from raw counts, Seurat annotations and perform clustering and cut tree for displaying on complex heatmaps.    

  
`InferCNV_make_counts_from_integrated_Seurat.R` can be used to export counts and annotations from an integrated Seurat object. The integrated Seurat object should contain counts after QC (i.e. after removal of low quality cells, doublets etc.) 
  
  

## 2. Installing BETSY on Apollo 
  
  1. Clone repository from Megatron to Apollo via Tikvah.
  2. Create a conda environment with python 2.7. 
  3. Activate conda environment and install required packages listed in README.
  4. Run `python setup.py build` and `python setup.py install`.
  5. Copy required Singularity containers from Megatron to Apollo via Tikvah.
  6. Follow instructions in README to set paths to the executables and Singularity containers in `~/.genomicoderc`  and `~/.betsyrc`.
  

  
  
  
