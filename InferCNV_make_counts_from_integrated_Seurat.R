
#### Make counts and annotations for InferCNV from integrated Seurat object ####

#f.seu is an integrated Seurat object with samples labelled with NEW_SAMPLE_ID

samples <- unique(f.seu$NEW_SAMPLE_ID)
counts <- GetAssayData(f.seu, assay = "RNA", slot = "counts")

dir.create("InferCNV")

set.seed(1234)
cIN.ref <- sample(which(f.seu$encode != "Epithelial cells"), size = 3000, replace = F) #f.seu$encode contains annotations from SinglerR encode 

for (i in 1:length(samples)) {
  cIN.epi <- which(f.seu$NEW_SAMPLE_ID == samples[i] & f.seu$encode == "Epithelial cells")
  split.counts <- data.frame(counts[, c(cIN.ref, cIN.epi)])
  fwrite(split.counts, file=paste("./InferCNV/", samples[i], ".ic.counts.txt", sep=""), row.names = T, col.names = T, quote=F, sep="\t")  
      
  A <- colnames(split.counts)
  B <- gsub("Epithelial cells", "malignant", f.seu$encode[c(cIN.ref, cIN.epi)])
  split.annot <- data.frame(A, B)
  fwrite(split.annot, file=paste("./InferCNV/", samples[i], ".ic.annot.txt", sep=""), row.names = F, col.names = F, quote=F, sep="\t")    
}

fwrite(data.frame(samples), file="samples.ic.txt", quote=F, col.names=F, sep="\n")
