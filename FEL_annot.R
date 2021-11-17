
FEL_anno1 <- read.delim(file="~/Desktop/Phylogeny/FEL001046_scRNA.metadata.clinical.txt", header=T)

Sample <- 'FEL012' 

colnames(FEL_anno1)

FEL_anno1$Outgroup <- ifelse(FEL_anno1$Celltype == "Cancer cells", "no", "yes")

table(FEL_anno1$Outgroup, FEL_anno1$orig.ident)

cIN <- which(FEL_anno1$orig.ident == Sample)
meta <- data.frame("Cell"=FEL_anno1$Cell.ID[cIN], 
                   "Category"=FEL_anno1$Timepoint[cIN], 
                   "Outgroup"=FEL_anno1$Outgroup[cIN])

write.table(meta, file=paste("~/Desktop/Phylogeny/",Sample,"_full_phylo_meta.txt", sep=""), quote=F, row.names = F, col.names = T, sep="\t")



# set.seed(1234)
# cSUB_1 <- sort(sample(which(meta$Outgroup == "no"), size = 50, replace = F))
# cSUB_2 <- sort(sample(which(meta$Outgroup == "yes"), size = 50, replace = F))

# sub.meta <- meta[c(cSUB_1, cSUB_2), ]

# write.table(sub.meta, file="~/Desktop/Phylogeny/FEL013_sub_phylo_meta.txt", quote=F, row.names = F, col.names = T, sep="\t")
