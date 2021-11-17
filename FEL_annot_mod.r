
library(data.table)

FEL_anno1 <- read.delim(file="~/Desktop/Phylogeny/Feline_metadata_101421.txt", header=T)
FEL_anno1$ARM
FEL_anno1 [which(FEL_anno1$Patient == 'FEL048'),]

Sample <- 'FEL048' 

colnames(FEL_anno1)

FEL_anno1$Outgroup <- ifelse(FEL_anno1$CancerNormal == "Cancer", "no", "yes")

FEL_anno1$orig.ident.short <- substr(FEL_anno1$orig.ident,1,6)

table(FEL_anno1$Outgroup, FEL_anno1$orig.ident.short)

cIN <- which(FEL_anno1$orig.ident.short == Sample)

Category <- if (as.integer(substr(Sample,4,6) ) >= 48) { substr(FEL_anno1$Cell[cIN], 10, 10) } else { substr(FEL_anno1$Cell[cIN], 8, 8) }
Cell <- if (as.integer(substr(Sample,4,6)) >= 48) {FEL_anno1$RandCell[cIN]} else {FEL_anno1$Cell[cIN]}

meta1 <- data.frame("Cell"=Cell, 
                   "Category"= Category, 
                   "Outgroup"=FEL_anno1$Outgroup[cIN])

write.table(meta, file=paste("~/Desktop/Phylogeny/",Sample,"/" ,Sample,"_full_phylo_meta.txt", sep=""), quote=F, row.names = F, col.names = T, sep="\t")


FEL <- fread(paste("~/Desktop/Phylogeny/",Sample,"/" ,Sample,"_full_phylo_meta.txt", sep= ""))
data.table::fwrite(list(FEL$Cell), paste("~/Desktop/Phylogeny/",Sample,"/" ,Sample,'_cells.txt', sep=""))
# set.seed(1234)
# cSUB_1 <- sort(sample(which(meta$Outgroup == "no"), size = 50, replace = F))
# cSUB_2 <- sort(sample(which(meta$Outgroup == "yes"), size = 50, replace = F))

# sub.meta <- meta[c(cSUB_1, cSUB_2), ]

# write.table(sub.meta, file="~/Desktop/Phylogeny/FEL013_sub_phylo_meta.txt", quote=F, row.names = F, col.names = T, sep="\t")