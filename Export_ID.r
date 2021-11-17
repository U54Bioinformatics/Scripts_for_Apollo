library(data.table)

args = commandArgs(trailingOnly = TRUE)
Sample <- as.character(args[1])
wd <- '/net/nfs-irwrsrchnas01/labs/abild/ibishara/Phylogeny/'

FEL_anno1 <- read.delim(file=paste(wd, "Feline_metadata_101421.txt", sep=""), header=T)

FEL_anno1$Outgroup <- ifelse(FEL_anno1$CancerNormal == "Cancer", "no", "yes")
FEL_anno1$orig.ident.short <- substr(FEL_anno1$orig.ident, 1, 6)
table(FEL_anno1$Outgroup, FEL_anno1$orig.ident.short)
cIN <- which(FEL_anno1$orig.ident.short == Sample)
Category <- if (as.integer(substr(Sample, 4, 6) ) >= 48) {substr(FEL_anno1$Cell[cIN], 10, 10)} else {substr(FEL_anno1$Cell[cIN], 8, 8)}
Cell <- if (as.integer(substr(Sample,4,6)) >= 48) {FEL_anno1$RandCell[cIN]} else {FEL_anno1$Cell[cIN]}
meta <- data.frame("Cell" = Cell,
                   "Category" = Category,
                   "Outgroup" = FEL_anno1$Outgroup[cIN])

write.table(meta, file = paste(wd, Sample, "/", Sample, "_full_phylo_meta.txt", sep=""), quote=F, row.names = F, col.names = T, sep="\t")
data.table::fwrite(list(meta$Cell), paste(wd, Sample, "/", Sample, '_cells.txt', sep=""))
