library(ComplexHeatmap)
library(circlize)
library(data.table)
library(fastcluster)
library(RColorBrewer)

counts <-fread("counts.txt", sep = "\t", header=T, check.names=F) # The raw counts file that were used in InferCNV
chrs <- fread("hg19.RefSeq.NM_pos_unique_sort.txt", header=F)
colnames(chrs) <- c("Gene", "Chr", "Start", "Stop")

CNA_infer <- read.table("infercnv.observations.txt", sep="", check.name=F, header=T) #Output from InferCNV
out_name <- "CNA_infer"
annot_df <- fread("metadata.UMAPcluster.SingleR_anno.txt", header=T)


#Annotation
annot_df  <- annot_df[match(colnames(CNA_infer), annot_df$Cell.ID),]

#Counts
gene_list <- c("PDPN", "THY1", "CD34", "CDH11", "PTPRC", "GNLY", "CD3E", "CD3D", "CD68", "CD14", "FCGR3A", "CD8A", "ESR1", "EPCAM", "MSLN", "KRT8", "KRT18", "THY1", "CD34", "CDH11")

counts_sub <- counts[counts$gene_id %in% gene_list, ]
genes <- counts_sub$gene_id
counts_sub <- counts_sub[,colnames(counts_sub) %in% colnames(CNA_infer),with=FALSE]
counts_sub <- setcolorder(counts_sub, colnames(CNA_infer))
counts_sub <- log10(counts_sub +1 )
counts_sub <- as.data.table(t(counts_sub))
colnames(counts_sub) <- genes

#Transpose 
CNA_infer_M <- as.matrix(CNA_infer)
CNA_infer_M <- t(CNA_infer_M)

###Col Annotation
##
#Chromosome top annot InferCNV
chr_sub <- chrs[chrs$Gene %in% colnames(CNA_infer_M), ]
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
palette <- getPalette(24)
chr_names <- unique(chr_sub$Chr)
names(palette) <- chr_names
hc1 = HeatmapAnnotation(Chr = chr_sub$Chr, show_annotation_name=F, show_legend = F,
                        col = list(Chr=palette))



###Row Annotations
hr1 = rowAnnotation(
  Sample = annot_df$Sample.fixed, 
  col = list(Sample=c("S0011050"="#1B9E77", "S0024135"="#D95F02", "S0025863"="#7570B3")),
  na_col = "grey", border = TRUE, show_annotation_name = F)

hr2 = rowAnnotation(
  Gene_Number = anno_barplot(annot_df$nFeature_RNA, gp = gpar(col = "#666666")),
  na_col = "grey", border = TRUE, show_annotation_name=F)


col_Fibro = colorRamp2(c(0, 3), c("white", "red"))
col_Imm = colorRamp2(c(0, 3), c("white", "blue"))
col_Epi = colorRamp2(c(0, 3), c("white", "Green"))

hr3 = rowAnnotation(
  PDPN = as.matrix(counts_sub$PDPN),
  THY1 = as.matrix(counts_sub$THY1),
  #CD34 = as.matrix(counts_sub$CD34),
  #CDH11 = as.matrix(counts_sub$CDH11),
  CD45 = as.matrix(counts_sub$PTPRC),
  CD3E = as.matrix(counts_sub$CD3E),
  #CD3D = as.matrix(counts_sub$CD3D),
  #CD8A = as.matrix(counts_sub$CD8A),
  CD68 = as.matrix(counts_sub$CD68),
  #CD14 = as.matrix(counts_sub$CD14),
  #CD16 = as.matrix(counts_sub$FCGR3A),
  EPCAM = as.matrix(counts_sub$EPCAM),
  ESR1 = as.matrix(counts_sub$ESR1),
  KRT8 = as.matrix(counts_sub$KRT8),
  KRT18 = as.matrix(counts_sub$KRT18),
  col = list(PDPN=col_Fibro, THY1=col_Fibro, 
             CD45=col_Imm, CD3E=col_Imm, CD68=col_Imm,
             EPCAM=col_Epi, ESR1=col_Epi, KRT8=col_Epi, KRT18=col_Epi),
  border = TRUE,  gap = unit(0, "mm"),
  simple_anno_size = unit(3, "mm"),
  show_legend =rep(FALSE, 9),
  annotation_name_gp =gpar(fontsize=10)
)

getPalette2 = colorRampPalette(brewer.pal(8, "Accent"))
annot_names <- unique(annot_df$Encode_main_type)
palette2 <- getPalette2(length(annot_names))
names(palette2) <- annot_names

hr4 = rowAnnotation(
  Encode = annot_df$Encode_main_type, simple_anno_size = unit(3, "mm"),
  na_col = "grey", border = TRUE, show_annotation_name = T, annotation_name_gp =gpar(fontsize=10),
  col= list(Encode=palette2))

if(out_name == "CNA_infer"){
  heat_cols = c("navy", "white", "darkred")
} 
if(out_name == "HMM_infer"){
  heat_cols =colorRamp2(c(.5, 1, 1.5, 2, 3), c("royalblue1", "white", "#fee0d2", "#ef3b2c", "#67000d"))
}

####No Clustering
ha <- Heatmap(CNA_infer_M, cluster_rows = F, cluster_columns = F, show_row_names = F, show_heatmap_legend = FALSE,
              show_column_names = F, heatmap_width= unit(20, "cm"), heatmap_height= unit(10, "cm"), col = heat_cols,
              top_annotation = hc1, border=T,
              column_split = factor(chr_sub$Chr, levels= unique(chr_sub$Chr)), cluster_column_slices = FALSE, column_gap = unit(0, "mm"),
              row_split = factor(annot_df$Sample.fixed, levels= unique(annot_df$Sample.fixed)), cluster_row_slices = FALSE, 
              column_title_gp = gpar(fontsize=8), row_title = NULL)

ht_list <- ha + hr1 +hr2 + hr3 +hr4

png(paste(out_name, "_NoClust.png",sep=""), width = 36, height = 12, units = "cm", res=300)
ht <- draw(ht_list, ht_gap = unit(c(2,2,1,1), "mm"))
dev.off()




#Cluster method
dend <- fastcluster::hclust(dist(CNA_infer_M), method="ward.D2")

####K3
ha <- Heatmap(CNA_infer_M, cluster_rows = dend, cluster_columns = F, show_row_names = F, show_heatmap_legend = FALSE,
              show_column_names = F, heatmap_width= unit(20, "cm"), heatmap_height= unit(10, "cm"), col = heat_cols,
              clustering_method_rows = "ward.D2", top_annotation = hc1,
              row_split = 3, row_gap = unit(0, "mm"), border=T,
              column_split = factor(chr_sub$Chr, levels= unique(chr_sub$Chr)), cluster_column_slices = FALSE, column_gap = unit(0, "mm"),
              column_title_gp = gpar(fontsize=8))

ht_list <- ha + hr1 +hr2 + hr3 +hr4
png(paste(out_name, "_k3.png",sep=""), width = 36, height = 12, units = "cm", res=300)
ht = draw(ht_list, ht_gap = unit(c(2,2,1,1), "mm"))
dev.off()


#Cluster method
dend <- fastcluster::hclust(dist(CNA_infer_M), method="ward.D2")

####K3
ha <- Heatmap(CNA_infer_M, cluster_rows = dend, cluster_columns = F, show_row_names = F, show_heatmap_legend = FALSE,
              show_column_names = F, heatmap_width= unit(20, "cm"), heatmap_height= unit(10, "cm"), col = heat_cols,
              clustering_method_rows = "ward.D2", top_annotation = hc1, 
              row_split = 3, row_gap = unit(0, "mm"), border=T,
              column_split = factor(chr_sub$Chr, levels= unique(chr_sub$Chr)), cluster_column_slices = FALSE, column_gap = unit(0, "mm"),
              column_title_gp = gpar(fontsize=8))

ht_list <- ha + hr1 +hr2 + hr3 +hr4
png(paste(out_name, "_k3.png",sep=""), width = 36, height = 12, units = "cm", res=300)
ht = draw(ht_list, ht_gap = unit(c(2,2,1,1), "mm"))
dev.off()

#Get clusters
out_list <- row_order(ht)
names(out_list) <- c("1", "2", "3")

for (i in 1:length(out_list)){
  if (i == 1) {
    print(i)
    print(names(out_list)[i][1])
    temp <- row.names(CNA_infer_M[out_list[[i]],])
    out <- cbind(temp, paste("Cluster", names(out_list)[i][1], sep=""))
    colnames(out) <- c("Cell.ID", "Cluster")
  } 
  else {
    print(i)
    print(names(out_list)[i][1])
    clu <- row.names(CNA_infer_M[out_list[[i]],])
    clu <- cbind(clu, paste("Cluster", names(out_list)[i][1], sep=""))
    colnames(clu) <- c("Cell.ID", "Cluster")
    out <- rbind(out, clu)
  }
}
#Annot k3
colnames(out) <- c("Cell.ID", paste(out_name, "k3",sep=""))
out <- as.data.table(out)
annot_df <- merge(annot_df, out, by = "Cell.ID")

####K4
ha <- Heatmap(CNA_infer_M, cluster_rows = dend, cluster_columns = F, show_row_names = F, show_heatmap_legend = FALSE,
              show_column_names = F, heatmap_width= unit(20, "cm"), heatmap_height= unit(10, "cm"), col = heat_cols,
              clustering_method_rows = "ward.D2", top_annotation = hc1, 
              row_split = 4, row_gap = unit(0, "mm"), border=T,
              column_split = factor(chr_sub$Chr, levels= unique(chr_sub$Chr)), cluster_column_slices = FALSE, column_gap = unit(0, "mm"),
              column_title_gp = gpar(fontsize=8))

ht_list <- ha + hr1 +hr2 + hr3 +hr4
png(paste(out_name, "_k4.png",sep=""), width = 36, height = 12, units = "cm", res=300)
ht = draw(ht_list, ht_gap = unit(c(2,2,1,1), "mm"))
dev.off()

out_list <- row_order(ht)
names(out_list) <- c("1", "2", "3", "4")

#Get clusters
for (i in 1:length(out_list)){
  if (i == 1) {
    print(i)
    print(names(out_list)[i][1])
    temp <- row.names(CNA_infer_M[out_list[[i]],])
    out <- cbind(temp, paste("cluster", names(out_list)[i][1], sep=""))
    colnames(out) <- c("Cell.ID", "Cluster")
  } 
  else {
    print(i)
    print(names(out_list)[i][1])
    clu <- row.names(CNA_infer_M[out_list[[i]],])
    clu <- cbind(clu, paste("cluster", names(out_list)[i][1], sep=""))
    colnames(clu) <- c("Cell.ID", "Cluster")
    out <- rbind(out, clu)
  }
}
#Annot k4
colnames(out) <- c("Cell.ID", paste(out_name, "k4",sep=""))
out <- as.data.table(out)
annot_df <- merge(annot_df, out, by = "Cell.ID")

####K5
ha <- Heatmap(CNA_infer_M, cluster_rows = dend, cluster_columns = F, show_row_names = F, show_heatmap_legend = FALSE,
              show_column_names = F, heatmap_width= unit(20, "cm"), heatmap_height= unit(10, "cm"), col = heat_cols,
              clustering_method_rows = "ward.D2", top_annotation = hc1, 
              row_split = 5, row_gap = unit(0, "mm"), border=T,
              column_split = factor(chr_sub$Chr, levels= unique(chr_sub$Chr)), cluster_column_slices = FALSE, column_gap = unit(0, "mm"),
              column_title_gp = gpar(fontsize=8))

ht_list <- ha + hr1 +hr2 + hr3 +hr4
png(paste(out_name, "_k5.png",sep=""), width = 36, height = 12, units = "cm", res=300)
ht <- draw(ht_list, ht_gap = unit(c(2,2,1,1), "mm"))
dev.off()

#Get clusters
out_list <- row_order(ht)
names(out_list) <- c("1", "2", "3", "4", "5")
for (i in 1:length(out_list)){
  if (i == 1) {
    print(i)
    print(names(out_list)[i][1])
    temp <- row.names(CNA_infer_M[out_list[[i]],])
    out <- cbind(temp, paste("cluster", names(out_list)[i][1], sep=""))
    colnames(out) <- c("Cell.ID", "Cluster")
  } 
  else {
    print(i)
    print(names(out_list)[i][1])
    clu <- row.names(CNA_infer_M[out_list[[i]],])
    clu <- cbind(clu, paste("cluster", names(out_list)[i][1], sep=""))
    colnames(clu) <- c("Cell.ID", "Cluster")
    out <- rbind(out, clu)
  }
}
#Annot k5
colnames(out) <- c("Cell.ID", paste(out_name, "k5",sep=""))
out <- as.data.table(out)
annot_df <- merge(annot_df, out, by = "Cell.ID")

fwrite(annot_df, "Seurat.annot.singleR.clust.txt", sep="\t", quote = F)
