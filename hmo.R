# Title: hmo.R
# Author: Renee Oles
# Date: Sep 21 2022
# Purpose: create heatmap from proteomics data for hmo project 
# Publication: 

# Load libraries
library(pheatmap)
library(paletteer)
library(tidyverse)
library(dplyr)

# Edit the format of the data frame based off hmo_proteomics.Rmd output 
rownames(dfANOVA_edit) <- rownames(dfANOVA)
dfANOVA_edit$ProteinID <- rownames(dfANOVA_edit)
dfANOVA_edit <- left_join(dfANOVA_edit,protein_id)
#colnames(proteins) <- "ProteinID"
proteins <- left_join(proteins,dfANOVA_edit)
proteins <- dfANOVA_edit

# Create a protein matrix by filtering for differentially abundant proteins
protein_matrix <- proteins[proteins$MaxLog2FC >=2 & proteins$pval_BH_adj < 0.05,c(1:3)]
rownames(protein_matrix) <- proteins[proteins$MaxLog2FC >=2 & proteins$pval_BH_adj < 0.05,]$Name
#protein_matrix <- protein_matrix[protein_matrix$TP1 > 9 & protein_matrix$TP2 > 9 & protein_matrix$TP3 > 9, ]
protein_matrix$intermediateLocusTag <- rownames(protein_matrix)
protein_matrix <- left_join(protein_matrix, tag)
other <- protein_matrix[!complete.cases(protein_matrix),]
protein_matrix <- na.omit(protein_matrix)
other <- other[,c(1:4)]
colnames(other)[4] <- "gene"
other <- left_join(other,tag)
protein_matrix <- protein_matrix[,c(1,2,3,6)]
other <- other %>% 
  mutate(originalLocusTag = coalesce(originalLocusTag,gene))
other <- other[,c(1:3,6)]
protein_matrix <- rbind(protein_matrix,other)
protein_matrix <- left_join(protein_matrix, tag)
select <- left_join(select, protein_matrix)

# Create a heatmap from the protein matrix
plot_matrix <- as.matrix(protein_matrix[,c(1:3)])
rownames(plot_matrix) <- protein_matrix$intermediateLocusTag
plot_matrix <- as.matrix(select[,c(2:4)])
rownames(plot_matrix) <- select$originalLocusTag
coul <- paletteer_c("ggthemes::Red-Blue-White Diverging", 30)
coul <- paletteer_c("ggthemes::Blue", 30)
tiff(file="proteomics_bf.tiff",res=300,height=12,width=4, units="in")
pheatmap(plot_matrix, col=coul, cluster_cols=F, 
         fontsize_col=10,fontsize_row=5, border_color=NA)
dev.off()
