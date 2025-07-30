# This script conducts DEG analysis using DESeq2, adjusting for batch (sample) and sequencing platform as covariates

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(gplots)
library(BiocManager)
library(EnhancedVolcano)
library(Matrix)
library(scales)
library(DESeq2)
library(data.table)

path_dir='~/chang_HS/'
load(paste0(path_dir, "20250106_3datasets/allSkin.merged.clustering.celltype.0106.Robj"))

######################################################### 
####### Psudobulk Differential Expression Analysis
#########################################################
##################### lesion vs normal ##################### 
#### load inclusive genes - protein coding gene list
p_df = read.csv(paste0(path_dir, "protein-coding_gene.txt"), sep='\t')
protein_coding_genes <- p_df$symbol

#### load exclusive genes - genes only express within single dataset
eg_df <- read.csv(paste0(path_dir, "20250106_3datasets/single_cohort_genes.txt"), sep='\t')
exclusive_genes <- eg_df$Gene

skin.merged.lesion.normal <- subset(skin.merged.obj, subset= Condition.New.V2 %in% c("Normal", "Lesion"))

#### prepare psudobulk RNA data
DefaultAssay(skin.merged.lesion.normal)

# meta data
sample_meta <- skin.merged.lesion.normal@meta.data[, c("Patient", "platform", "Condition", "Condition.New.V2", "Cohort")]
sample_meta_unique <- sample_meta[!duplicated(sample_meta$Patient), ]
rownames(sample_meta_unique) <- sample_meta_unique$Patient

cell_list = c(
  'Myofibroblast',
  'T cells',
  'B cells',
  'Myeloid cells',
  'Fibroblasts',
  'Keratinocytes',
  'Endothelial cells',
  'Plasma cells',
  'Sweat gland & Myoepithelial cells'
)

main_dir <- paste0(path_dir, "20250106_3datasets")
sub_dir <- "DEGs_lesion_vs_normal.psudobulk"
output_dir <- file.path(main_dir, sub_dir)
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
} else {
  print("Dir already exists!")
}

for (c in cell_list){
  print(paste0(c, "-------------------------------------"))
  subset.obj = subset(skin.merged.lesion.normal, subset = Cell_Type == c)
  cnts <- AggregateExpression(subset.obj, 
                              group.by = c('Patient'),
                              assays = 'RNA',
                              #slot = 'counts',
                              return.seurat = FALSE)
  
  cnts <- as.matrix(cnts$RNA) + 1
  cnts <- cnts[!rownames(cnts) %in% exclusive_genes, ]
  cnts <- cnts[rownames(cnts) %in% protein_coding_genes, ]
  
  sample_meta_unique_sorted <- sample_meta_unique[match(colnames(cnts), sample_meta_unique$Patient), ]
  
  dds <- DESeqDataSetFromMatrix(countData = cnts,
                                colData = sample_meta_unique_sorted,
                                design = ~ platform + Condition.New.V2)
  
  #keep <- rowSums(counts(dds)) >= 10
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep,]
  
  dds$Condition.New.V2 <- factor(dds$Condition.New.V2, levels=c('Normal', 'Lesion'))
  
  dds <- DESeq(dds)
  
  res <- results(dds)
  
  res_df <- as.data.frame(res)
  
  out_csv <- paste(c, '_filtered.csv', sep = '')
  write.csv(res_df, file.path(main_dir, sub_dir, out_csv), row.names=TRUE)
  
  
  out_idr <- paste(c, '.pdf', sep = '')
  out_idr <- file.path(main_dir, sub_dir, out_idr)
  volcanoplot <- EnhancedVolcano(res_df, lab = rownames(res_df), x = 'log2FoldChange', y = 'padj',
                                 pCutoff = 0.05, FCcutoff = 0.5, 
                                 pointSize = 3.0,labSize = 6.0,
                                 legendLabSize = 16,
                                 legendIconSize = 5.0,
                                 title = c)
  ggsave(out_idr, plot = volcanoplot, width=10, height = 12)
}


# ######################################################### 
# ####### Differential Expression Analysis - Considering cells as samples
# ####### Using DESeq2 package
# #########################################################
#### load inclusive genes - protein coding gene list
p_df = read.csv(paste0(path_dir, "protein-coding_gene.txt"), sep='\t')
protein_coding_genes <- p_df$symbol
#### load exclusive genes - genes only express within single dataset
eg_df <- read.csv(paste0(path_dir, "20250106_3datasets/single_cohort_genes.txt"), sep='\t')
exclusive_genes <- eg_df$Gene

cell_list = c(
  'Myofibroblast',
  'T cells',
  'B cells',
  'Myeloid cells',
  'Fibroblasts',
  'Keratinocytes',
  'Endothelial cells',
  'Plasma cells',
  'Sweat gland & Myoepithelial cells'
)

main_dir <- paste0(path_dir, "20250106_3datasets")
sub_dir.2 <- "DEGs_lesion_vs_normal.DESeq2"
output_dir <- file.path(main_dir, sub_dir.2)
if (!dir.exists(output_dir)){
  dir.create(output_dir,recursive=T)
} else {
  print("Dir already exists!")
}

skin.merged.lesion.normal <- subset(skin.merged.obj, subset= Condition.New.V2 %in% c("Normal", "Lesion"))

for (c in cell_list){
  print(paste0(c, "-------------------------------------"))
  subset.obj = subset(skin.merged.lesion.normal, subset = Cell_Type == c)
  # meta data
  cell_meta <- subset.obj@meta.data[, c("cell", "Patient", "platform", "Condition", "Condition.New.V2", "Cohort", "Cell_Type")]
  #View(cell_meta)

  # get raw counts
  raw_counts <- as.matrix(subset.obj@assays$RNA$counts) + 1
  raw_counts <- raw_counts[!rownames(raw_counts) %in% exclusive_genes, ]
  raw_counts <- raw_counts[rownames(raw_counts) %in% protein_coding_genes, ]
  
  dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                                colData = cell_meta,
                                design = ~ platform + Condition.New.V2)

  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]

  dds$Condition.New.V2 <- factor(dds$Condition.New.V2, levels=c('Normal', 'Lesion'))#MZ: change Condition to Condition.New.V2. change to captial

  dds <- DESeq(dds, test='LRT', reduced = ~ platform,
               useT=TRUE, minmu=1e-6, minReplicatesForReplace=Inf, fitType = "glmGamPoi")

  res <- results(dds)

  res_df <- as.data.frame(res)

  out_csv <- paste(c, '_filtered.csv', sep = '')
  write.csv(res_df, file.path(main_dir, sub_dir.2, out_csv), row.names=TRUE)


  out_idr <- paste(c, '.pdf', sep = '')
  out_idr <- file.path(main_dir, sub_dir.2, out_idr)
  volcanoplot <- EnhancedVolcano(res_df, lab = rownames(res_df), x = 'log2FoldChange', y = 'padj',
                                 pCutoff = 0.05, FCcutoff = 0.5,
                                 pointSize = 3.0,labSize = 6.0,
                                 legendLabSize = 16,
                                 legendIconSize = 5.0,
                                 title = c)
  ggsave(out_idr, plot = volcanoplot, width=10, height = 12)
  gc()
}


