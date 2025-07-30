# Integrate:
#   
#   Skin ONLY!!!!!!!!
#
#   - NYU data
#   - external: 
#       - GSE173706 - HC
#   - external.3
#       - GSE195452 - Normal skins

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(gplots)
library(BiocManager)
library(EnhancedVolcano)
library(Matrix)
library(scales)
library(umap)
library(devtools)
library(iCellR)
library(tidyverse)
library(Matrix)
library(harmony)


########################################################################
# Load data
########################################################################
## load and process GSE195452 first because it is sparse !!!
path_dir <- '~/chang_HS/'#MZ changed path
load(paste0(path_dir, "GSE195452.normal.merged.Robj"))

# QC
GSE195452.normal.merged[["percent.mt"]] <- PercentageFeatureSet(GSE195452.normal.merged, pattern = "^MT-")
# filtering
GSE195452.normal.merged <- subset(GSE195452.normal.merged, subset = nFeature_RNA > 500 & percent.mt < 20)

# After QC, some layers will become empty
# Use the following code to remove empty layers
rna_assay <- GSE195452.normal.merged@assays[["RNA"]]

empty_layers <- sapply(rna_assay@layers, function(layer) {
  all(layer == 0) || is.null(layer) || sum(dim(layer)) == 0
})

empty_layer_names <- names(which(empty_layers))

if (length(empty_layer_names) > 0) {
  rna_assay@layers <- rna_assay@layers[!names(rna_assay@layers) %in% empty_layer_names]
  rna_assay@default <- rna_assay@default - length(empty_layer_names)
}

GSE195452.normal.merged@assays[["RNA"]] <- rna_assay

## load other data
load(paste0(path_dir, "nyu.HTO_all.noNegative.merged.Nov22_2024.Robj"))#MZ changed file name
load(paste0(path_dir, "GSE173706.NS.merged.Robj"))

skin.merged.obj <- merge(x = NYU.HTO_all.noNegative.merged, 
                    y = list(
                      GSE173706.NS.merged,
                      GSE195452.normal.merged
                    ))


# update condition info 
a <- skin.merged.obj@meta.data
a["Condition.New"][a['Condition.New'] == 'Normal'] <- 'normal'
a["Condition"][a['Condition'] == 'Normal'] <- 'normal'
skin.merged.obj@meta.data <- a

# view data statistics
table(skin.merged.obj$Patient, skin.merged.obj$Condition.New)
table(skin.merged.obj$Patient, skin.merged.obj$Condition)
table(skin.merged.obj$sampleType)
View(skin.merged.obj@meta.data)
table(skin.merged.obj$Condition.New, skin.merged.obj$Cohort)

# drop blood cells
skin.merged.obj <- subset(x = skin.merged.obj, subset = Condition.New != "HS-blood")
skin.merged.obj <- subset(x = skin.merged.obj, subset = Condition.New != "normal−blood")

table(skin.merged.obj$Condition.New, skin.merged.obj$Cohort)

Idents(skin.merged.obj) <- "Condition.New"
skin.merged.obj <- RenameIdents(skin.merged.obj,
                                'normal' = 'Normal',
                                'HS-lesion' = 'Lesion',
                                'HS-perilesion' = 'Perilesion',
                                'HS-lesion-maybe-w-tunnel' = 'Lesion',
                                'HS-lesion-nodule' = 'Lesion',
                                'HS-perilesion-tunnel' = 'Lesion',
                                'HS-lesion-tunnel' = 'Lesion')

skin.merged.obj$Condition.New.V2 <- Idents(skin.merged.obj)
table(skin.merged.obj$Condition.New.V2, skin.merged.obj$Cohort)

########################################################################
# Preprocessing
########################################################################

# QC
skin.merged.obj[["percent.mt"]] <- PercentageFeatureSet(skin.merged.obj, pattern = "^MT-")
#VlnPlot(skin.merged.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# filtering
skin.merged.obj <- subset(skin.merged.obj, subset = nFeature_RNA > 500 & percent.mt < 20)
skin.merged.obj=JoinLayers(skin.merged.obj)#MZ added
skin.merged.obj=split(skin.merged.obj,f=skin.merged.obj$Patient)#split by Patient

# processing
skin.merged.obj <- NormalizeData(skin.merged.obj, normalization.method = "LogNormalize", scale.factor = 10000)
skin.merged.obj <- FindVariableFeatures(skin.merged.obj, selection.method = "vst", nfeatures = 2000)
skin.merged.obj <- ScaleData(object = skin.merged.obj)

# dimension reduction
skin.merged.obj <- RunPCA(object = skin.merged.obj, features = VariableFeatures(object = skin.merged.obj))
skin.merged.obj <- RunUMAP(object = skin.merged.obj, dims = 1:30, min.dist = 0.4, n.neighbors = 50)


###################################################################
# Integration with Harmony for batch effect removal
# https://satijalab.org/seurat/articles/seurat5_integration
###################################################################
skin.merged.obj <- IntegrateLayers(object = skin.merged.obj, method = HarmonyIntegration, orig.reduction = "pca",
                              new.reduction = "harmony", verbose = FALSE)

skin.merged.obj <- RunUMAP(skin.merged.obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")

skin.merged.obj <- RunTSNE(object = skin.merged.obj, reduction = "harmony", dims = 1:30, perplexity=50, reduction.name = "tsne.harmony")


###################################################################
# Clustering
###################################################################
skin.merged.obj <- FindNeighbors(skin.merged.obj, reduction = "harmony", dims = 1:30)

skin.merged.obj <- FindClusters(skin.merged.obj, resolution = .6, cluster.name = "harmony_clusters")

skin.merged.obj <- FindClusters(skin.merged.obj, resolution = .8, cluster.name = "harmony_clusters.08")

###################################################### 
###### Identify top markers for each cluster
###################################################### 
skin.merged.obj <- JoinLayers(skin.merged.obj)
clusters <- unique(skin.merged.obj$harmony_clusters.08)

for (c in clusters){
  print(c)
  
  markers <- FindMarkers(skin.merged.obj, ident.1 = c, group.by = "harmony_clusters.08",
                         verbose = T, only.pos = T, logfc.threshold = 0.25)
  df <- data.frame(markers)
  write.csv(df, file=paste0(path_dir, "20250106_3datasets/_Top_markers_cluster_", c, ".csv"))
}

###################################################################
# Cluster Annotation
###################################################################
Idents(skin.merged.obj) <- "harmony_clusters.08"

skin.merged.obj <- RenameIdents(skin.merged.obj, 
                           '2' = 'Endothelial cells', '9' = 'Endothelial cells')

skin.merged.obj <- RenameIdents(skin.merged.obj, 
                           '11' = 'T cells', '14' = 'T cells')

skin.merged.obj <- RenameIdents(skin.merged.obj, 
                           '0' = 'Keratinocytes', '3' = 'Keratinocytes', '6' = 'Keratinocytes', 
                           '16' = 'Keratinocytes', '17' = 'Keratinocytes')

skin.merged.obj <- RenameIdents(skin.merged.obj, 
                           '5' = 'Plasma cells', '18' = 'Plasma cells', '19' = 'Plasma cells', '29' = 'Plasma cells')

skin.merged.obj <- RenameIdents(skin.merged.obj, 
                           '4' = 'Fibroblasts', '8' = 'Fibroblasts',
                           '12' = 'Fibroblasts', '15' = 'Fibroblasts', '28' = 'Fibroblasts')

skin.merged.obj <- RenameIdents(skin.merged.obj, 
                           '1' = 'Myofibroblast', '10'='Myofibroblast')

skin.merged.obj <- RenameIdents(skin.merged.obj, 
                           '13' = 'B cells')

skin.merged.obj <- RenameIdents(skin.merged.obj,
                           '7' = 'Myeloid cells')

skin.merged.obj <- RenameIdents(skin.merged.obj,
                           '22' = "Sweat gland & Myoepithelial cells")

skin.merged.obj$Cell_Type <- Idents(skin.merged.obj)

pdf(paste0(path_dir, "20250106_3datasets/__cell_type_umap.pdf"), width = 10, height = 8)
DimPlot(skin.merged.obj, reduction = 'umap.harmony', label = T, group.by = 'Cell_Type', 
        cols = c('T cells'='#ff758f', 'B cells'='#ff83c1', 
                 'Myeloid cells' = '#d5a19c',
                 'Keratinocytes'='#0096c7','Proliferating cells'='#9cafb7', 
                 'Fibroblasts'='#F6BD60', 'Sweat gland & Myoepithelial cells'='#a663cc',
                 'Myofibroblast' = '#87d5c6',
                 'Endothelial cells'='#CCB1F1',
                 'Plasma cells'='#f7a072',
                 'Melanocytes'='#eddeae'),
        raster.dpi = c(1024, 1024))

dev.off()


table(skin.merged.obj$Cell_Type, skin.merged.obj$Condition)
table(skin.merged.obj$Cell_Type, skin.merged.obj$Condition.New.V2)
table(skin.merged.obj$Cell_Type, skin.merged.obj$Cohort)
