library(Seurat)
library(dplyr)
library(ggplot2)
source("~/kwanho/src/seurat_tools.R")

seur1 = readRDS("/stanley/levin_dr/kwanho/projects/vahbiz/redo_all/data/obj_3_no_doublet_no_astro_neurons.rds")
seur2 = readRDS("/stanley/levin_dr/kwanho/projects/vahbiz/analysis_20211210/new_data_analysis/seurat/post_QC/seur_Oligo_rep2_CB_noDoublet_noAstro_v2.rds")

seur1$replicate = "rep1"
seur2$replicate = "rep2"

seur1$sample = paste0(seur1$age, '_', seur2$layer)

mat1 = seur1@assays$RNA@counts
mat2 = seur2@assays$RNA@counts
common.genes = intersect(rownames(mat1), rownames(mat2))
mat = cbind(mat1[common.genes,], mat2[common.genes,])

meta1 = seur1@meta.data
meta2 = seur2@meta.data
my.cols = c("nCount_RNA","nFeature_RNA","age","layer","sample","replicate","percent.mt","percent.rb")
meta = rbind(meta1[,my.cols], meta2[,my.cols])

seur <- CreateSeuratObject(counts=mat, project='Oligo', meta.data=meta)

seur <- seur %>%
        SCTransform(vars.to.regress=c('percent.mt', 'age', 'replicate')) %>%
        RunPCA() %>%
        FindNeighbors(dims=1:30) %>%
        RunUMAP(dims=1:30) %>%
        FindClusters(resolution=0.3)

saveRDS(seur, "seur_Oligo_merged_processed_v2.rds")

pdf("umap_Oligo_merged_v2.pdf", height=9, width=10)
DimPlot(seur, reduction='umap', pt.size=.5) + NoAxes()
DimPlot(seur, reduction='umap', pt.size=.5, group.by='replicate') + NoAxes()
DimPlot(seur, reduction='umap', pt.size=.5, group.by='sample') + NoAxes()
DimPlot(seur, reduction='umap', pt.size=.5, group.by='age') + NoAxes()
DimPlot(seur, reduction='umap', pt.size=.5, group.by='layer') + NoAxes()
dev.off()

