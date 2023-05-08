library(Seurat)
library(ggplot2)
library(dplyr)
library(pracma)

seur <- readRDS("seur_initial_CB.rds")

# remove predicted doublets (134 cells)
meta = readRDS("../../scrublet/metadata_oligo_rep2_scrublet_res.rds")
rm.cells = rownames(meta)[which(meta$doublet_pred>0)]

seur <- subset(seur, cells = rm.cells, invert=T)

seur$age <- ordered(as.factor(seur$age), levels=c("P30", "P90"))

pdf("qc_CB_filtered.pdf")
VlnPlot(seur, features=c('nFeature_RNA', 'nCount_RNA', 'percent.mt', 'percent.rb'), group.by='age', ncol=2, pt.size=0) + NoLegend()
FeatureScatter(seur, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by='age')
FeatureScatter(seur, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by='age')
dev.off()

pdf("qc_CB_hists.pdf")
hist(seur$nFeature_RNA[seur$age=='P30'], col='skyblue3', breaks=linspace(0,8000, 100))
hist(seur$nFeature_RNA[seur$age=='P90'], col='darkorange', breaks=linspace(0,8000, 100))
hist(seur$nCount_RNA[seur$age=='P30'], col='skyblue3', breaks=linspace(0,max(seur$nCount_RNA), 100))
hist(seur$nCount_RNA[seur$age=='P90'], col='darkorange', breaks=linspace(0,max(seur$nCount_RNA), 100))
dev.off()


# increased max nFeature cutoff from 3500 to 6500
cat("QC filtering -- nFeature_RNA > 500 & nFeature_RNA < 6500 & nCount_RNA > 1500 & percent.mt < 10 \n")
sseur <- subset(seur, subset = nFeature_RNA > 500 & nFeature_RNA < 6500 & nCount_RNA > 1500 & percent.mt < 15)
pdf("qc_CB_filtered_v2.pdf")
VlnPlot(sseur, features=c('nFeature_RNA', 'nCount_RNA', 'percent.mt', 'percent.rb'), group.by='age', ncol=2, pt.size=0) + NoLegend()
FeatureScatter(sseur, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by='age')
FeatureScatter(sseur, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by='age')
dev.off()

sseur <- sseur %>%
	SCTransform(vars.to.regress=c('percent.mt', 'percent.rb', 'age')) %>% 
	RunPCA() %>% 
	FindNeighbors(dims=1:30) %>% 
	RunUMAP(dims=1:30) %>% 
	FindClusters(resolution=0.3)

saveRDS(sseur, "seur_Oligo_rep2_processed_v2.rds")

pdf("umap_Oligo_rep2_processed_v2.pdf", height=9, width=10)
DimPlot(sseur, reduction='umap', pt.size=.5, label=T) + NoAxes()
DimPlot(sseur, reduction='umap', pt.size=.5, group.by='age') + NoAxes()
DimPlot(sseur, reduction='umap', pt.size=.5, group.by='layer') + NoAxes()
dev.off()

# remove astro = cluster11
sseur <- subset(sseur, idents=11, invert=T)

sseur <- sseur %>%
	SCTransform(vars.to.regress=c('percent.mt', 'percent.rb', 'age')) %>% 
	RunPCA() %>% 
	FindNeighbors(dims=1:30) %>% 
	RunUMAP(dims=1:30) %>% 
	FindClusters(resolution=0.3)

saveRDS(sseur, "seur_Oligo_rep2_CB_noDoublet_noAstro_v2.rds")

