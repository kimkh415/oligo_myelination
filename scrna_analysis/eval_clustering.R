library(Seurat)
library(dplyr)
library(scSHC)
source("~/kwanho/src/seurat_tools.R")




mol.file = paste0("seur_MOL_subset_processed_final.rds")
seur <- readRDS(mol.file)


cls.res = 0.3
print(cls.res)
seur <- FindClusters(seur, resolution=cls.res)
pdf(paste0("umap_mol_clusters_res",cls.res,".pdf"))
MyDimPlot(seur, alpha=0.7)
dev.off()

dat = seur@assays$RNA@counts[VariableFeatures(seur),]
res <- testClusters(dat, as.character(Idents(seur)), batch=as.character(seur$rep_age), num_features=nrow(dat), num_PCs=20, cores=8)
saveRDS(res, paste0("scshs_cluster_sig_res_", cls.res, ".rds"))
print(table(res[[1]], Idents(seur)))
print(res[[2]])

seur$scshs_clusters = res[[1]]
pdf(paste0("umap_mol_scshs_clusters_res",cls.res,".pdf"))
MyDimPlot(seur, alpha=0.7, group.by='scshs_clusters')
dev.off()

sseur <- subset(seur, age=='P14', invert=T)
sseur$sample = droplevels(sseur$sample)
PropsPlot(sseur, my.group='seurat_clusters', my.sample='sample',name.prefix=paste0('barplot_props_MOL_seur_res',cls.res))
PropsPlot(sseur, my.group='scshs_clusters', my.sample='sample',name.prefix=paste0('barplot_props_MOL_scshs_res',cls.res))





cls.res = 1
print(cls.res)
seur <- FindClusters(seur, resolution=cls.res)
pdf(paste0("umap_mol_clusters_res",cls.res,".pdf"))
MyDimPlot(seur, alpha=0.7)
dev.off()

dat = seur@assays$RNA@counts[VariableFeatures(seur),]
res <- testClusters(dat, as.character(Idents(seur)), batch=as.character(seur$rep_age), num_features=nrow(dat), num_PCs=20, cores=8)
saveRDS(res, paste0("scshs_cluster_sig_res_", cls.res, ".rds"))
print(table(res[[1]], Idents(seur)))
print(res[[2]])

seur$scshs_clusters = res[[1]]
pdf(paste0("umap_mol_scshs_clusters_res",cls.res,".pdf"))
MyDimPlot(seur, alpha=0.7, group.by='scshs_clusters')
dev.off()

sseur <- subset(seur, age=='P14', invert=T)
sseur$sample = droplevels(sseur$sample)
PropsPlot(sseur, my.group='seurat_clusters', my.sample='sample',name.prefix=paste0('barplot_props_MOL_seur_res',cls.res))
PropsPlot(sseur, my.group='scshs_clusters', my.sample='sample',name.prefix=paste0('barplot_props_MOL_scshs_res',cls.res))


# Try feeding it harmony


print('DONE')

