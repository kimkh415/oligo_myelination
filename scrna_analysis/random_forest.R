library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(pheatmap)
library(dplyr)
library(harmony)
source("~/kwanho/src/RF_utils.R")
source("~/kwanho/src/seurat_tools.R")


##############################################################################
# Prepare data
print("Loading data")
ref <- readRDS("/stanley/levin_dr/kwanho/projects/vahbiz/redo_all/data/ref_Branco_cortex_only_MOL.rds")

Idents(ref) <- 'cell_type'
mols <- paste0("MOL", 1:6)
Idents(ref) <- ordered(Idents(ref), levels=mols)
ref <- subset(ref, idents='MOL2', invert=T)

v='final'  # final
mol.file = paste0("seur_MOL_subset_processed_",v,".rds")
if (!file.exists(mol.file)) {
seur <- readRDS("/stanley/levin_dr/kwanho/projects/vahbiz/analysis_20211210/new_data_analysis/seurat/rca_merge_reps/seur_oligo_final_rca_v2.rds")

seur <- subset(seur, idents='MOL')
seur@active.assay = 'RNA'
seur <- seur %>%
	#SCTransform(vars.to.regress=c('percent.mt','nFeature_RNA','rep_age','layer')) %>%
	NormalizeData() %>%
	FindVariableFeatures() %>%
	ScaleData(vars.to.regress=c('nCount_RNA','rep_age','layer')) %>%
	RunPCA() %>%
	RunHarmony('age') %>%
	FindNeighbors(reduction = "harmony", dims=1:20) %>%
	FindClusters(resolution=0.1) %>%
	RunUMAP(reduction = "harmony", dims=1:20)

print(table(Idents(seur), seur$sample))
#s.genes <- cc.genes$s.genes
#g2m.genes <- cc.genes$g2m.genes
#s.genes = human2mouse(s.genes)
#g2m.genes = human2mouse(g2m.genes)
#
#pdf("umap_mol_CC.pdf")
#DimPlot(seur, group.by='Phase')
#FeaturePlot(seur, features='S.Score')
#FeaturePlot(seur, features='G2M.Score')
#dev.off()

pdf(paste0("umap_mol_subset.pdf"), height=7, width=8)
DimPlot(seur, group.by=c('seurat_clusters', 'replicate', 'age', 'layer'), nc=2)
FeaturePlot(seur, 'nCount_RNA')
FeaturePlot(seur, 'nFeature_RNA')
dev.off()

#meta = readRDS("../../../our_clusters/metadata_our_3_clusters.rds")
#rep1 = subset(seur, replicate=='rep1')
#meta=meta[intersect(names(meta), colnames(rep1))]
#rep1$prev_MOL_subtype = "MOL"
#rep1$prev_MOL_subtype[names(meta)] = as.character(meta)
#print(table(rep1$prev_MOL_subtype))
#pdf(paste0("umap_rep1_prev_mol_subtype_v_new.pdf"))
#DimPlot(rep1, group.by='prev_MOL_subtype') + ggtitle("Old MOL subtypes")
#DimPlot(rep1) + ggtitle("New MOL clusters")
#dev.off()

#sseur <- subset(seur, age=='P14', invert=T)
#sseur$sample = droplevels(sseur$sample)
#PropsPlot(sseur, my.group='seurat_clusters', my.sample='sample',name.prefix=paste0('barplot_props_MOL_subtypes'))

cls.res = 0.3
print(cls.res)
seur <- FindClusters(seur, resolution=cls.res)
pdf(paste0("umap_mol_clusters_res",cls.res,".pdf"))
MyDimPlot(seur, alpha=0.7)
dev.off()

seur <- subset(seur, idents='6', invert=T)
print(table(Idents(seur)))

dat = seur@assays$RNA@counts[VariableFeatures(seur),]
res <- testClusters(dat, as.character(Idents(seur)), batch=as.character(seur$rep_age), num_features=nrow(dat), num_PCs=20, cores=8)
saveRDS(res, paste0("scshs_cluster_sig_res_", cls.res, ".rds"))
print(table(res[[1]], Idents(seur)))
print(res[[2]])

seur$scshs_clusters = res[[1]]
pdf(paste0("umap_mol_scshs_clusters_res",cls.res,".pdf"))
MyDimPlot(seur, alpha=0.7, group.by='scshs_clusters')
dev.off()

seur$MOL_type = res[[1]]
seur$MOL_type = factor(seur$MOL_type, levels=c('new1','new2','new3'))
levels(seur$MOL_type) = paste0('MOL_', LETTERS[1:3])
saveRDS(seur$MOL_type, "metadata_MOL_subtype.rds")

cols.mol = c("#FFA500", "#FF5800", "#883000")
pdf("main_fig_1f_umap_MOL_subtype_15505_cells.pdf")
MyDimPlot(seur, alpha=0.7, group.by='MOL_type', cols=cols.mol, pt.size=.5)
dev.off()


sseur <- subset(seur, age=='P14', invert=T)
sseur$sample = droplevels(sseur$sample)
#PropsPlot(sseur, my.group='seurat_clusters', my.sample='sample',name.prefix=paste0('barplot_props_MOL_seur_res',cls.res))
PropsPlot(sseur, my.group='scshs_clusters', my.sample='sample',name.prefix=paste0('barplot_props_MOL_scshs_res',cls.res))



saveRDS(seur, mol.file)

} else {
seur <- readRDS(mol.file)
}

cat("Reference data clusters:\n")
print(table(Idents(ref)))
cat("MOL data current clusters:\n")
print(table(Idents(seur)))

seur@active.assay = 'RNA'
ref@active.assay = 'RNA'
#ref <- ref %>% NormalizeData() %>% FindVariableFeatures()
#seur <- seur %>% NormalizeData() %>% FindVariableFeatures()

#all_genes <- intersect(rownames(ref), rownames(seur))
#all_genes <- all_genes[-grep("^Rp[ls]", all_genes, value=F)]
comb_var_genes = unique(c(VariableFeatures(ref), VariableFeatures(seur)))
genes <- intersect(comb_var_genes, rownames(ref))
genes <- intersect(genes, rownames(seur))

# Sample cells to have equal number of cells for each cluster
num <- min(table(Idents(ref)))
cat(paste0("Number of cells to be sampled from each cluster: ", num, "\n"))
ref<-subset(ref,downsample=num)

# Save sampled cells
write.table(colnames(ref), file="sampled_cells.tsv", sep='\t', quote=F, row.names=F, col.names=F)

##############################################################################
# For plotting
#cols <- brewer.pal(nlevels(ref), "Set1")
cols = rev(viridis(5))
names(cols) <- levels(ref)
cols["MOL1"] <- "#FEC20C"

dat1 <- data.frame(t(as.matrix(ref@assays$RNA@data[genes,])))
dat1["CellType"]=ref@active.ident
test.dat <- data.frame(t(as.matrix(seur@assays$RNA@data[genes,])))

# Split training data for validation
train.size <- floor(0.8*nrow(dat1))
sampled.rows <- sample(seq_len(nrow(dat1)), size=train.size)
train.dat <- dat1[sampled.rows,]
val.dat <- dat1[-sampled.rows,]
print("Training dataset size")
print(dim(train.dat))
print("Validation dataset size")
print(dim(val.dat))
print("Test dataset size")
print(dim(test.dat))

null.train <- train.dat
labels <- levels(train.dat$CellType)
nLabel <- nlevels(train.dat$CellType)
sampled.labels <- sample(1:nLabel, nrow(train.dat), replace=T)
null.train$CellType <- as.factor(labels[sampled.labels])

prefix='rerun'

# Training
print("Training!")
model <- fitMethod(train.dat)
saveRDS(model, paste0(prefix, "_RF_model.rds"))

# Training null model
print("Null model training!")
null_model <- fitMethod(null.train)
saveRDS(null_model, paste0(prefix, "_null_model.rds"))

# Validation
print("validation!")
#colors <- readRDS("~/microglia/fezf2_wt_ko/ml_DEG_without_Cd63/heatmap_colors.rds")
valMethod(model, val.dat, cols, model_name="RF")
valMethod(null_model, val.dat, cols, model_name="null")

# Test
print("Test!")
res <- testMethod(model, test.dat)
out <- res[[1]]
prob <- res[[2]]
saveRDS(out, paste0(prefix, "_prediction_output.rds"))
saveRDS(prob, paste0(prefix, "_prediction_probs.rds"))

seur$RF_prediction <- out
pdf(paste0(prefix, "_prediction_result.pdf"))
DimPlot(seur, reduction='umap', pt.size=.1, group.by="RF_prediction", label=TRUE, cols=cols) + NoAxes() + ggtitle("Prediction")
DimPlot(seur, reduction='umap', pt.size=.1, label=TRUE) + NoAxes() + ggtitle("Original")
dev.off()

plotRes(out, prob, cols, prefix)

pdf("supp_fig_3f_uamp_RF_prediction_15505_cells.pdf")
MyDimPlot(seur, alpha=.8, group.by='RF_prediction', cols=cols, legend.row=1, pt.size=.5)
dev.off()

saveRDS(seur, mol.file)

##############################################################################

print("DONE")

