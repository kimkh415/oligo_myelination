# run in r_4

library(Seurat)
library(ggplot2)
library(dplyr)
library(Hmisc)
source("~/kwanho/src/dimReduce_new.r")
source("~/kwanho/src/Glmer.R")
source("~/kwanho/src/seurat_tools.R")


seur <- readRDS("seur_Oligo_merged_processed_v2.rds")
ref <- readRDS("../../../data/processed_data/ProjectOnto.RDS")

# save initial clustering before comparing to the reference
seur$clusts_before_rca = Idents(seur)

# perform Sean's implementation of RCA-based analysis
seur <- RCA_Seurat(seur, ref)
print(table(Idents(seur)))

seur <- FindNeighbors(seur, reduction='cor', dims=1:dim(seur@reductions$cor@cell.embeddings)[2])
seur <- FindClusters(seur, resolution=0.5)
seur = RunUMAP(seur, reduction='cor', dims=1:dim(seur@reductions$cor@cell.embeddings)[2])

oligo.col = c(OPC="#A6CEE3",COP="mediumorchid3",NFO="#028F1E",MFO="#FF474C",MOL="#FFB16D")
seur$Corr_type = factor(seur$Corr_type, levels=names(oligo.col))

seur <- NormalizeData(seur, assay='RNA')

saveRDS(seur, "seur_oligo_merged_rca_v2.rds")

pdf("vis_rca_initial_v2.pdf")
MyDimPlot(seur, 1, reduction='umap', group.by='Corr_type', label=F, cols=oligo.col)
FeaturePlot(seur, reduction='umap', features='max_Cor')
MyDimPlot(seur, 1, reduction='umap', group.by='seurat_clusters', label=T)
#MyDimPlot(seur, 1, reduction='umap', group.by='sample', label=F)
MyDimPlot(seur, 1, reduction='umap', group.by='replicate', label=F)
MyDimPlot(seur, 1, reduction='umap', group.by='age', label=F)
#MyDimPlot(seur, 1, reduction='umap', group.by='layer', label=F)
dev.off()

oligo.genes = c('Pdgfra','Bmp4','Tmem2','Opalin','Ptgds','Plp1')
pdf("violin_genes_rca_initial_v2.pdf")
myStackedVlnPlot(seur, oligo.genes, assay='RNA')
dev.off()

plot_feature2(seur, oligo.genes, reduction='umap', nc=2, size=3, filename='featureplot_oligo_genes_v2.pdf')

# annotation
anno = c('OPC','OPC','OPC','MOL','MOL','MOL','OPC','subcluster','subcluster','OPC','MOL','subcluster','subcluster')
names(anno) = 0:12

seur <- RenameIdents(seur, anno)


# refine annotation by subclustering
sseur <- subset(seur, idents='subcluster')
sseur <- sseur %>%
        FindNeighbors(sseur, reduction='cor', dims=1:dim(seur@reductions$cor@cell.embeddings)[2]) %>%
        RunUMAP(reduction='cor', dims=1:dim(seur@reductions$cor@cell.embeddings)[2]) %>%
        FindClusters(resolution=0.1)

pdf("umap_rca_v2_subcluster.pdf")
MyDimPlot(sseur, 1, reduction='umap', group.by='seurat_clusters', label=T)
MyDimPlot(sseur, 1, reduction='umap', group.by='replicate')
MyDimPlot(sseur, 1, reduction='umap', group.by='sample')
dev.off()

pdf("violin_genes_rca_initial_v2_subcluster.pdf")
myStackedVlnPlot(sseur, oligo.genes, assay='RNA')
dev.off()

plot_feature2(sseur, oligo.genes, reduction='umap', nc=2, size=3, filename='featureplot_oligo_genes_v2_subcluster.pdf')

sub.anno = c('COP','MOL','MFOL','NFOL')
names(sub.anno) = 0:3


# Merge annotation
oligo.col = c(OPC="#A6CEE3",COP="mediumorchid3",NFOL="#028F1E",MFOL="#FF474C",MOL="#FFB16D")

seur$OligoType = as.character(Idents(seur))
seur$OligoType[colnames(sseur)] = as.character(Idents(sseur))
seur$OligoType = factor(seur$OligoType, levels=c(names(oligo.col)))

Idents(seur) = 'OligoType'

pdf("umap_rca_OligoType_v2.pdf")
MyDimPlot(seur, 1, reduction='umap', group.by='OligoType', label=F, cols=oligo.col)
MyDimPlot(seur, 1, reduction='umap', group.by='replicate', label=F)
MyDimPlot(seur, 1, reduction='umap', group.by='sample', label=F)
#MyDimPlot(seur, 1, reduction='umap', group.by='age', label=F)
#MyDimPlot(seur, 1, reduction='umap', group.by='layer', label=F)
dev.off()

pdf("supp_fig_S2B_v2.pdf")
myStackedVlnPlot(seur, oligo.genes, group.by='OligoType', assay='RNA', cols=oligo.col)
dev.off()

# supp fig 2c
pdf("supp_fig_S2C_age_v2.pdf", height=10, width=10)
DimPlot(seur, reduction='umap', group.by='OligoType', split.by='age', nc=2, cols=oligo.col)
dev.off()

# supp fig 2c
pdf("supp_fig_S2C_rep_v2.pdf", height=5, width=10)
DimPlot(seur, reduction='umap', group.by='OligoType', split.by='replicate',nc=2, cols=oligo.col)
dev.off()

# main fig 1b
pdf("main_fig_1B.pdf")
MyDimPlot(seur,.8, group.by='OligoType', cols=oligo.col, pt.size=.5)
dev.off()

seur$rep_age = paste0(seur$replicate, '_', seur$age)
seur$rep_age = factor(seur$rep_age, levels=c('rep1_P7','rep1_P14','rep1_P30','rep2_P30','rep1_P90','rep2_P90'))

pdf("main_fig_1C_top_v2.pdf", height=5, width=18)
MyDimPlot(seur,.8, group.by='OligoType', split.by='rep_age', nc=6, cols=oligo.col, pt.size=.5)
dev.off()

# supp fig 7d
seur$layer = gsub("II-III", "II-IV", seur$layer)
seur$sample = paste0(seur$replicate, '_', seur$age, '_', seur$layer)
seur$sample = factor(seur$sample, levels=c('rep1_P7_LayerII-IV','rep1_P7_LayerV','rep1_P7_LayerVI','rep1_P14_LayerII-IV','rep1_P14_LayerV','rep1_P14_LayerVI','rep1_P30_LayerII-IV','rep1_P30_LayerV','rep1_P30_LayerVI','rep2_P30_LayerII-IV','rep2_P30_LayerV','rep2_P30_LayerVI','rep1_P90_LayerII-IV','rep1_P90_LayerV','rep1_P90_LayerVI','rep2_P90_LayerII-IV','rep2_P90_LayerV','rep2_P90_LayerVI'))
seur = MyModuleScore(seur, gene.list=list(DL_Oligo=c("Lgr4", "Fgfr1", "Fgfr2", "Tyro3")), filename="metadata_module_score_DL_OL_receptor_4_genes_v2.rds")
MyViolinPlot(seur, group.by='sample', features='DL_Oligo', filename="supp_fig7d.pdf", height.multiplier=1.5)

# props plot
PropsPlot(seur, 'OligoType','sample', name.prefix="main_fig_1C_barplot_props_OligoType", cols=oligo.col, hei=5, wid=8)

# Oligo module exp
gl = list(
opc=c("Cspg5","Fabp7","Pdgfra","Cdo1","Ptprz1","Ednrb","Serpine2","Meg3","Lhfpl3","Pcp4","Mt3","Rlbp1","Ddah1","Cntn1","Rpl13","Zfp36l1","Nnat","Sparcl1","Gpm6a","Atp1a2"),
cop=c("Fyn","Epb41l2","Bcas1","Map1b","Lims2","Tmsb4x","Frmd4a","Pdcd4","Dynll2","Sirt2","Pfn2","Chn2","Enpp6","Tns3","Bmp4","Mpzl1","Gpr17","Mycl","Slc44a1","Cyfip2"),
nfol=c("Mbp","Gsto1","Nfasc","Reep5","Arpc1b","Cpox","Rras2","Mobp","Cldn11","Peli1","Gm15440","Tmem141","Secisbp2l","Idh1","Tcf7l2","Elovl7","Tmem2","Tppp","Ninj1"),
mfol=c("Ctsl","Mobp","Opalin","Tspan2","Mog","Cldn11","H2afj","Tmem141","Hist1h2bc","Pdlim2","Mag","Hist1h1c","Mal","Grb14","Ppp1r14a","Birc2","Gamt","Ugt8a","Plp1"),
mol=c("Arpc1b","Tmeff2","Ptgds","Trf","Car2","Sept4","Mal","Cryab","Csrp1","Glul","Apod","Qdpr","Aplp1","Ndrg1","Ermn","Fth1","Plp1","Ppp1r14a","Mog","Mag","Pex5l","Cldn11")
)

seur <- MyModuleScore(seur, gl, filename="metadata_OligoType_module_score.rds")

# supp 2e
cols.feature = colorRampPalette(c("lightgray", "red"))(50)
plot_feature2(seur, features=names(gl), cols=cols.feature, alpha=1, size=3, filename="supp_2e_OligoType_module_score.pdf", min.cutoff='q10')


saveRDS(seur$OligoType, "metadata_OligoType_v2.rds")
saveRDS(seur, "seur_oligo_final_rca_v2.rds")



