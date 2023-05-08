# r_4
# plot in pdf - easier to work with in AI


library(Seurat)
library(ggplot2)
source("~/kwanho/src/seurat_tools.R")

# Colors
cols.age = c("#FFA07A", "#FF0000" ,"#A63A7A", "#67349C")
names(cols.age) = c('P7','P14','P30','P90')
cols.age.rep2 = cols.age[c('P30','P90')]

cols.feature = colorRampPalette(c("lightgray", "red"))(50)

# Modules
gl = list(
astro=c('Aldoc','Aldh1l1','S100b','Sox9'),
neuronal=c('Slc17a7','Neurod2','Grin2b','Grin1','Gad2','Gad1','Pax6','Pvalb','Sst','Npy'),
opc=c("Cspg5","Fabp7","Pdgfra","Cdo1","Ptprz1","Ednrb","Serpine2","Meg3","Lhfpl3","Pcp4","Mt3","Rlbp1","Ddah1","Cntn1","Rpl13","Zfp36l1","Nnat","Sparcl1","Gpm6a","Atp1a2"),
cop=c("Fyn","Epb41l2","Bcas1","Map1b","Lims2","Tmsb4x","Frmd4a","Pdcd4","Dynll2","Sirt2","Pfn2","Chn2","Enpp6","Tns3","Bmp4","Mpzl1","Gpr17","Mycl","Slc44a1","Cyfip2"),
nfol=c("Mbp","Gsto1","Nfasc","Reep5","Arpc1b","Cpox","Rras2","Mobp","Cldn11","Peli1","Gm15440","Tmem141","Secisbp2l","Idh1","Tcf7l2","Elovl7","Tmem2","Tppp","Ninj1"),
mfol=c("Ctsl","Mobp","Opalin","Tspan2","Mog","Cldn11","H2afj","Tmem141","Hist1h2bc","Pdlim2","Mag","Hist1h1c","Mal","Grb14","Ppp1r14a","Birc2","Gamt","Ugt8a","Plp1"),
mol=c("Arpc1b","Tmeff2","Ptgds","Trf","Car2","Sept4","Mal","Cryab","Csrp1","Glul","Apod","Qdpr","Aplp1","Ndrg1","Ermn","Fth1","Plp1","Ppp1r14a","Mog","Mag","Pex5l","Cldn11")
)




#########################################
## Rep2 Pre QC
#########################################
# rep2 after cellbender before filtering cells Seurat object
seur <- readRDS("/stanley/levin_dr/kwanho/projects/vahbiz/analysis_20211210/new_data_analysis/seurat/post_QC/seur_initial_CB.rds")

seur <- RunUMAP(seur, dims=1:30)
mito.genes <- grep("^mt-", rownames(seur), value=T)
percent.mito <- 100 * Matrix::colSums(GetAssayData(seur, slot='counts')[mito.genes, ]) / Matrix::colSums(GetAssayData(seur, slot="counts"))
seur[['percent.mito']] <- percent.mito
ribo.genes <- grep("Rp[sl]", rownames(seur), value=T)
percent.ribo <- 100 * Matrix::colSums(GetAssayData(seur, slot='counts')[ribo.genes, ]) / Matrix::colSums(GetAssayData(seur, slot="counts"))
seur[['percent.ribo']] <- percent.ribo


# supp 1b umap age
pdf("supp_1b_rep2_umap_preQC_ages_17052_cells.pdf", height=6, width=5)
MyDimPlot(seur, alpha=1, group.by='age', cols=cols.age.rep2, legend.pos='right')
dev.off()

# supp 1c pre-QC violin - nGene, nUMI, percent.mito, percent.ribo
pdf("supp_1c_rep2_violin_preQC.pdf", height=6, width=12)
VlnPlot(seur, group.by='age', cols=cols.age.rep2, features=c('nFeature_RNA','nCount_RNA','percent.mito','percent.ribo'), ncol=4, pt.size=0)
dev.off()

seur <- MyModuleScore(seur, gl, filename="metadata_rep2_preQC_cell_type_module_score.rds")

# supp 1def featureplot neuronal, astro and oligo markers (gray to red colors)
plot_feature2(seur, features=names(gl), cols=cols.feature, alpha=1, size=3, filename="supp_1def_rep2_preQC_cell_type_module_score.pdf", min.cutoff='q10')

# removed 
# 134 doublets
# 502 from QC metrics
# 109 Astrocytes


#########################################
## Rep2 Post QC
#########################################
# rep2 post QC Seurat object
seur <- readRDS("/stanley/levin_dr/kwanho/projects/vahbiz/analysis_20211210/new_data_analysis/seurat/post_QC/seur_Oligo_rep2_CB_noDoublet_noAstro_v2.rds")



# supp 1g umap age
cols.age.rep2 = cols.age[c('P30','P90')]
pdf("supp_1g_rep2_umap_preQC_ages_16307_cells.pdf", height=6, width=5)
MyDimPlot(seur, alpha=1, group.by='age', cols=cols.age.rep2, legend.pos='right')
dev.off()

# supp 1h pre-QC violin - nGene, nUMI, percent.mito, percent.ribo
pdf("supp_1h_rep2_violin_preQC.pdf", height=6, width=12)
VlnPlot(seur, group.by='age', cols=cols.age.rep2, features=c('nFeature_RNA','nCount_RNA','percent.mt','percent.rb'), ncol=4, pt.size=0)
dev.off()

seur <- MyModuleScore(seur, gl, filename="metadata_rep2_postQC_cell_type_module_score.rds")

# supp 1i - module score feature plot post qC
plot_feature2(seur, features=c('opc','cop','nfol','mfol','mol'), cols=cols.feature, alpha=1, size=3, filename="supp_1i_rep2_postQC_cell_type_module_score.pdf", min.cutoff='q10')













