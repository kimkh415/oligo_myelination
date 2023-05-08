library(Seurat)
library(SeuratDisk)
source("~/kwanho/src/seurat_tools.R")
#library(scales)

seur <- readRDS("../rca_merge_reps/seur_oligo_final_rca_v2.rds")
mol = readRDS("../random_forest/seur_MOL_subset_processed_final.rds")

cols.mol = c("#FFA500", "#FF5800", "#883000")

maturation.gl = list(
myelin=c('Mag','Cnp','Ptgds','Qdpr','Mal','Mobp'),
NFOL=c('Enpp6','Peli1','Map1b','Tmsb4x','Epb41l2','Bcas1')
)



mol = MyModuleScore(mol, maturation.gl, filename="metadata_MOL_module_score_maturation.rds")



pdf("supp_fig_3a_violin_maturation_MOL_subtypes_left.pdf")
myStackedVlnPlot(mol, features=maturation.gl[['myelin']], group.by='MOL_type', cols=cols.mol)
dev.off()

pdf("supp_fig_3a_violin_maturation_MOL_subtypes_right.pdf")
myStackedVlnPlot(mol, features=maturation.gl[['NFOL']], group.by='MOL_type', cols=cols.mol)
dev.off()


MyViolinPlot(mol, features=names(maturation.gl), group.by='MOL_type', cols=cols.mol, filename='supp_fig3a_violin_module_scores.pdf', hei=5, wid=3)



smol = subset(mol, age=="P14", invert=T)
smol$sample = droplevels(smol$sample)
PropsPlot(smol, 'MOL_type','sample',name.prefix='props_MOL_subtype_final', cols=cols.mol, hei=5, wid=8)



seur$final_name = as.character(seur$OligoType)
seur$final_name[colnames(mol)] = as.character(mol$MOL_type)
seur <- subset(seur, final_name=="MOL", invert=T)  # we removed samll MOL subcluster (52 cells) that had neuronal marker genes
seur$final_name = factor(seur$final_name, levels=c('OPC','COP','NFOL','MFOL','MOL_A','MOL_B','MOL_C'))


seur2 = CreateSeuratObject(counts=seur@assays$RNA@counts, meta.data=seur@meta.data[,c('rep_age','replicate','age','layer','sample','OligoType','final_name')])
out.file = "ref_OL_types.h5Seurat"
SaveH5Seurat(seur2, filename=out.file)
Convert(out.file, dest='h5ad')

print("DONE")



