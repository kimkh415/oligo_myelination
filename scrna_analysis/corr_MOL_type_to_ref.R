library(Seurat)
library(corrplot)
library(dplyr)


ref = readRDS("/stanley/levin_dr/kwanho/projects/vahbiz/redo_all/data/ref_Branco_cortex_only_MOL.rds")
Idents(ref) <- 'cell_type'
mols <- paste0("MOL", 1:6)
Idents(ref) <- ordered(Idents(ref), levels=mols)

ref = subset(ref, idents='MOL2', invert=T)

#deg.dir = "../../DEG/MOL_subtype/"
#deg.list = list.files(path=deg.dir, pattern='.rds', full.names=T)
#degs = list()
#for (i in 1:length(deg.list)) {
#res = readRDS(deg.list[i]) %>% filter(padj < 0.05 & Estimate > .3) %>% arrange(desc(Estimate))
#print(head(res))
#degs[[i]] = res$Gene
#}
#un = unname(unlist(degs))
#dup.genes = un[which(duplicated(un))]
#
#for (i in 1:length(degs)) {
#degs[[i]] = setdiff(degs[[i]], dup.genes)
#}
#
#all.genes = unique(unname(unlist(degs)))
#length(all.genes)

query = readRDS("../seur_MOL_subset_processed_final.rds")
#query = subset(query, idents='4', invert=T)

ref@active.assay='RNA'
query@active.assay='RNA'

ref.ae = AverageExpression(ref)
query.ae = AverageExpression(query)

#query <- FindVariableFeatures(query, assay='RNA')
#ref <- FindVariableFeatures(ref, assay='RNA')
query.vg = VariableFeatures(query)
ref.vg = VariableFeatures(ref)
all.vg = unique(c(ref.vg, query.vg))
#common.genes = intersect(rownames(query), ref.vg)
#common.genes = intersect(rownames(query), intersect(rownames(ref), all.vg))
#common.genes = intersect(all.genes, rownames(ref))
#qae = query.ae$RNA[common.genes,]
#rae = ref.ae$RNA[common.genes,]
#levels(query)
#levels(ref)
#cor.mat = cor(rae, qae, method="spearman")
#saveRDS(cor.mat, "corr_mat_degs.rds")
#
#pdf("corrplot_pie_degs.pdf", height=5, width=5)
#corrplot(t(cor.mat), method='pie', tl.col="black", col=colorRampPalette(c("gray","gray"))(200), cl.pos='n', tl.srt=90, cl.length=5, addCoef.col='black', number.cex=.7)
#dev.off()

#common.genes = intersect(rownames(query), ref.vg)
common.genes = intersect(rownames(ref), query.vg)
#common.genes = intersect(rownames(ref), common.genes)
print(length(common.genes))
qae = query.ae$RNA[common.genes,]
rae = ref.ae$RNA[common.genes,]
levels(query)
levels(ref)
cor.mat = cor(rae, qae, method="spearman")
saveRDS(cor.mat, "corr_mat_ref_var_genes_final.rds")

pdf("corrplot_pie_ref_var_genes_final.pdf", height=5, width=5)
corrplot(t(cor.mat), method='pie', tl.col="black", col=colorRampPalette(c("gray","gray"))(200), cl.pos='n', tl.srt=90, cl.length=5, addCoef.col='black', number.cex=.7)
dev.off()


deg.a = readRDS("../../DEG/MOL_subtype/nb_res_MOL_clusterMOL_A.rds")
deg.a$cluster = "MOL_A"
deg.b = readRDS("../../DEG/MOL_subtype/nb_res_MOL_clusterMOL_B.rds")
deg.b$cluster = "MOL_B"
deg.c = readRDS("../../DEG/MOL_subtype/nb_res_MOL_clusterMOL_C.rds")
deg.c$cluster = "MOL_C"
deg.all = do.call(rbind, list(deg.a, deg.b, deg.c)) %>% filter(padj<0.05 & Estimate>0.5)
print(table(deg.all$cluster))
deg = deg.all %>% group_by(cluster) %>% top_n(n=200, wt=Estimate) %>% distinct(Gene) %>% pull(Gene)
#top100 = deg %>% group_by(cluster) %>% top_n(n=100, wt=Estimate)

common.genes = intersect(rownames(ref), deg)
print(length(common.genes))
qae = query.ae$RNA[common.genes,]
rae = ref.ae$RNA[common.genes,]
levels(query)
levels(ref)
cor.mat = cor(rae, qae, method="spearman")
saveRDS(cor.mat, "corr_mat_DE_genes_final.rds")

pdf("corrplot_pie_DE_genes_final.pdf", height=5, width=5)
corrplot(t(cor.mat), method='pie', tl.col="black", col=colorRampPalette(c("gray","gray"))(200), cl.pos='n', tl.srt=90, cl.length=5, addCoef.col='black', number.cex=.7)
dev.off()


