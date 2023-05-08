library(data.table)
source("~/kwanho/src/seurat_tools.R")


mat = as.matrix(fread('Puck_200306_02.digital_expression.txt'), rownames=1)

coords = read.csv("BeadLocationsForR copy.csv", row.names=1)
colnames(coords) = c('spatial_1','spatial_2')

seur <- CreateSeuratObject(mat, project='slide-seq_puck_200306_02')

seur[['spatial']] = CreateDimReducObject(embeddings=as.matrix(coords), key='spatial_')

seur <- NormalizeData(seur)

plot_feature2(seur, alpha=0.5, reduction='spatial', features="nCount_RNA", filename='featurplot_nUMI.pdf', cols=rev(magma(50)), size=6)

genes = c('Neurod2','Slc17a7','Cux1','Cux2','Satb2','Rorb','Fezf2','Bcl11b','Tcerg1l','Tle4','Olig1','Olig2','Plp1','Pdgfra','Mbp')
plot_feature2(seur, alpha=0.5, reduction='spatial', features=genes, filename='featurplot_genes.pdf', cols=rev(magma(50)), size=4)

# rotate puck (45 degrees = .785 rad)
library(spdep)
rot_coords = Rotation(coords, angle=.6)
colnames(rot_coords) = c('x', 'y')
write.csv(rot_coords, "data/rotated_coordinates.csv", quote=F)

colnames(rot_coords) = c('rotated_1','rotated_2')
seur[['spatial_rot']] = CreateDimReducObject(embeddings=as.matrix(rot_coords), key='rotated_')

plot_feature2(seur, alpha=0.5, reduction='spatial_rot', features="nCount_RNA", filename='featurplot_nUMI_rotated.pdf', cols=rev(magma(50)), size=6)

plot_feature2(seur, alpha=0.5, reduction='spatial_rot', features=genes, filename='featurplot_genes_rotated.pdf', cols=rev(magma(50)), size=4)

p <- FeaturePlot(seur, reduction='spatial_rot', features='Slc17a7', cols=rev(magma(50)))
pdf("slice_selection_in_puck.pdf")
p
p +	geom_rect(aes(xmin=250, xmax=1500, ymin=4500, ymax=6900), color='red', fill='transparent')
dev.off()

saveRDS(seur, "seur_init_rotated.rds")

