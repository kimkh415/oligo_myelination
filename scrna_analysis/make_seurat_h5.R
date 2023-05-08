library(Seurat)
library(stringr)
library(Matrix)
source("~/kwanho/src/readCB_h5.R")


#loads all 10X lanes from a given directort
dir10X<-function(h5.list,outdir="",makeSeurat=T,minGenes=500,regress=c("nCount_RNA", "percent.mt"), outname='CB')
{
        print(h5.list)
	mat = NULL

        for (i in 1:length(h5.list)) {
	  h5.file = h5.list[i]
	  print(h5.file)

	  samp.name = gsub("_filtered.h5", "", basename(h5.file))

	  dat = ReadCB_h5(h5.file)
	  cell.name = colnames(dat)
	  new.cell.name = paste0(samp.name, '_', cell.name)
	  colnames(dat) <- new.cell.name

          print(samp.name)
	  print(paste("Count matrix dims: ",toString(dim(dat))))
	  
	  if (is.null(mat)) {
	    mat = dat
	  } else {
	    mat = cbind(mat, dat)
	  }
	}

	print("Number of cells from each sample")
	print(table(str_split(colnames(mat), '_', simplify=T)[,1], str_split(colnames(mat), '_', simplify=T)[,2]))

	if(!makeSeurat){return(mat)}

        print("Make object!")
        # don't subset by minGenes yet. Must have all the cells to add another assay
	seur<-CreateSeuratObject(mat,"Seurat", min.features=minGenes)
	seur$age = str_split(colnames(seur), '_', simplify=T)[,1]
	seur$layer = str_split(colnames(seur), '_', simplify=T)[,2]
	seur$sample = paste0(seur$age, '_', seur$layer)
	print(table(seur$age, seur$layer))

        #QC from tutorial
        # The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.
        mito.genes <- grep(pattern = "^mt-", x = rownames(x = seur), value = TRUE)
        percent.mt <- 100 * Matrix::colSums(GetAssayData(seur, slot='counts')[mito.genes, ]) / Matrix::colSums(GetAssayData(seur, slot="counts"))
        seur[['percent.mt']] <- percent.mt
        ribo.genes <- grep("^Rp[s,l]",rownames(seur), value = TRUE)
        percent.rb <- 100 * Matrix::colSums(GetAssayData(seur, slot='counts')[ribo.genes, ]) / Matrix::colSums(GetAssayData(seur, slot="counts"))
        seur[['percent.rb']] <- percent.rb
        pdf(file.path(outdir,paste0('QC_plots_',outname,'.pdf')))
        par(mfrow = c(1, 3))
        print(FeatureScatter(object = seur, feature1 = "nCount_RNA", feature2 = "percent.mt"))
        print(FeatureScatter(object = seur, feature1 = "nCount_RNA", feature2 = "percent.rb"))
        print(FeatureScatter(object = seur, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
        dev.off()

        print("Normalize data!")
        seur = NormalizeData(seur, normalization.method = "LogNormalize")

        print("Get variable genes!")
        seur<-FindVariableFeatures(seur)

        print("Regress out!")
        seur<-ScaleData(seur,features=VariableFeatures(seur),vars.to.regress=regress)

        print("Run PCA!")
        seur<-RunPCA(seur,features = VariableFeatures(seur),verbose=F)
	
        print("Save object!")
        saveRDS(seur, file = file.path(outdir,paste0("seur_initial_",outname,".rds")))

        pdf(file.path(outdir,paste0('pca_plots_',outname,'.pdf')))
        print(VizDimLoadings(object = seur, dims = 1:2))
        print(DimPlot(object = seur, dims = c(1,2)))
        DimHeatmap(object = seur, dims = 1:9, cells = 500, balanced = TRUE)
        print(ElbowPlot(object = seur, 30))
        dev.off()
}


#args = commandArgs(trailingOnly=T)
#data.dir = args[1]
data.dir = "/stanley/levin_dr/kwanho/projects/vahbiz/analysis_20211210/new_data_analysis/cellbender"
h5.list = list.files(path=data.dir, pattern="filtered.h5", full.names=T)
outdir='./'

dir10X(h5.list, outdir=outdir, minGenes=500)

