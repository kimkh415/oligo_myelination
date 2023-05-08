library(MASS)
library(Seurat)
library(openxlsx)
source("~/kwanho/src/Glmer.R")

args = commandArgs(trailingOnly=T)

if (args[1] == 'OL_layer') {
print(args[1])

dir.create(args[1])

olg <- readRDS("../rca_merge_reps/seur_oligo_final_rca_v2.rds")

print(table(Idents(olg)))

res=list()

for(age in levels(olg$age))
{
  print(age)
  seur=subset(olg,cells=colnames(olg)[olg$age==age])
  seur@meta.data["RIBO"]=log(seur@meta.data[,"percent.rb"]+1)
  seur@meta.data["MT"]=log(seur@meta.data[,"percent.mt"]+1)

  lst<-summary(Idents(seur))
  lst<-names(lst)[lst>100]

  for(i in lst)
  {
    print(i)
    for(j in names(table(seur$layer)))
    {
      print(j)
      res.file = paste0(args[1],"/nb_res_",age,'_',i,'_',j,".rds")
      if (!file.exists(res.file)) {
        seur@meta.data["ko_vs_wt"]="k"
        seur@meta.data[seur@meta.data$layer==j,"ko_vs_wt"]="w"
        mrk<-TestGlmNB(seur,i,formula1=Expression~offset(logUMI)+ko_vs_wt+MT+RIBO)
        saveRDS(mrk, res.file)
      } else {
        mrk <- readRDS(res.file)
	res[[paste(age,i,j,sep="_")]]=mrk
      }

      sig.mrk = mrk[mrk$padj<0.05,]
      write.xlsx(sig.mrk, file=paste0(args[1],'/nb_res_sig_only_',age,'_',i,'_',j,'.xlsx'), overwrite=T)
    }
  }
}

saveRDS(res,file="OligoType.layer.NB.rds")

}


if (args[1] == 'MOL_subtype') {
print(args[1])

dir.create(args[1])

mol <- readRDS("../random_forest/seur_MOL_subset_processed_final.rds")
mol@meta.data["RIBO"]=log(mol@meta.data[,"percent.rb"]+1)
mol@meta.data["MT"]=log(mol@meta.data[,"percent.mt"]+1)

print(table(Idents(mol)))

res=list()

for (x in levels(mol)) {
  print(x)
  res.file = paste0(args[1], "/nb_res_MOL_cluster", x, ".rds")
  if (!file.exists(res.file)) {
    mol@meta.data['ko_vs_wt'] = 'k'
    mol@meta.data[Idents(mol)==x, 'ko_vs_wt'] = 'w'
    mrk<-TestGlmNB(mol,levels(mol),formula1=Expression~offset(logUMI)+ko_vs_wt+MT+RIBO)
    saveRDS(mrk, res.file)
  } else {
    mrk <- readRDS(res.file)
    res[[paste0("cluster", x)]] = mrk
  }
  sig.mrk = mrk[mrk$padj<0.05,]
  write.xlsx(sig.mrk, file=paste0(args[1],'/nb_res_sig_only_MOL_cluster', x, '.xlsx'), overwrite=T)
}

saveRDS(res, file="MOL_cluster_NB.rds")

}


if (args[1] == 'OL_type') {
print(args[1])

dir.create(args[1])

olg <- readRDS("../rca_merge_reps/seur_oligo_final_rca_v2.rds")
olg@meta.data["RIBO"]=log(olg@meta.data[,"percent.rb"]+1)
olg@meta.data["MT"]=log(olg@meta.data[,"percent.mt"]+1)

print(table(Idents(olg)))

res=list()

for (x in levels(olg)) {
  print(x)
  res.file = paste0(args[1], "/nb_res_OL_type", x, ".rds")
  if (!file.exists(res.file)) {
    olg@meta.data['ko_vs_wt'] = 'k'
    olg@meta.data[Idents(olg)==x, 'ko_vs_wt'] = 'w'
    mrk<-TestGlmNB(olg,levels(olg),formula1=Expression~offset(logUMI)+ko_vs_wt+MT+RIBO)
    saveRDS(mrk, res.file)
  } else {
    mrk <- readRDS(res.file)
    res[[paste0("cluster", x)]] = mrk
  }
  sig.mrk = mrk[mrk$padj<0.05,]
  write.xlsx(sig.mrk, file=paste0(args[1],'/nb_res_sig_only_OL_type', x, '.xlsx'), overwrite=T)
}

saveRDS(res, file="OL_type_NB.rds")

}


