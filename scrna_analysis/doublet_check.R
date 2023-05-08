library(Seurat)
library(Matrix)

args = commandArgs(trailingOnly=T)
seur.file = args[1]
save.dir = args[2]
if (endsWith(save.dir, '/')) gsub("/$", "", save.dir)
split.by = args[3]

code.dir = "~/kwanho/src"

if (!dir.exists(save.dir)) {dir.create(save.dir)}

seur <- readRDS(seur.file)

seur.list <- SplitObject(seur, split.by=split.by)

doub.score <- c()
doub.pred <- c()

for (i in 1:length(seur.list))
{
  nam = names(seur.list)[i]
  print(nam)
  
  obj = seur.list[[nam]]
  print(obj)

  mat.file = paste(save.dir,"/MM.",nam,".for.scrublet.txt",sep="")
  if (!file.exists(mat.file)) {
    writeMM(obj@assays$RNA@counts, mat.file)
  }
  
  # Run Scrublet
  system(paste("/stanley/levin_dr/kwanho/anaconda3/envs/kb/bin/python ",code.dir,"/run_scrublet.py ", mat.file, " ", save.dir,"/scrublet_res_", nam, sep=""))

  # Collect results
  scrub.s <- scan(paste(save.dir,"/scrublet_res_", nam, "_score.txt", sep=""))
  scrub.p <- scan(paste(save.dir,"/scrublet_res_", nam, "_pred.txt", sep=""))

  # Add cell names
  cn = names(doub.score)
  doub.score <- c(doub.score, scrub.s)
  doub.pred <- c(doub.pred, scrub.p)
  names(doub.score) <- c(cn, colnames(obj))
  names(doub.pred) <- c(cn, colnames(obj))
}

seur$doublet_score <- doub.score[colnames(seur)]
seur$doublet_pred <- doub.pred[colnames(seur)]
saveRDS(seur@meta.data[,c('doublet_score','doublet_pred')], paste0(save.dir, "/metadata_scrublet_res_for_seurat.rds"))



