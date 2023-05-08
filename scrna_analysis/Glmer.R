#library(lme4)
library(MASS)
library(Seurat)


TestGlmNB <- function(seur,celltype,formula1=Expression~offset(logUMI)+ko_vs_wt,formula2=Expression~offset(logUMI),cutoff=.05)
{
print("Get Cell Type(s)!")
seur<-subset(seur,cells=names(Idents(seur))[Idents(seur) %in% celltype])

meta=seur@meta.data
meta["logUMI"]=log(meta[,"nCount_RNA"])

meta["nUMI"]=scale(meta[,"nCount_RNA"])

print("Filter Genes and Get Counts")

#mn=rowMeans(seur@raw.data[,names(seur@ident)]>0)
mn=rowMeans(seur@assays$RNA@counts>0)

#dat=seur@raw.data[mn>cutoff,names(seur@ident)]
dat=seur@assays$RNA@counts[mn>cutoff,]
dat=data.frame(as.matrix(dat))

print("Test")
val=0
print(dim(dat))
#test<-lapply(1:dim(dat)[1],function(i){if(i%%10==0){print(i)};x=dat[i,];temp=meta;temp=data.frame(temp);temp["Expression"]=as.numeric(x);res1<-glm(formula1,data=temp,family=poisson);res2<-glm(formula2,data=temp,family=poisson);return(anova(res1,res2))})
test<-lapply(1:dim(dat)[1], function(i){
	x=dat[i,];
	temp=meta;
	temp=data.frame(temp);
	temp["Expression"]=as.numeric(x);
	coef<-tryCatch({res1<-glm.nb(formula1,data=temp);
		coef=summary(res1)$coefficients;
		coef=data.frame(coef);
		coef["Gene"]=rownames(dat)[i];
		coef["Param"]=rownames(coef);
		return(coef)
		}, error = function(err) {
			print(rownames(dat)[i]);
			return(NULL)
			}
		);
	return(coef)})

#pvals<-lapply(test,function(x){print(x);return(x$"Pr(>Chisq)"[2])})
#padj<-p.adjust(as.numeric(pvals),"fdr")
#ret=cbind(rownames(dat)[1:50],pvals,padj)
ret=do.call(rbind,test)
#print(head(ret))
#print(colnames(ret)[4])
colnames(ret)[4]="pval"
ret=ret[grep("^ko_vs_wt",ret[,"Param"]),]
ret["padj"]=p.adjust(ret[,"pval"],"fdr")
ret=ret[,c(5,4,7,1,2,3)]
rownames(ret)=NULL
#print(head(ret))

print("Return!")
return(ret)

}

