setwd("~/Repositories/ploidyEvolution/R")
library("TCGAbiolinks")
library("TCGAutils")
library("curatedTCGAData")
library("xlsx"); 
library(matlab); 
devtools::source_url("https://github.com/noemiandor/Utils/blob/master/Pathways/getGenesInvolvedIn.R?raw=TRUE")
devtools::source_url("https://github.com/noemiandor/Utils/blob/master/Pathways/getAllPathways.R?raw=TRUE")
devtools::source_url("https://github.com/noemiandor/Utils/blob/master/DatabaseInterfaces.R?raw=TRUE")
library(Biobase)
library(GSVA)
library(GEOquery)
library(sva)
library("M3C")

## Read phenotypic info from Sam Bakhoum
MN=read.xlsx("../data/Micronuclei SingleCell and BulkRNAseq GSE98183.xlsx",sheetName = "Micronucleus_FromSamBakhoum",stringsAsFactors=F)
MN$Micronuclei=as.numeric(MN$Micronuclei)
MN$Lagging.Chromosome=as.numeric(MN$Lagging.Chromosome)
MN$mouse.number=readr::parse_number(MN$Sample.Name)
MN$metastasis=0
MN$metastasis[grep("metastasis",MN$sample.type)] = 1
rownames(MN)=MN$GEO.ID

################
### GEO data ###
setwd("~/Downloads/")
geo_id="GSE98183"
unlink(geo_id, recursive = TRUE)
gset <- getGEO(geo_id, GSEMatrix =F, getGPL=T)
samples=lapply(gset@gsms, function(x) c(x@header$characteristics_ch1, x@header$title))
A=do.call(rbind,samples[sapply(samples,length)==4])
B=do.call(rbind,samples[sapply(samples,length)==5])
colnames(A) = c("cell line:", "originating cell line: ", "sample type: ","title")
A=as.data.frame(A)
AA=loadGEOdataset(geo_id, loadRaw=T, header = T,sep = ",")
names(AA$data) = sapply(names(AA$data), function(x) fileparts(x)$name)

## Parse IDs from GEO
Z = AA$data$GSE98183_counts.geneSymbols
rownames(Z)=Z[,1]
Z=Z[,-1]
Z = Z[,colnames(Z) %in% A$title]
colnames(Z) = rownames(A)[match(colnames(Z),A$title)]


###################
#### TCGA data ####
cancerList = c("BRCA","HNSC","COAD");#,"ESCA","READ","CESC","SKCM","LUAD","DLBC")
allmat = list()
paths = list()
groups=list()
for (can in cancerList){
  print(paste0("working on ",can))
  X <- curatedTCGAData(diseaseCode = can, assays = "RNASeq2*", dry.run = FALSE, version = "2.0.1")
  # Choose raw data
  i = grep("Norm",names(assays(X)),invert = T)
  tmp=as.matrix(assays(X)[[i]])
  sampleType=sapply(colnames(tmp), function(x) strsplit(x, "-")[[1]][4])
  ctrl=tmp[,grep("11",sampleType)]
  tmr=tmp[,grep("01",sampleType)]
  ## limit samples per can type
  tmr = tmr[,1:min(200,ncol(tmr))]
  Y = cbind(ctrl, tmr)
  groups[[can]] = c(rep("control",ncol(ctrl)), rep(can, ncol(tmr)))
  met=tmp[,grep("06",sampleType),drop=F]
  if(ncol(met)>0){
    Y = cbind(Y,met)
    groups[[can]] = c(groups[[can]], rep("met", ncol(met)))
  }
  allmat[[can]] = Y
}


###############################################################################
### Integrate TCGA & GEO data adjusting for Batch Effects using comBat-seq ####
###############################################################################
allmat$GEO = Z
groups$GEO = rep("BRCA",ncol(Z))
groups$GEO[grep("met",samples[colnames(Z)])] = "met"
## Genes expressed in all datasets
fr=plyr::count(unlist(lapply(allmat,rownames)))
fr=fr[fr$freq==length(allmat),]
### Combine all our gene expression data together into the combined data frame
combined = sapply(allmat, function(x) x[fr$x,])
combined = do.call(cbind,combined)
combined =as.matrix(combined)
### Here we need to assign batch numbers to each batch in our data
batchMat<-sapply(allmat,ncol)
batch<-list()
for (i in 1:length(batchMat)){
  batch[[i]]<-rep(i, batchMat[i])
}
batch<-unlist(batch)
# ## Collapse all TCGA batches into one
# batch[batch!=max(batch)]=1
### Now we apply the ComBat_seq function to adjust for batch effects
adjusted <- ComBat_seq(combined, batch=batch, group=group)


###################################
### VISUALIZE GENE/PATHWAY DATA ###
###################################

### visualize combined gene expression data
pLcombined<-list()
colnames(combined)<-batch
for (bat in unique(batch)){
  pLcombined[[bat]]<-as.numeric(combined[,bat==batch])
  pLcombined[[bat]]=pLcombined[[bat]][pLcombined[[bat]]>0.1]
}
boxplot(pLcombined, log="y")
title("Combined")
### dimensionality reduction
tsne(combined,labels = paste(unlist(groups),batch), perplex = 70)
title("Combined")
umap(combined, labels = paste(unlist(groups),batch))
title("Combined")

### visualize adjusted gene expression data
pLadjusted<-list()
for (bat in unique(batch)){
  pLadjusted[[bat]]<-as.numeric(adjusted[,batch==bat])
  pLadjusted[[bat]] = pLadjusted[[bat]][pLadjusted[[bat]]>0.1]
}
boxplot(pLadjusted, log="y")
title("Adjusted")
### dimensionality reduction
tsne(adjusted,labels = paste(unlist(groups),batch),perplex = 5)
title("adjusted")
umap(adjusted, labels  = paste(unlist(groups),batch) )
title("adjusted")

# ### visualize pathway data
# plotList<-list()
# for (bat in unique(batch)){
#   plotList[[bat]]<-as.numeric(pq[,bat])
# }
# boxplot(plotList)
# tsne(pq,labels = batch)


##############################
## Fit & apply linear model ##
##############################
##Pathway quantification
gs=getAllPathways(include_genes=T, loadPresaved = T);
gs=gs[sapply(gs, length)>=5]
all_pq <- gsva(as.matrix(adjusted), gs, kcdf="Poisson", mx.diff=T, verbose=FALSE, parallel.sz=2, min.sz=10)
pq = all_pq[,paste0("GEO.",colnames(Z))]
colnames(pq) = gsub("GEO.","",colnames(pq))
# pq <- gsva(as.matrix(Z), gs, kcdf="Poisson", mx.diff=T, verbose=FALSE, parallel.sz=2, min.sz=10)

## Compare RNA-seq to micro-nuclei
MN = MN[colnames(pq),]
r=sapply(rownames(pq), function(x) cor(pq[x,rownames(MN)], MN$Lagging.Chromosome,use = "pairwise.complete.obs"))
head(r[sort(abs(r),index.return=T,decreasing = T)$ix])
r[grep("feron",names(r))]
r[grep("IFN",names(r))]

## Fit linear model
poi="Interferon Signaling"
inp=cbind(MN[,c("Lagging.Chromosome","metastasis")],t(pq[poi,rownames(MN),drop=F]))
m = lm(inp$Lagging.Chromosome ~inp$`Interferon Signaling`)
summary(m)
## Visualize
plot(inp$`Interferon Signaling`,inp$Lagging.Chromosome,col=2+inp$metastasis,pch=20,cex=1.5)
o<-predict(m, inp, interval="confidence")
lines(inp$`Interferon Signaling`,o[,"fit"])
legend("bottomleft",paste("metastasis",unique(inp$metastasis)),fill=unique(c(2+inp$metastasis)),bty="n")

## Use model to predict MS rate in TCGA cancers
pq = all_pq[,grep("GEO.",colnames(all_pq), invert = T)]
inp = as.data.frame(matrix(NA,ncol(pq),1));
rownames(inp) = colnames(pq)
colnames(inp) = poi
inp[,poi]=pq[poi,] 
o<-predict(m, inp, interval="confidence")
boxplot(o[,"fit"]~sapply(strsplit(rownames(o),".",fixed = T),"[[",1))
