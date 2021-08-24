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
Z = AA$data$GSE98183_fpkm.geneSymbols
rownames(Z)=Z[,1]
Z=Z[,-1]
Z = Z[,colnames(Z) %in% A$title]
colnames(Z) = rownames(A)[match(colnames(Z),A$title)]

##Pathway quantification
gs=getAllPathways(include_genes=T, loadPresaved = T);     
gs=gs[sapply(gs, length)>=5]
pqT <- gsva(as.matrix(Z), gs, kcdf="Poisson", mx.diff=T, verbose=FALSE, parallel.sz=2, min.sz=10)

## Compare RNA-seq to micro-nuclei
MN = MN[colnames(pqT),]
r=sapply(rownames(pqT), function(x) cor(pqT[x,rownames(MN)], MN$Lagging.Chromosome,use = "pairwise.complete.obs"))
head(r[sort(abs(r),index.return=T,decreasing = T)$ix])
r[grep("feron",names(r))]
r[grep("IFN",names(r))]


###################
#### TCGA data ####
cancerList = c("HNSC","ESCA","COAD","READ","BRCA","CESC","SKCM","LUAD","DLBC")
allmat = list()
paths = list()
for (can in cancerList){
  print(paste0("working on ",can))
  X <- curatedTCGAData(diseaseCode = can, assays = "RNASeq2*", dry.run = FALSE, version = "2.0.1")
  # sampleTables(X)
  # tmp = splitAssays(X, sampleTables(X))
  # Y = tmp[[grep(paste0("01_",can,"_RNASeq2"),
  #               names(tmp))]]@assays$data@listData[[1]]
  # tmp=as.matrix(assays(X)[["HNSC_RNASeq2GeneNorm-20160128"]])
  tmp=as.matrix(assays(X)[[2]])
  sampleType=sapply(colnames(tmp), function(x) strsplit(x, "-")[[1]][4])
  tmp=tmp[,grep("01",sampleType)]
  Y = tmp[,1:min(20,ncol(tmp))]
  allmat[[can]] = Y
}

# tmp=as.matrix(assays(X)[["HNSC_RNASeq2GeneNorm-20160128"]])
# sampleType=sapply(colnames(tmp), function(x) strsplit(x, "-")[[1]][4])
# tmp=tmp[,grep("01",sampleType)]

###############################################################################
### Integrate TCGA & GEO data adjusting for Batch Effects using comBat-seq ####
###############################################################################
combined<-data.frame(ncol=1)
for (name in names(allmat)){
  combined=cbind(combined, allmat[[name]])
}
combined = cbind(combined, Z[rownames(allmat$HNSC),])
combined <- combined[,-c(1)]
combined =as.matrix(combined)
batch <- c(rep(1, 20), rep(2, 20),rep(3, 20),rep(4, 20),rep(5, 20),rep(6, 20),rep(7, 20),rep(8, 20), rep(9, 20), rep(10, 39))
adjusted <- ComBat_seq(combined, batch=batch, group=NULL)
pq <- gsva(as.matrix(adjusted), gs, kcdf="Poisson", mx.diff=T, verbose=FALSE, parallel.sz=2, min.sz=10)
colnames(pq)<-batch
### add together pq objects for plotting in one boxplot
plotList<-list()
for (bat in unique(batch)){
  plotList[[bat]]<-as.numeric(pq[,bat])
  plotList[[bat]]<-plotList[[bat]]+1
}
boxplot(plotList)

## Fit linear model
poi="Interferon Signaling"
inp=cbind(MN[,c("Lagging.Chromosome","metastasis")],t(pqT[poi,rownames(MN),drop=F]))
m = lm(inp$Lagging.Chromosome ~inp$`Interferon Signaling`)
summary(m)
## Visualize
plot(inp$`Interferon Signaling`,inp$Lagging.Chromosome,col=2+inp$metastasis,pch=20,cex=1.5)
o<-predict(m, inp, interval="confidence")
lines(inp$`Interferon Signaling`,o[,"fit"])
legend("bottomleft",paste("metastasis",unique(inp$metastasis)),fill=unique(c(2+inp$metastasis)),bty="n")


