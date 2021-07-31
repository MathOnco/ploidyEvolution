setwd("~/Projects/PMO/HighPloidy_DoubleEdgedSword/code/MathModel/GrowOrGo_PloidyEnergyRobustness/R")
library(TCGAbiolinks)
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
X = AA$data$GSE98183_fpkm.geneSymbols
rownames(X)=X[,1]
X=X[,-1]
X = X[,colnames(X) %in% A$title]
colnames(X) = rownames(A)[match(colnames(X),A$title)]

##Pathway quantification
gs=getAllPathways(include_genes=T, loadPresaved = T);     
gs=gs[sapply(gs, length)>=5]
pq <- gsva(as.matrix(X), gs, kcdf="Poisson", mx.diff=T, verbose=FALSE, parallel.sz=2, min.sz=10)

## Compare RNA-seq to micro-nuclei
MN = MN[colnames(pq),]
r=sapply(rownames(pq), function(x) cor(pq[x,rownames(MN)], MN$Lagging.Chromosome,use = "pairwise.complete.obs"))
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
  X <- curatedTCGAData(diseaseCode = can, assays = "RNASeq2*", dry.run = FALSE)
  # sampleTables(X)
  tmp = splitAssays(X, c("01"))
  Y = tmp[[grep(paste0("01_",can,"_RNASeq2"),
                names(tmp))]]@assays$data@listData[[1]]
  Y = Y[,1:min(100,ncol(Y))]
  allmat[[can]] = Y
}


#########################################
### @TODO: integrate TCGA & GEO data ####
## Follow MetaIntegrator package manual: https://cran.r-project.org/web/packages/MetaIntegrator/vignettes/MetaIntegrator.html
## Use GSVA on integrated dataset

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


