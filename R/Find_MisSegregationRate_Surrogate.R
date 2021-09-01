setwd("~/Repositories/ploidyEvolution/R")
library("TCGAbiolinks")
library("TCGAutils")
library("curatedTCGAData")
library("xlsx"); 
library(matlab); 
library(Biobase)
library(GSVA)
library(GEOquery)
library(sva)
library(wesanderson)
library("M3C")
source("Utils.R")
devtools::source_url("https://github.com/noemiandor/Utils/blob/master/grpstats.R?raw=TRUE")
devtools::source_url("https://github.com/noemiandor/Utils/blob/master/DatabaseInterfaces.R?raw=TRUE")
devtools::source_url("https://github.com/noemiandor/Utils/blob/master/Pathways/getGenesInvolvedIn.R?raw=TRUE")
devtools::source_url("https://github.com/noemiandor/Utils/blob/master/Pathways/getAllPathways.R?raw=TRUE")
devtools::source_url("https://github.com/noemiandor/Utils/blob/master/intersect_MatlabV.R?raw=TRUE")
cancerList = c("BRCA","HNSC","COAD","ESCA","READ","CESC","SKCM","LUSC","DLBC")

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
od=getwd()
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
setwd(od)

## Parse IDs from GEO
Z = AA$data$GSE98183_counts.geneSymbols
rownames(Z)=Z[,1]
Z=Z[,-1]
Z = Z[,colnames(Z) %in% A$title]
colnames(Z) = rownames(A)[match(colnames(Z),A$title)]


###################
#### TCGA data ####
allkaryo <- allrna <- list()
paths = list()
groups = list()
for (can in cancerList) {
  print(paste0("working on ", can))
  # CNV data
  allkaryo[[can]] = as.data.frame(getCNVsFromTCGA(can))
  # Choose raw RNAseq data
  X <-
    curatedTCGAData(
      diseaseCode = can,
      assays = "RNASeq2*",
      dry.run = FALSE,
      version = "2.0.1"
    )
  # Choose raw data
  i = grep("Norm", names(assays(X)), invert = T)
  tmp = as.matrix(assays(X)[[i]])
  sampleType = sapply(colnames(tmp), function(x)
    strsplit(x, "-")[[1]][4])
  tmr = tmp[, grep("01", sampleType)]
  ## limit non BRCA samples per can type
  if (can != "BRCA") {
    tmr = tmr[, 1:min(100, ncol(tmr))]
  }
  # gather control for BRCA else not
  if (can == "BRCA") {
    ctrl = tmp[, grep("11", sampleType), drop = F]
    Y = cbind(ctrl, tmr)
    groups[[can]] = c(rep("control", ncol(ctrl)), rep(can, ncol(tmr)))
  } else {
    Y = tmr
    groups[[can]] = c(rep(can, ncol(tmr)))
  }
  met = tmp[, grep("06", sampleType), drop = F]
  if (ncol(met) > 0 && can == "BRCA") {
    Y = cbind(Y, met)
    groups[[can]] = c(groups[[can]], rep("met", ncol(met)))
  }
  allrna[[can]] = Y
}


###############################################################################
### Integrate TCGA & GEO data adjusting for Batch Effects using comBat-seq ####
###############################################################################
allrna$GEO = Z
groups$GEO = rep("BRCA",ncol(Z))
groups$GEO[grep("met",samples[colnames(Z)])] = "met"
## Genes expressed in all datasets
fr=plyr::count(unlist(lapply(allrna,rownames)))
fr=fr[fr$freq==length(allrna),]
### Combine all our gene expression data together into the combined data frame
combined = sapply(allrna, function(x) x[fr$x,])
combined = do.call(cbind,combined)
combined =as.matrix(combined)
### Here we need to assign batch numbers to each batch in our data
batchMat<-sapply(allrna,ncol)
batch<-list()
for (i in 1:length(batchMat)){
  batch[[i]]<-rep(i, batchMat[i])
}
batch<-unlist(batch)

### read in table w/gene lengths to normalize raw counts
geneInfo<-annotateFromBioMart(as.data.frame(1:nrow(combined),row.names=rownames(combined)), join_id='hgnc_symbol')
geneLengths<-as.data.frame(geneInfo$endpos-geneInfo$startpos, row.names = geneInfo$hgnc_symbol)
colnames(geneLengths)<-"Length"
## Transform counts to TPM
combined = combined[rownames(combined) %in% rownames(geneLengths),]
combined = 100*sweep(combined, MARGIN=1, geneLengths[rownames(combined), "Length"],FUN = "/")

### Now we apply the ComBat_seq function to adjust for batch effects
adjusted <- ComBat_seq(combined, batch=batch!=max(batch), group=unlist(groups))




###################################
### VISUALIZE GENE/PATHWAY DATA ###
###################################
lab=paste(unlist(groups),batch); 
fr=plyr::count(lab);
rownames(fr) = fr$x
for(x in fr$x){ 
  lab= gsub(paste0(x,"$"),paste0(x," (n=",fr[x,"freq"],")"),lab)
}
pdf("~/Downloads/MS_surrogate_batchCorrection_TPM.pdf")
### visualize combined gene expression data
pLcombined<-list()
for (bat in unique(batch)){
  pLcombined[[bat]]<-as.numeric(combined[,bat==batch])
  pLcombined[[bat]]=pLcombined[[bat]][pLcombined[[bat]]>0.1]
}
umap(combined, labels  = lab,legendtitle = "combined",dotsize=2)

### visualize adjusted gene expression data
pLadjusted<-list()
for (bat in unique(batch)){
  pLadjusted[[bat]]<-as.numeric(adjusted[,batch==bat])
  pLadjusted[[bat]] = pLadjusted[[bat]][pLadjusted[[bat]]>0.1]
}
umap(adjusted, labels  = lab,legendtitle = "adjusted",dotsize=2)
dev.off()


##############################
## Fit & apply linear model ##
##############################
##Pathway quantification
gs=getAllPathways(include_genes=T, loadPresaved = T);
gs=gs[sapply(gs, length)>=5]
all_pq <- gsva(as.matrix(adjusted), gs, kcdf="Poisson", mx.diff=T, verbose=FALSE, parallel.sz=4, min.sz=10)
pq = all_pq[,paste0("GEO.",colnames(Z))]
colnames(pq) = gsub("GEO.","",colnames(pq))

## Compare RNA-seq to micro-nuclei
MN = MN[colnames(pq),]
r=sapply(rownames(pq), function(x) cor(pq[x,rownames(MN)], MN$Lagging.Chromosome,use = "pairwise.complete.obs"))
head(r[sort(abs(r),index.return=T,decreasing = T)$ix])
r[grep("feron",names(r))]
r[grep("IFN",names(r))]

## Fit linear model
poi="Interferon Signaling"
inp=cbind(MN[,c("Lagging.Chromosome","Micronuclei","metastasis")],t(pq[poi,rownames(MN),drop=F]), "cell.line"=MN$cell.line)
inp = inp[sort(inp[,poi],index.return=T)$ix,]
inp$Lagging.Chromosome = log(inp$Lagging.Chromosome)
m = lm(inp$Lagging.Chromosome ~inp[[poi]])
summary(m)
## Visualize
## Set up colors
dd <- unique(inp$cell.line)
dd.col <- wes_palette(n=length(dd), "Darjeeling2", type = "discrete")
names(dd.col)  <- dd
pdf("~/Downloads/MS_predictions_9Cancers_TPM_.pdf")
plot(inp[[poi]],exp(inp$Lagging.Chromosome),col=dd.col,xlab=poi,pch=ifelse(inp$metastasis>0, 8,20),cex=2,main=paste("R^2 =",summary(m)$adj.r.squared))
o<-predict(m, inp, interval="confidence")
lines(inp[[poi]],exp(o[,"fit"]))
# legend("topright",paste("metastasis",unique(inp$metastasis)),fill=unique(c(3+inp$metastasis)),bty="n",cex=1.5)
legend("topright",names(dd.col),fill=dd.col,bty="n",cex=1.5)


## Use model to predict MS rate in TCGA cancers
pq = all_pq#[,grep("GEO.",colnames(all_pq), invert = T)]
inp = as.data.frame(matrix(NA,ncol(pq),1));
rownames(inp) = colnames(pq)
colnames(inp) = poi
inp[,poi]=pq[poi,] 
inp$Lagging.Chromosome = exp(predict(m, inp, interval="confidence")[,"fit"])
inp$batch=batch
inp$group=unlist(groups)
boxplot(inp$Lagging.Chromosome~sapply(strsplit(rownames(inp),".",fixed = T),"[[",1),las=2,xlab="")
boxplot(inp$Lagging.Chromosome~paste(inp$group,inp$batch),las=2,xlab="",ylab="% Lagging Chromosome")
boxplot(inp$Lagging.Chromosome~inp$group,las=2,xlab="",ylab="% Lagging Chromosome")
dev.off()

## Save output
inp = inp[inp$batch!=max(inp$batch) & inp$group!="control",]
write.table(inp, file="~/Downloads/Predicted_LaggingChrPct_TCGA.txt",sep="\t",quote = F)



# ##########################
# ## Karyotype saturation ##
# coi=names(allkaryo)
# par(mfrow=c(3,3))
# karyotypes<-list()
# for (can in coi){
#   print(paste0("working on ",can))
#   # CNV data
#   N=nrow(allkaryo[[can]])
#   karyo_u = dplyr::distinct(round(allkaryo[[can]]))
#   karyotypes[[can]] = c(nrow(karyo_u),nrow(allkaryo[[can]]))
#   ## Sample to estimate how much we miss:
#   saturation = as.data.frame(matrix(NA, N,2))
#   colnames(saturation) = c("all","unique")
#   for(x in 2:N){
#     ij = sample(1:N,x)
#     tmp = dplyr::distinct(round(allkaryo[[can]][ij,]))
#     saturation[x,] = c(x,nrow(tmp))
#   }
#   plot(saturation$all,saturation$unique,pch=20,main=can,ylim=quantile(saturation$all,c(0,1),na.rm=T))
# }
# karyotypes = sapply(karyotypes,function(x) x)
# rownames(karyotypes) = c("unique","all")
# pct=1-karyotypes["unique",]/karyotypes["all",]
# print(paste("Saturation:", names(pct),pct))
# nkaryo = karyotypes["unique",]*1/pct
# print(paste("Expected Karyotypes:", names(pct),nkaryo))
# ## Candidates for MIE?
# quantile(karyotypes["unique",])
# coi=coi[sort(nkaryo[coi],index.return=T)$ix]
# karyotypes[,coi]
# hist(nkaryo,40,col="blue",border="white");
