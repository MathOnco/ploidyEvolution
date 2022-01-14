library(matlab)
library("xlsx")
library("GSVA")
library("TCGAutils")
library("curatedTCGAData")
library(RMySQL)
library(cloneid)
devtools::source_url("https://github.com/noemiandor/Utils/blob/master/grpstats.R?raw=TRUE")
devtools::source_url("https://github.com/MathOnco/ploidyEvolution/blob/main/R/Utils.R?raw=TRUE")
setwd("../Data")

###########################################
#### Cell lines: bulk expression & CNVs ###
###########################################
cellpass = read.csv("/mnt/ix1/Resources/data/Databases/CellPassports/model_list_20200825.csv")
cn=read.table("/mnt/ix1/Resources/data/Databases/CCLE/CCLE_copynumber_2017-06-22.nocnv.seg.txt",sep="\t", header=T)
ex=read.table("/mnt/ix1/Resources/data/Databases/CCLE/CCLE.rpkm.2016-06-17f.gct", skip = 2, header = T, check.names = F, stringsAsFactors = F, sep="\t")
appCL=read.table("/mnt/ix1/Resources/data/Databases/Cmap_LINCS/L1000/Cell_app_export.txt",sep="\t",check.names = F, stringsAsFactors = F, header = T, quote = "", comment.char = "")
appCL = appCL[appCL$`CCLE name` %in% cellpass$CCLE_ID,]
appCL$ploidy = cellpass$ploidy[match(appCL$`CCLE name`,cellpass$CCLE_ID)]

## Format copy numbers:
cn$CN_Estimate=1*2^cn$Segment_Mean
cn$seglength = 1 + cn$End - cn$Start
colnames(cn)[2:4]=c("chr","startpos","endpos")
cn$CCLE_name=as.character(cn$CCLE_name)
cn=cn[cn$chr %in% 1:22,]; 
cn$chr=as.numeric(cn$chr)

## Add ploidy information to copy numbers:
cn$ploidy = NA
for(x in intersect(appCL$`CCLE name`,cn$CCLE_name)){
  ii = which(cn$CCLE_name==x)
  cn$ploidy[ii] = appCL$ploidy[appCL$`CCLE name`==x]
  cn$CN_Estimate[ii] = round(cn$ploidy[ii]) * cn$CN_Estimate[ii]
  cn$ploidy[ii] = round(cn$ploidy[ii])
}
cn=lapply(unique(cn$CCLE_name), function(x) cn[cn$CCLE_name==x,])
names(cn)= sapply(cn, function(x) x[1,1])

## Test calculation of mis-segregation imprints:
test  = calcMSimprints(cn$T98G_CENTRAL_NERVOUS_SYSTEM, p = 0, cnColumn = "CN_Estimate", precision = 1)
print(test$landscape[c("7","10","9")])




##################################################
#### Cell lines: single-cell expression & CNVs ###
##################################################
cellLine = "SNU-16"
## scDNA-seq derived clone profiles
clones_d = getSubProfiles(cellLine, whichP = "GenomePerspective")
## scRNA-seq derived clone profiles
clones_r = getSubProfiles(cellLine, whichP = "TranscriptomePerspective")

## Load single cell data and count karyotypes
load('listOfSeurats_G0G1_d.RData')
X=sapply(listOfSeurats, function(x) x$X)
X=sapply(X, function(x) x[grep(":",rownames(x),fixed=T),])
X=sapply(X, function(x) cbind(cloneid::parseLOCUS(rownames(x)),x))
K=calculateKaryotypes(X)
## How many karyotypes are present across multiple CLs
Kunique=sapply(K, function(a) unique( apply(a,2,paste,collapse=",")))
Kunique=do.call(c,Kunique)
fr=plyr::count(Kunique)
fr=fr[order(fr$freq,decreasing = T),]
print(paste("Karyotypes are present across >1 CLs:",sum(fr$freq>1)))
## Make sure to use only chromosomes represented by all cell lines
fr=plyr::count(unlist(lapply(K,rownames)))
fr=fr[fr$freq==length(X),]
K=lapply(K, function(x) x[fr$x,])
## Merge across CLs
K=do.call(cbind,K)
K=apply(K,2,paste,collapse=",")
fr=plyr::count(K)
fr=fr[order(fr$freq,decreasing = T),]
ggplot(data.frame(fr$freq), aes(fr.freq)) + geom_histogram(bins = 30) +
  scale_x_log10() +scale_y_log10() + ggtitle(paste(length(K),"cells:",nrow(fr),"unique karyotypes"))



##############################
#### primary tumors: CNVs ####
##############################
panDir = "/mnt/ix1/Resources/data/PANCAN_WholeGenomes/"
tables1 = read.xlsx(paste0(panDir,"supplementary Tables/Supplementary Table 1.xlsx"), sheetIndex=1, stringsAsFactors=F,startRow=3)
clinical = read.xlsx(paste0(panDir,"Clinical/pcawg_donor_clinical_August2016_v9.xlsx"), sheetIndex=1, stringsAsFactors=F)
clinical = clinical[!is.na(clinical$project_code),]
ploidy = read.table(paste0(panDir,"CNV/consensus.20170218.purity.ploidy.txt"), sep="\t", check.names=F, stringsAsFactors=F, header=T)
ploidy$icgc_donor_id = tables1$icgc_donor_id[ match(ploidy$samplename,tables1$tumour_specimen_aliquot_id) ]
PANCAN = list()
for(i in 1:nrow(ploidy)){
  f = list.files(paste0(panDir,"CNV/somatic.cna.tcga"), full.names = T, pattern = ploidy$samplename[i])
  PANCAN[[ploidy$samplename[i]]] = read.table(f, sep="\t", check.names=F, stringsAsFactors=F, header=T)
}

####################################
#### primary tumors: Expression ####
####################################
cancers = unique(sapply(strsplit(clinical$project_code,"-"),"[[",1))
TCGA = list()
for(can in intersect(cancers,c("GBM","STAD"))){
  TCGA[[can]]  <- curatedTCGAData(diseaseCode = can, assays = c("CNASeq", "RNASeq2*", "mRNAArray"), dry.run = FALSE)
}






############################################
#### GBM IHC and H&E whole slide images ####
############################################
d=list.dirs("/mnt/ix1/Projects/M004_HighPloidy_DoubleEdgedSword_2018/data/Konstantin_GBM", full.names = T,recursive = F)
d=d[grep("Slides",d)]
cells_vessels=list()
for (d_ in d){
  ## Read QuPath outputs
  f=list.files(d_,pattern=".txt",full.names = T)
  v = grep("vessels",f,value=T)
  f = grep("H_E",f,value=T)
  if(isempty(v)){
    next;
  }
  cells_vessels[[fileparts(f)$name]] = read.table(v,header = T,check.names = F,stringsAsFactors = F,sep="\t")
}
