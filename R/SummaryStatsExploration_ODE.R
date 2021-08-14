library(gplots)
library(matlab)
# INDIR = "~/Projects/PMO/HighPloidy_DoubleEdgedSword/results/aneuploidyEvolutionSims/"
# BIRTHLANDSCAPEDIR = paste0(INDIR, "../../code/MathModel/GrowOrGo_PloidyEnergyRobustness/data/")
INDIR = "/mnt/ix2/gkimmel/ploidyEvolution/"
BIRTHLANDSCAPEDIR = '/mnt/ix1/Projects/M004_HighPloidy_DoubleEdgedSword_2018/code/MathModel/GrowOrGo_PloidyEnergyRobustness/data/'
inVitroCN = read.table(paste0(BIRTHLANDSCAPEDIR,"CNs brain cancer CLs.txt"), sep="\t", check.names = F, stringsAsFactors = F)
inVitroGR = read.table(paste0(BIRTHLANDSCAPEDIR, "GrowthRate brain cancer CLs.txt"), sep="\t", check.names = F, stringsAsFactors = F)
chr = paste0("Chr",setdiff(colnames(inVitroCN), "ploidy"))
DATEFLAG="2021-05-03_run"
# DATEFLAG="2021-04-30_run"

# f=list.files("~/Dropbox/15_PloidyContinuum/data/testruns", full.names = T)
f=list.files(INDIR, full.names = T, pattern = paste0("output_",DATEFLAG))
f = grep("_params.csv", f, invert = T, fixed = T, value = T)
LASTDAY = sapply(f, function(x) colnames(read.csv(x, header = T, sep="\t", check.names = F, stringsAsFactors = F)))
LASTDAY = unique(sapply(LASTDAY, function(x) x[length(x)]))
## Read simulation parameters
param = sapply(f, function(x) read.csv(gsub(".csv","_params.csv", x), sep=",", check.names = F, stringsAsFactors = F) )

## Prepare stats summary object
stats=matrix(NA, length(f), 4);
colnames(stats) = c("pct_50","pct_70","pct_80","pct_90")

cnv = matrix(NA, length(f), length(chr))
colnames(cnv) = chr

icstats = matrix(NA, length(f),1)
colnames(icstats) = "N_compartments_visited"

before_after = matrix(NA, length(f), 3);
colnames(before_after) = c("end of 1st passage","start of 2nd passage", "day of passaging")

rownames(icstats) <- rownames(before_after) <- rownames(cnv) <- rownames(stats) <- colnames(param) <- sapply(f, function(x) fileparts(x)$name)
allcnv <- allstats <- list()
dois = as.character(paste0(seq(0,as.numeric(LASTDAY), by=5),".0"))
for(doi in dois){
  allstats[[doi]] = stats
  allcnv[[doi]] = cnv
}

## Gather stats
pdf(paste0("~/Downloads/CNVevolution_day45_stats_",DATEFLAG,".pdf"))
par(mfrow=c(3,3))
for(x in f){
  y = read.csv(x, header = T, sep="\t", check.names = F, stringsAsFactors = F)
  # colnames(y)[1:length(chr)]=chr
  # colnames(y)[(length(chr)+1):ncol(y)] = paste0("day_",1:(ncol(y)-length(chr))-1)
  rownames(y) = apply(y[,chr], 1, paste, collapse=",")
  
  for(doi in dois){
    ## sort by comparmtnet size on day of interest (doi)
    y = y[sort(y[,doi], decreasing = T, index.return=T)$ix,]
    pct = sapply(1:nrow(y), function(x) sum(y[1:x,doi]))
    pct = pct/sum(y[,doi])
    
    ## compartments making up 90% of population on day of interest
    for(thr in colnames(allstats[[doi]])){
      thr_ = as.numeric(gsub("pct_","",thr))/100
      y_ = y[pct<=thr_,]
      allstats[[doi]][fileparts(x)$name,thr] = nrow(y_)
    }
    ## Consensus CNV profile
    allcnv[[doi]][fileparts(x)$name, chr] = sapply(chr, function(chr_) sum(y[,chr_]*(y[,doi]/sum(y[,doi]))))
    # heatmap.2(as.matrix(y_[,chr]))
  }
  
  ## survival time of transient compartments
  doi = "350.0"
  y = y[sort(y[,doi], decreasing = T, index.return=T)$ix,]
  y_ = y[,grep("Chr",colnames(y), invert = T)]
  y_ = sweep(y_, 2, STATS = apply(y_, 2, sum), FUN = "/")
  plot(100*as.numeric(y_[1,1:200]), main=paste(y[rownames(y_)[1],1:length(chr)], collapse = ","), pch=20, xlab="day", ylab="% cells")
  
  ## how many compartments does a patient visit on its way to the optimum?
  icstats[fileparts(x)$name,"N_compartments_visited"] = sum(apply(y_[,1:which(colnames(y_)==doi)]>0.01, 1, any))
}
heatmap.2(as.matrix(allcnv[["300.0"]]), margins = c(10,15), labRow = paste(param["misRate",rownames(allcnv[["300.0"]])], rownames(allcnv[["300.0"]])))
par(mfrow=c(2,1))
hist(icstats,20,col='cyan', xlab=paste("# compartments with >1% cell representation at",doi),ylab="# patients", main="")
dev.off()



## Visualize CNVs
## Final state for each patient:
lastCNV = as.matrix(allcnv[[LASTDAY]])
## The unique state(s) all patients evolve to:
targetCNVs = plyr::count(apply(round(lastCNV,1), 1, paste, collapse=","))

targetCNVs = as.numeric(strsplit(as.character(targetCNVs$x[which.max(targetCNVs$freq)]),",")[[1]])
## days to target state for each patient:
dtt = rep(NA,nrow(cnv))
names(dtt) = rownames(cnv)
## patients of interest (those that are evolving)
poi = c()
pdf(paste0("~/Downloads/CNV_evolution_",DATEFLAG,".pdf"))
for(patient in rownames(cnv)){
   x = sapply(allcnv, function(x) x[patient,] )
   x = x[,1:which(colnames(x)==doi)]
   ##Test if patient's cells are evolving in-vitro:
   if(min(apply(x, 1, function(x) length(unique(x))))>1){
     poi = c(poi, patient)
   }
   try(heatmap.2(as.matrix(x),trace = "n",Colv = NULL, dendrogram = "row",main = patient))
   dtt[patient] = names(which(apply(x, 2, function(y) all(round(y,1)==round(targetCNVs,1))))[1])
}
vioplot::vioplot(sort(as.numeric(dtt)), xlab="", ylab="days to target CNV")
dev.off()


## Visualize number of compartments visited for evolving patients
par(mfrow=c(2,2))
for(pct in colnames(stats)){
  x = sapply(allstats, function(s) s[poi,pct,drop=F] )
  boxplot(x, ylab=paste0("# compartments making up ",gsub("pct_","",pct), "%"), col="cyan")
}

pdf(paste0("~/Downloads/PopulationSize_",DATEFLAG,".pdf"))
par(mfrow=c(3,3))
for(x in f){
  y=read.csv(x, header = T, sep="\t", check.names = F, stringsAsFactors = F)
  plot(apply(y[,grep("Chr",colnames(y), invert = T)],2,sum), xlab="day",log="x",ylab="Population size", main=fileparts(x)$name)
}
dev.off()
