library(gplots)
library(RColorBrewer)
library(chebpol)
library(matlab)
BIRTHLANDSCAPEDIR = '/mnt/ix1/Projects/M004_HighPloidy_DoubleEdgedSword_2018/code/MathModel/GrowOrGo_PloidyEnergyRobustness/data/'
setwd(paste0(BIRTHLANDSCAPEDIR,filesep,"../R"))
source("Utils.R")

chr = c("5","9","13","14","20")
inVitroCN = read.table(paste0(BIRTHLANDSCAPEDIR,"CNs brain cancer CLs.txt"), sep="\t", check.names = F, stringsAsFactors = F)
inVitroGR = read.table(paste0(BIRTHLANDSCAPEDIR, "GrowthRate brain cancer CLs.txt"), sep="\t", check.names = F, stringsAsFactors = F)
r = inVitroGR$V2
names(r) = inVitroGR$V1

o = interpolateBirthRateLandscape(cnv = inVitroCN[names(r),chr], r = r, ndim = 4,  usePCA = F, save="~/Downloads/BirthRateLandscape_Test.pdf"); 

