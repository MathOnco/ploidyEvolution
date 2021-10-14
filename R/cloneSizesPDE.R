setwd("C:/Users/4473331/Documents/projects/ploidyEvolution-main/Julia")
ff <- list.files()
ff <- ff[grepl("_population",ff)]
ff <- ff[1:300]
#ff <- ff[round(1:300*(length(ff)/100))]

proc <- function(fname,nchrom=5){
  i <- as.numeric(unlist(strsplit(fname,split="_"))[1])
  x <- read.table(fname,sep=",")
  cn <- sapply(1:nrow(x), function(i) as.character(interaction(x[i,1:nchrom])))
  x <- x[,(1+nchrom):ncol(x)]
  x$cn<-cn
  colnames(x)[1:2] <- c("br","n")
  x$t <- i/10
  
  x <- x[x$n>10,]
  
  return(x)
}
library(parallel)
n.cores<-4
clust <- makeCluster(n.cores)
df <- parLapply(clust,ff,proc)
df <- do.call(rbind,df)
library(ggplot2)
p <- ggplot(df,aes(x=t,y=n,group=cn))+
  geom_line(aes(color=br))+
  geom_hline(yintercept=1000,color="green")+
  scale_y_log10("total cells")+
  scale_x_continuous("days")+
  scale_color_viridis_c("birthrate",option="inferno")+
  theme_bw()
p
