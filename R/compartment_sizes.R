setwd("C:/Users/4473331/Documents/projects/ploidyEvolution-main/Julia")
ff <- list.files()
ff <- ff[grepl("_population",ff)]

proc <- function(fname,nchrom=5){
  i <- as.numeric(unlist(strsplit(fname,split="_"))[1])
  x <- read.table(fname,sep=",")
  x <- x[,(nchrom+1):ncol(x)]
  colnames(x)[1:2] <- c("br","n")
  
  ind <- x$n>1000
  
  x <- data.frame(type=c("stoch","pde"),
             n=c(sum(x$n[!ind]),sum(x$n[ind])),
             br=c(sum(x$n[!ind]*x$br[!ind])/sum(x$n[!ind]),sum(x$n[ind]*x$br[ind])/sum(x$n[ind])))
 
  x$t <- i/10
  return(x)
}
library(parallel)
n.cores<-4
clust <- makeCluster(n.cores)
df <- parLapply(clust,ff,proc)
df <- do.call(rbind,df)
library(ggplot2)

p <- ggplot(df,aes(x=br,y=n,color=type))+
  geom_path()+
  geom_hline(yintercept=1000,color="green")+
  scale_y_log10("total cells")+
  scale_x_continuous("avg compartment birthrate")+
  scale_color_viridis_d("")+
  theme_bw()
p

p <- ggplot(df,aes(x=t,y=n,color=type))+
  geom_line()+
  geom_hline(yintercept=1000,color="green")+
  scale_y_log10("total cells")+
  scale_x_continuous("days")+
  scale_color_viridis_d("birthrate",option="inferno")+
  theme_bw()
p



