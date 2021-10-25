setwd("C:/Users/4473331/Documents/projects/ploidyEvolution-main/Julia")
library(stringr)
fill_label <- "cells"

ff <- sapply(seq(0,60,5), function(i) paste(str_pad(i,4,"left","0"),"_pdestates.csv",sep=""))


proc_spat <- function(fname,dir){
  i <- as.numeric(unlist(strsplit(fname,split="_"))[1])
  x <- read.table((paste(dir,fname,sep="")),sep=",")
  if (grepl("pde",fname)&i>0){
    t <- unlist(strsplit(fname,split="_"))[1]
    y <- read.table((paste(dir,t,"_stochstates.csv",sep="")),sep=",")
    y <- y[,c(ncol(y)-1,ncol(y))]
    y <- ceiling(y)
    for(j in 1:nrow(y)){
      x[y[j,1],y[j,2]]<-x[y[j,1],y[j,2]]+1
    }
  }
  colnames(x) <- 1:ncol(x)
  x$y <- 1:nrow(x)
  x$y <- x$y%%(ncol(x)-1)
  x$y[x$y==0] <- (ncol(x)-1)
  
  x <- reshape2::melt(x,id.vars="y")
  x <- aggregate(x$value,by=list(x=x$variable,y=x$y),sum)
  colnames(x) <- c("x","y","v")
  x$y <- as.numeric(x$y)
  x$x <- as.numeric(x$x)
  x$y <- max(x$y)-x$y
  x$x <- max(x$x)-x$x
  
  x$i <- i/10
  return(x)
}

x <- do.call(rbind,lapply(ff,proc_spat,dir="C:/Users/4473331/Documents/projects/ploidyEvolution-main/Julia/"))

library(gganimate)
library(ggplot2)

p <- ggplot(x,aes(x=x,y=y,fill=v))+
  geom_raster()+

  theme_classic(base_size=16)+  
  scale_fill_viridis_c(fill_label,trans="log")+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())+
  transition_time(i)+
  labs(title = 'day: {frame_time}')
p

anim_save("cells.gif", animation = last_animation())
