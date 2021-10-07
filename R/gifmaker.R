setwd("C:/Users/4473331/Documents/projects/ploidyEvolution-main/Julia/test_output/00_full_tests/misrate_0_1/")
fill_label <- "cells"
ff <- list.files()
ff <- ff[grepl("_states",ff)]
ff <- ff[round(1:100*(length(ff)/100))] #downsample

proc <- function(fname){
  i <- as.numeric(unlist(strsplit(fname,split="_"))[1])
  x <- read.table(fname,sep=",")
  colnames(x) <- 1:ncol(x)
  x$y <- 1:nrow(x)
  
  x <- reshape2::melt(x,id.vars="y")
  colnames(x) <- c("y","x","v")
  x$t <- i/10
  return(x)
}

x <- do.call(rbind,lapply(ff,proc))

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
  transition_time(t)+
  labs(title = 'day: {frame_time}')
p

anim_save("cells.gif", animation = last_animation())