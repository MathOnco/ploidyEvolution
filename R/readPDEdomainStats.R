devtools::source_url("https://github.com/noemiandor/Utils/blob/master/grpstats.R?raw=TRUE")
library(matlab)
Plot_ConcaveHull <- function(xx, yy,lcolor="black", alpha=0.4, add=T, level=0.5/length(xx)) {
  library(MASS) 
  ##Remove outliers
  hQ=0.975; lQ=0.025
  iK1=which(xx<=quantile(xx,hQ) & xx>=quantile(xx,lQ))
  iK2=which(yy<=quantile(yy,hQ) & yy>=quantile(yy,lQ))
  iK=intersect(iK1,iK2)
  xx=xx[iK]; yy=yy[iK]; 
  
  ##Contour
  dens2 <- kde2d(xx, yy, lims=c(min(xx)-sd(xx), max(xx)+sd(xx),   
                                min(yy)-sd(yy), max(yy)+sd(yy)),n=55  )
  contour(dens2$x, dens2$y, dens2$z,col=lcolor, add=add, alpha=alpha); #,drawlabels=F,lwd=2 
  return(dens2)
}

OUTD="~/Projects/PMO/HighPloidy_DoubleEdgedSword/data/Konstantin_GBM/VesselCoordinates"
## Statistic to evaluate vessel thickness
for(STAT in c("Nucleus: Area","Nucleus: Perimeter","Nucleus: Eosin OD mean","Nucleus: Hematoxylin OD mean")){
  
  ## Stats objects
  d=list.dirs("~/Projects/PMO/HighPloidy_DoubleEdgedSword/data/Konstantin_GBM/", full.names = T,recursive = F)
  d=d[grep("Slides",d)]
  stats=matrix(0,length(d),3+5)
  rownames(stats) = sapply(d, function(x) fileparts(x)$name)
  colnames(stats) =  c("width","height","nCells",paste("vessels",c("0%","25%","50%","75%","100%")))
  for (d_ in d){
    ## Read QuPath outputs
    f=list.files(d_,pattern=".txt",full.names = T)
    vtest = grep("vesselMember",f,value=T)
    v = grep("vessels",f,value=T)
    f = grep("H_E",f,value=T)
    nuclei = read.table(f,header = T,check.names = F,stringsAsFactors = F,sep="\t")
    
    size_xy = apply(nuclei[,c("Centroid X µm","Centroid Y µm")], 2,quantile, c(0,1))
    stats[fileparts(d_)$name,c("width","height")] = size_xy[2,]-size_xy[1,]
    stats[fileparts(d_)$name,"nCells"]=nrow(nuclei)
    
    if(isempty(v)){
      next;
    }
    cells_vessels = read.table(v,header = T,check.names = F,stringsAsFactors = F,sep="\t")
    plot(-cells_vessels$`Centroid X µm`,cells_vessels$`Centroid Y µm`,pch=20,main=fileparts(v)$name)
    vessels=cells_vessels[cells_vessels$Class=="Positive",]
    points(-vessels$`Centroid X µm`,vessels$`Centroid Y µm`,pch=20,col="red")
    stats[fileparts(d_)$name,grep("vessel",colnames(stats))] = quantile(vessels$`Nucleus: Perimeter`)
    
    ## @TODO: use to test clustering below
    if(!isempty(vtest)){
      vesselTest = read.table(vtest,header = T,check.names = F,stringsAsFactors = F,sep="\t")
      vesselTest=vesselTest[vesselTest$`Num Positive`>0,]
      points(-vesselTest$`Centroid X µm`,vesselTest$`Centroid Y µm`,pch=20,col="cyan",cex=log(vesselTest$`Num Positive`)+0.5)
    }
    
    
    ## Output file format: every row is a coordinate within your grid occupied by blood vesseld
    ## (i.e. all coordinates not listed in the file are to be considered negative for blood vessel)
    o=dbscan::dbscan(vessels[,c("Centroid X µm","Centroid Y µm")],eps = 20,minPts = 1)
    vessels$VesselID = o$cluster
    fr=plyr::count(o$cluster)
    fr=fr[fr$freq>=6,]
    vessels=vessels[vessels$VesselID %in% fr$x,]
    # points(-vessels$`Centroid X µm`,vessels$`Centroid Y µm`,pch=20,col=o$cluster+1)
    coord=list()
    for(i in unique(vessels$VesselID)){
      ii=which(vessels$VesselID==i)
      onevessel=Plot_ConcaveHull(-vessels$`Centroid X µm`[ii],vessels$`Centroid Y µm`[ii],lcolor = i+1)
      
      xy=which(onevessel$z>quantile(onevessel$z,0.7), arr.ind=TRUE)
      
      coord_=as.data.frame(t(apply(xy, 1, function(kl) c(onevessel$x[kl[1]],onevessel$y[kl[2]]))))
      colnames(coord_)=c("Centroid X µm","Centroid Y µm")
      coord_$VesselID=i
      # plot(coord_$`Centroid X µm`,coord_$`Centroid Y µm`,pch=20)
      coord[[as.character(i)]]=coord_
      points(-vessels$`Centroid X µm`[ii],vessels$`Centroid Y µm`[ii],pch=20,col=i+1)
    }
    allcoord=do.call(rbind,coord)
    points(allcoord$`Centroid X µm`,allcoord$`Centroid Y µm`,pch=20,cex=0.4,col="brown")
    
    ## Remove noise
    vessels=vessels[vessels$VesselID!=0,]
    ## Evaluate various vessel thickness metrics
    vesselstats=grpstats(vessels[,c("VesselID",grep("Nucleus",colnames(vessels),value=T))],vessels$VesselID,statscols = c("mean","numel"))
    vesselstats=as.data.frame(vesselstats$mean)
    vesselstats[,STAT]=vesselstats[,STAT]/max(vesselstats[,STAT])
    vesselstats[,STAT]=round(1+100*(vesselstats[,STAT]-min(vesselstats[,STAT])))
    vesselstats=vesselstats[sort(vesselstats[,STAT],index.return=T)$ix,]
    col=rainbow(1.2*max(vesselstats[,STAT]))
    col=col[vesselstats[,STAT]]
    names(col)=as.character(vesselstats$VesselID)
    pdf(paste0(OUTD,filesep,fileparts(d_)$name,"_xy_",STAT,".pdf"),width = 10,height = 5)
    par(mfrow=c(1,2))
    plot(-cells_vessels$`Centroid X µm`,cells_vessels$`Centroid Y µm`,pch=20,main=fileparts(v)$name)
    points(allcoord$`Centroid X µm`,allcoord$`Centroid Y µm`,pch=20,cex=0.4,col=col[as.character(allcoord$VesselID)])
    plot(1,main="vessel thickness",xaxt = "n",yaxt = "n",axes=F,xlab="",ylab="",col="white"); 
    legend("topleft",as.character(unique(vesselstats[,STAT])),fill=unique(col),cex=0.75)
    dev.off()
    
    ## Save output
    allcoord$`Centroid X µm`=-allcoord$`Centroid X µm`
    write.table(allcoord,file=paste0(OUTD,filesep,fileparts(d_)$name,"_xy.txt"),sep="\t",quote = F,row.names = F)
    
  }
  
  print(apply(stats,2,quantile,c(0,1)))
}


