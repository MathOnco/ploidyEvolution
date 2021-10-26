devtools::source_url("https://github.com/noemiandor/Utils/blob/master/grpstats.R?raw=TRUE")
library(matlab)
Plot_ConcaveHull <- function(xx, yy,lcolor="black", alpha=0.4, add=T, level=55) {
  library(MASS) 
  ##Remove outliers
  hQ=0.975; lQ=0.025
  iK1=which(xx<=quantile(xx,hQ) & xx>=quantile(xx,lQ))
  iK2=which(yy<=quantile(yy,hQ) & yy>=quantile(yy,lQ))
  iK=intersect(iK1,iK2)
  xx=xx[iK]; yy=yy[iK]; 
  
  ##Contour
  dens2 <- kde2d(xx, yy, lims=c(min(xx)-sd(xx), max(xx)+sd(xx),   
                                min(yy)-sd(yy), max(yy)+sd(yy)),n=level  )
  contour(dens2$x, dens2$y, dens2$z,col=lcolor, add=add, alpha=alpha); #,drawlabels=F,lwd=2 
  return(dens2)
}
minmax_x=c(200,840)
minmax_y=c(200,850)
OUTD="~/Projects/PMO/HighPloidy_DoubleEdgedSword/data/Konstantin_GBM/VesselCoordinates"
d=list.dirs("~/Projects/PMO/HighPloidy_DoubleEdgedSword/data/Konstantin_GBM/", full.names = T,recursive = F)
d=d[grep("Slides",d)]
## Statistic to evaluate vessel thickness
for(STAT in c("Nucleus: 6069-20-H 2 mean")){ #},"Nucleus: Area","Nucleus: Perimeter","Nucleus: Eosin OD mean","Nucleus: Hematoxylin OD mean")){
  
  ## Stats objects
  stats=matrix(0,length(d),3+5)
  rownames(stats) = sapply(d, function(x) fileparts(x)$name)
  colnames(stats) =  c("width","height","nCells",paste("vessels",c("0%","25%","50%","75%","100%")))
  for (d_ in d){
    ## Read QuPath outputs
    f=list.files(d_,pattern=".txt",full.names = T)
    vtest = grep("vesselMember",f,value=T)
    v = grep("vessels",f,value=T)
    f = grep("H_E",f,value=T)
    cells = read.table(f,header = T,check.names = F,stringsAsFactors = F,sep="\t")
    
    size_xy = apply(cells[,c("Centroid X µm","Centroid Y µm")], 2,quantile, c(0,1))
    stats[fileparts(d_)$name,c("width","height")] = size_xy[2,]-size_xy[1,]
    stats[fileparts(d_)$name,"nCells"]=nrow(cells)
    
    if(isempty(v)){
      next;
    }
    vessels = read.table(v,header = T,check.names = F,stringsAsFactors = F,sep="\t")
    vessels$Class="Positive"
    cells$Class="Negative"
    cells$`Nucleus: 6069-20-H 2 mean`=NA
    jj=intersect(colnames(cells),colnames(vessels))
    cells_vessels=rbind(vessels[,jj],cells[,jj])
    cells_vessels[,c("Centroid X µm","Centroid Y µm")]=apply(cells_vessels[,c("Centroid X µm","Centroid Y µm")],2,function(x) x-min(x))
    # rect(-minmax_x[2],minmax_y[1],-minmax_x[1],minmax_y[2],col="purple",density = 0,lwd = 5)
    # cells_vessels = cells_vessels[cells_vessels$`Centroid X µm`>minmax_x[1] & 
    #                                 cells_vessels$`Centroid X µm`<minmax_x[2] & 
    #                                 cells_vessels$`Centroid Y µm`>minmax_y[1] & 
    #                                 cells_vessels$`Centroid Y µm`<minmax_y[2],]
    vessels=cells_vessels[cells_vessels$Class=="Positive",]
    cells=cells_vessels[cells_vessels$Class!="Positive",]
    plot(cells_vessels$`Centroid X µm`,-cells_vessels$`Centroid Y µm`,pch=20,main=fileparts(v)$name,cex=0.2)
    points(vessels$`Centroid X µm`,-vessels$`Centroid Y µm`,pch=20,col="red")
    stats[fileparts(d_)$name,grep("vessel",colnames(stats))] = quantile(vessels$`Nucleus: Perimeter`)
    
    ## @TODO: use to test clustering below
    if(!isempty(vtest)){
      vesselTest = read.table(vtest,header = T,check.names = F,stringsAsFactors = F,sep="\t")
      vesselTest=vesselTest[vesselTest$`Num Positive`>0,]
      points(-vesselTest$`Centroid X µm`,vesselTest$`Centroid Y µm`,pch=20,col="cyan",cex=log(vesselTest$`Num Positive`)+0.5)
    }
    
    
    ## Output file format: every row is a coordinate within your grid occupied by blood vesseld
    ## (i.e. all coordinates not listed in the file are to be considered negative for blood vessel)
    # o=dbscan::dbscan(vessels[,c("Centroid X µm","Centroid Y µm")],eps = 20,minPts = 1)
    o=dbscan::dbscan(vessels[,c("Centroid X µm","Centroid Y µm")],eps = 40,minPts = 1)
    vessels$VesselID = o$cluster
    fr=plyr::count(o$cluster)
    # fr=fr[fr$freq>=5,]
    fr=fr[fr$freq>=10,]
    vessels=vessels[vessels$VesselID %in% fr$x,]
    # points(-vessels$`Centroid X µm`,vessels$`Centroid Y µm`,pch=20,col=o$cluster+1)
    coord=list()
    for(i in unique(vessels$VesselID)){
      ii=which(vessels$VesselID==i)
      onevessel=try(Plot_ConcaveHull(-vessels$`Centroid X µm`[ii],vessels$`Centroid Y µm`[ii],lcolor = i+1))
      if(class(onevessel)=="try-error"){
        next;
      }
      q=quantile(onevessel$z,(1:20)/20)
      xy=which(onevessel$z>q["85%"], arr.ind=TRUE)
      edge_xy=which(onevessel$z<=q["85%"] & onevessel$z>q["75%"], arr.ind=TRUE)
        
      coord_=as.data.frame(t(apply(xy, 1, function(kl) c(onevessel$x[kl[1]],onevessel$y[kl[2]]))))
      edge_coord=as.data.frame(t(apply(edge_xy, 1, function(kl) c(onevessel$x[kl[1]],onevessel$y[kl[2]]))))
      colnames(coord_) <- colnames(edge_coord) <-c("Centroid X µm","Centroid Y µm")
      coord_$Boundary = 0
      edge_coord$Boundary = 1
      coord_ = rbind(coord_,edge_coord)
      coord_$VesselID=i
      # plot(coord_$`Centroid X µm`,coord_$`Centroid Y µm`,pch=20)
      coord[[as.character(i)]]=coord_
      points(-vessels$`Centroid X µm`[ii],vessels$`Centroid Y µm`[ii],pch=20,col=i+1)
    }
    allcoord=do.call(rbind,coord)
    points(allcoord$`Centroid X µm`,allcoord$`Centroid Y µm`,pch=20,cex=0.4,col=3+allcoord$Boundary)
    
    ## Remove noise
    vessels=vessels[vessels$VesselID!=0,]
    ## Evaluate various vessel thickness metrics
    vesselstats=grpstats(vessels[,c("VesselID",grep("Nucleus",colnames(vessels),value=T))],vessels$VesselID,statscols = c("mean","numel"))
    vesselstats=as.data.frame(vesselstats$mean)
    vesselstats[,STAT]=vesselstats[,STAT]/max(vesselstats[,STAT])
    vesselstats[,STAT]=round(1+100*(vesselstats[,STAT]-min(vesselstats[,STAT])))
    vesselstats=vesselstats[sort(vesselstats[,STAT],index.return=T)$ix,]
    col=1:(1.2*max(vesselstats[,STAT]))
    col=col[vesselstats[,STAT]]
    names(col)=as.character(vesselstats$VesselID)
    pdf(paste0(OUTD,filesep,fileparts(d_)$name,"_xy_",STAT,".pdf"),width = 10,height = 5)
    par(mfrow=c(1,2))
    plot(cells$`Centroid X µm`,cells$`Centroid Y µm`,pch=20,main=fileparts(v)$name,cex=0.1)
    points(-allcoord$`Centroid X µm`,allcoord$`Centroid Y µm`,pch=20,cex=0.4,col=col[as.character(allcoord$VesselID)]+allcoord$Boundary)
    plot(1,main="vessel thickness",xaxt = "n",yaxt = "n",axes=F,xlab="",ylab="",col="white"); 
    legend("topleft",as.character(unique(vesselstats[,STAT])),fill=unique(col),cex=0.75)
    dev.off()
    
    
    ## Save output
    allcoord$`Centroid X µm`=-allcoord$`Centroid X µm`
    cellscoord=cells[,intersect(colnames(allcoord),colnames(cells))]
    write.table(allcoord,file=paste0(OUTD,filesep,fileparts(d_)$name,"_vessels_xy.txt"),sep="\t",quote = F,row.names = F)
    cellscoord$karyotype=paste(rep(2,3),collapse=",")
    cellscoord$concentration=1
    
    ## Exclude negative cells overlapping blood vessels
    dtc=flexclust::dist2(cellscoord[,1:2],allcoord[allcoord$Boundary==0,1:2])
    dbv=flexclust::dist2(allcoord[allcoord$Boundary==1,1:2],allcoord[allcoord$Boundary==0,1:2])
    threshold = apply(dbv,2,min)
    dst=apply(dtc,2,min)
    dtc=dtc[,dst<threshold]
    threshold=threshold[dst<threshold]
    exclude=apply(dtc, 1, function(x) any(x<threshold))
    cellscoord = cellscoord[!exclude,]
    write.table(cellscoord,file=paste0(OUTD,filesep,fileparts(d_)$name,"_cells_xy.txt"),sep="\t",quote = F,row.names = F)
  }
  
  print(apply(stats,2,quantile,c(0,1)))
}


