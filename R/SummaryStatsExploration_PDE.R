library(matlab)
setwd("~/Downloads/sim/")
name="cn_3_seq3"
excludeInitialComp = F
coi="#2:2:2"
f=list.files(".",pattern=name)
xy=c(grep("x.csv",f,value = T),grep("y.csv",f,value = T))
xy=sapply(xy,function(x) read.table(x)$V1)
xy =expand.grid(1:nrow(xy),1:nrow(xy))

## Read cell states for each timepoint
f=grep("states_time",f,value = T)
All=list()
Energy=list()
compartments = c()
for(csv in f){
  time=strsplit(fileparts(csv)$name,"_states_")[[1]][2]
  ## Read energy
  Energy[[time]]=read.csv(gsub("_states_","_Energy_",csv))
  ##Test only: for row-column consistency
  la=read.csv(csv,comment.char = "#")
  ## Read compartment info
  a=read.csv(csv,fill = T,header = F)
  idx = grep("#",a$V1)
  cptm = a$V1[idx]
  idx = c(idx, nrow(a)+1); ## last comparment not surrounded by two compartment names
  idx = lapply(2:length(idx), function(i) (idx[i-1]+1):(idx[i]-1) )
  cells_per_cptm = lapply(idx,function(i) apply(a[i,],2,as.numeric))
  names(cells_per_cptm) = cptm
  image(cells_per_cptm[[coi]])
  sapply(cells_per_cptm,dim)
  if(excludeInitialComp){
    cells_per_cptm=cells_per_cptm[names(cells_per_cptm) != coi ]
  }
  All[[time]]=cells_per_cptm
  compartments = c(compartments,names(cells_per_cptm))
}
compartments = unique(compartments)

# ## Populate empty compartments with 0:
# for (time in names(All)){
#   missing= setdiff(compartments,names(All[[time]]))
#   for (compartment in missing){
#     All[[time]][[compartment]]=zeros(max(xy$Var1),max(xy$Var2))
#   }
# }

## Sort by time
times = as.numeric(gsub("_",".",gsub("time_","",names(All))))
All=All[sort(times,index.return=T)$ix]
times = as.numeric(gsub("_",".",gsub("time_","",names(All))))

## population growth rates
total = sapply(compartments, function(x) sapply(All, function(y) sum(sum(y[[x]],na.rm=T))))
pdf("~/Downloads/ploidyModel_Total_ByCompartment.pdf")
par(mfrow=c(3,3)); 
tmp=sapply(colnames(total),function(x) plot(times,total[,x],xlab="time",ylab=x))
dev.off()

##Choose what (e.g. which compartment(s) or sum-stats) to visualize 
pop = list()
for(n in 3:length(All)){
  t = names(All)[n]
  t_ = names(All)[n-1]
  ## guarantee consistent compartment order
  All[[t]] = All[[t]][names(All[[1]])]
  dom = NA*ones(max(xy$Var1),max(xy$Var2))
  dom = list(total=dom, max=dom,mean=dom, nPop=dom, maxGrowth=dom, howManyGrow=dom)
  for(state in names(All[[1]])){
    dom[[state]] = dom$total
  }
  dom[[paste0("N_",coi)]] = dom$total 
  for(k in 1:nrow(xy)){
    i = xy[k,1]
    j = xy[k,2]
    focal= sapply(All[[t]], function(x) as.numeric(x[i,j]))
    focal_= sapply(All[[t_]], function(x) as.numeric(x[i,j]))
    ## Register growth of each compartment
    for(state in names(focal)){
      dom[[state]][i,j] = log(focal[state]/ focal_[state])
      if(focal_[state]==0){
        dom[[state]][i,j] =NA 
      }
    }
    dom$total[i,j] = sum(focal,na.rm=T)
    dom$mean[i,j] = mean(focal,na.rm=T)
    ## How many compartments co-exist at this gridpoint
    dom$nPop[i,j] = sum(focal>0,na.rm=T)
    ## How many representatives of our compartment of interest
    dom[[paste0("N_",coi)]][i,j] = focal[coi]
    if(any(focal>0)){
      ## Index of compartment which grows fastest
      dom$maxGrowth[i,j] =  which.max(sapply(names(focal), function(x) dom[[x]][i,j]))
      ## How many compartments grow significantly
      dom$howManyGrow[i,j] =  sum(sapply(names(focal), function(x) dom[[x]][i,j]>0.05),na.rm=T)
      ## Index of compartment which wins here:
      dom$max[i,j] = which.max(focal) 
    }
  }
  pop[[t]] = dom
}


##Visualize overall
col=rainbow(100)[1:80]
f <- function(dat, what,minmax) {
  image(dat,main=what,col=col,zlim=minmax)
  oceanmap::set.colorbar(cbx=c(0, 1), cby=c(-.05, -.08),cex.cb.ticks = 0.74,dat[!is.na(dat)],pal=col,zlim = minmax);#,breaks=seq(minmax[1],minmax[2],by=(minmax[2]-minmax[1])/length(col)))
}
show<-function(pop,what){
  minmax=sapply(pop,function(x) quantile(as.numeric(x[[what]]), c(0,1),na.rm = T));
  plot(minmax[2,],log="y",ylab=what)
  minmax=quantile(minmax,c(0,1),na.rm=T)
  tmp=sapply(names(pop), function(x) f(pop[[x]][[what]],paste(what,x),minmax));
  return(minmax)
}
pdf("~/Downloads/ploidyModel_PDEeval.pdf",width = 6,height = 5)
for(what in setdiff(names(pop$time_0_2),compartments)){ #c("total","nPop","growth")
  par(mfrow=c(1,1))
  try(show(pop,what))
}
minmax=quantile(unlist(Energy),c(0,1))
minmax[minmax==0] = 1E-3
tmp=sapply(names(Energy), function(x) f(as.matrix(log(Energy[[x]])),paste("log Energy",x),minmax = log(minmax)))
dev.off()

##Visualize by compartment
minmax=sapply(pop, function(x) quantile(unlist(x[compartments]),c(0,1),na.rm=T))
minmax = quantile(minmax,c(0,1))
pdf("~/Downloads/ploidyModel_PDEeval_byComp.pdf",width = 7,height = 8)
for(t in names(pop)){
  par(mfrow=c(3,3))
  sapply(compartments, function(x) try(f(pop[[t]][[x]],what=paste(x,t),minmax)))
}
dev.off()

# QOI: plot quantity of interest
gifski::save_gif(show(pop,"howManyGrow"),gif_file = "~/Downloads/Firsttest.gif",height = 600,width = 600)


# Compare low vs. high energy total growth
times=as.numeric(gsub("_",".",gsub("time_","",names(pop))))
tmp=as.matrix(Energy$time_0_1); 
ii = which(Energy$time_0_1<1)
iMax = ii[which.max(tmp[ii])]
ii = which(Energy$time_0_1>0.05)
iMin = ii[which.min(tmp[ii])]
## Energy @ selected locations
tmp[iMax]=NA; tmp[iMin]=NA;
image(log(tmp))
for(what in names(pop$time_0_2)){
  ## Cells @ selected locations:
  high_low= sapply(names(pop), function(x) pop[[x]][[what]][c(iMax,iMin)])
  if(all(is.na(high_low))){
    next;
  }
  plot(times,high_low[1,], ylim=quantile(high_low,c(0,1),na.rm=T), pch=20, col="red",ylab=what)
  points(times,high_low[2,], pch=20, col="blue")
  legend("topleft",c("high energy","low energy"), fill=c("red","blue"))
}