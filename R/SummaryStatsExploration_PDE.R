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
}

## population growth rates
total = sapply(names(All$time_0_0), function(x) sapply(All, function(y) sum(sum(y[[x]],na.rm=T))))
par(mfrow=c(3,3)); 
tmp=sapply(colnames(total),function(x) plot(total[,x],xlab="time",ylab=x))


##Choose what (e.g. which compartment(s) or sum-stats) to visualize 
pop = list()
for(t in names(All)){
  ## guarantee consistent compartment order
  All[[t]] = All[[t]][names(All[[1]])]
  dom = NA*ones(max(xy$Var1),max(xy$Var2))
  dom = list(total=dom, max=dom,mean=dom, nonZeroMin=dom, nPop=dom)
  if(!excludeInitialComp){
    dom[[coi]] = dom$total
  }
  for(k in 1:nrow(xy)){
    i = xy[k,1]
    j = xy[k,2]
    focal= sapply(All[[t]], function(x) as.numeric(x[i,j]))
    print(focal)
    if(!excludeInitialComp){
      dom[[coi]][i,j] = focal[coi]
    }
    dom$total[i,j] = sum(focal,na.rm=T)
    dom$mean[i,j] = mean(focal,na.rm=T)
    dom$nPop[i,j] = sum(focal>0,na.rm=T)
    if(any(focal>0)){
      dom$nonZeroMin[i,j] = which.min(focal[focal>0])
      dom$max[i,j] = which.max(focal) 
    }
  }
  pop[[t]] = dom
}

mstest=any(sapply(pop, function(x) sum(!is.na(x$max) && x$max!=x$nonZeroMin)))
print(paste("Mis-segregations observed:",mstest))

##Visualize 
pdf("~/Downloads/ploidyModel_PDEeval.pdf",width = 10,height = 8)
for(what in names(dom)){
  par(mfrow=c(1,1))
  minmax=sapply(pop,function(x) quantile(as.numeric(x[[what]]), c(0,1),na.rm = T));
  plot(minmax[2,],log="y",ylab=what)
  minmax=quantile(minmax,c(0,1),na.rm=T)
  tmp=sapply(names(pop), function(x) image(pop[[x]][[what]],main=paste(what,x),col=rainbow(100)));#zlim=minmax,
}
tmp=sapply(Energy, function(x) image(as.matrix(x),main="Energy"))
dev.off()
