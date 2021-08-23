setwd("~/Downloads/ploidyEvolution/Julia")
name="cn_5_seq2.csv"
coi="#1:1:1:1:1"
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
  image(cells_per_cptm$`#1:1:1:1:1`)
  sapply(cells_per_cptm,dim)
  All[[time]]=cells_per_cptm
}

##Choose what (e.g. which compartment(s) or sum-stats) to visualize 
pop = list()
for(t in names(All)){
  ## guarantee consistent compartment order
  All[[t]] = All[[t]][names(All[[1]])]
  dom = ones(max(xy$Var1),max(xy$Var2))
  dom = list(total=dom, max=dom)
  dom[[coi]] = dom$total
  for(k in 1:nrow(xy)){
    i = xy[k,1]
    j = xy[k,2]
    focal= sapply(All[[t]], function(x) as.numeric(x[i,j]))
    dom$max[i,j] = which.max(focal) 
    dom[[coi]][i,j] = focal[coi]/sum(focal)
    dom$total[i,j] = sum(focal)
  }
  pop[[t]] = dom
}

##Visualize
pdf("~/Downloads/ploidyModel_PDEeval.pdf",width = 10,height = 8)
for(what in names(dom)){
  par(mfrow=c(2,3))
  sapply(pop, function(x) image(x[[what]],main=what))
  sapply(Energy, function(x) image(as.matrix(x),main="Energy"))
}
dev.off()

