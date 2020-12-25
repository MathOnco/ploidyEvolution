.connect2DB <-function(){
  tmp = suppressWarnings(try(lapply( dbListConnections( dbDriver( drv = "MySQL")), dbDisconnect)))
  yml = yaml::read_yaml(paste0(system.file(package='cloneid'), '/config/config.yaml'))
  mydb = dbConnect(MySQL(), user=yml$mysqlConnection$user, password=yml$mysqlConnection$password, dbname=yml$mysqlConnection$database,host=yml$mysqlConnection$host, port=as.integer(yml$mysqlConnection$port))
  return(mydb)
}


#######################################################
### Interpolate and visualize birth rate landscape ####
interpolateBirthRateLandscape <- function(cnv, r, split_train_test=1, ndim = 4, save=NULL, nrep = 1, usePCA=T, ipolmethod = "polyharmonic", markState=NULL){
  ## PCA
  pcs <- FactoMineR::PCA(cnv, graph = T,scale.unit = T,ncp = ncol(cnv))
  pcs$x = pcs$ind$coord
  # print(factoextra::get_eigenvalue(pcs)[,"variance.percent"])
  barplot(sapply(1:ncol(pcs$x), function(x) sum(factoextra::get_eigenvalue(pcs)[1:x,"variance.percent"])), ylab="cummulative variance explained", names.arg = 1:ncol(pcs$x))
  pct_explained = sweep(pcs$var$contrib, 2, FUN="*", factoextra::get_eigenvalue(pcs)[,"variance.percent"]/100)
  pct_explained = apply(pct_explained[,1:3], 1, sum); 
  pct_explained = pct_explained[names(sort(pct_explained, decreasing = T))]
  par(mai=c(0.95,2.75,0.25,0.25))
  barplot(t(pcs$var$contrib[names(pct_explained),1:3]), las=2, beside = F, horiz = T, xlab="% variance explained");
  legend("topright", c("PC1", "PC2","PC3"), fill = gray.colors(3) )
  ## Scatter plot
  scatterplot3d::scatterplot3d(pcs$x[names(r),1], pcs$x[names(r),2], r, pch=20)
  # rgl::plot3d(pcs$x[,1], pcs$x[,2], r, pch=20)
  
  ## Interpolate PC components or raw cnv data?
  data = pcs$x[,1:ndim]
  if(!usePCA){
    data = cnv[,1:ndim]
  }
  
  ## Interpolation
  te = matrix(NA,nrep,2)
  colnames(te) = c("r","P")
  for(m in 1:nrep){
    ##Test vs. training
    ij <- ii <- names(r)
    if(split_train_test<1){
      ii = sample(ii, length(ii)*split_train_test); # train
      ij = setdiff(names(r), ii)                    # test
    }
    
    ## Multivariate interpolation
    if(ndim >2 ){
      ipol_birth <- ipol(r[ii], knots = t(data[ii,]), method=ipolmethod)
      r_pred = try(ipol_birth(t(data[ij,])))
      
      ## Prep for visualization:
      RES = 2*length(ii)-1
      grid <- apply(data[ii,], 2, function(x) seq(min(x),max(x),by=(range(x)[2]-range(x)[1])/RES))
      grid = split(grid, rep(1:ncol(grid), each = nrow(grid)))
      ## Modify real data with different versions of altered CNs:
      mat=list()
      cnt = 1;
      for(dim1 in 1:(ncol(data)-1)){
        for(dim2 in (dim1+1):ncol(data)){
          mat_=list(x=grid[[dim1]], y=grid[[dim2]], z=matrix(0,RES+1,RES+1), xlab=colnames(data)[dim1], ylab=colnames(data)[dim2])
          lm = expand.grid(1:length(mat_$x), 1:length(mat_$y))
          xy = as.matrix(expand.grid(mat_$x, mat_$y))
          for(k in 1:nrow(lm)){
            synthetic = data[ii,]
            synthetic[,c(dim1,dim2)] = repmat(xy[k,], nrow(synthetic),1)
            mat_$z[lm[k,1],lm[k,2]] = max(apply(synthetic, 1, ipol_birth))
          }
          
          mat[[cnt]] = mat_
          ## Stop at 9
          cnt = cnt+1;
          if(cnt>=9){
            break;
          }
        }
      }
    } else{
      ## 3D Interpolation
      mat = akima::interp(data[ii,1], data[ii,2], r[ii], duplicate = "median")
      ## Goodness of fit:
      idx_x = apply(flexclust::dist2(data[ij,1],mat$x), 1, which.min) 
      idx_y = apply(flexclust::dist2(data[ij,2],mat$y), 1, which.min)
      idx = cbind(idx_x, idx_y)
      r_pred = apply(idx, 1, function(i) mat$z[i[1],i[2]])
      mat$xlab=colnames(data)[1];
      mat$ylab=colnames(data)[2];
      mat = list(mat); ##wrap
    }
    if(!is.null(save)){
      pdf(save);
    }
    par(mfrow=c(3,3))
    for(mat_ in mat){
      fields::image.plot(mat_, xlab=mat_$xlab, ylab=mat_$ylab,legend.lab = "birth rate", legend.cex = 1, horizontal = T)
      ## Mark initial condition on tensor:
      if(!is.null(markState)){
        points(markState[,mat_$xlab], markState[,mat_$ylab], pch=18, col="white", cex=4)
      }
    }
    if(!is.null(save)){
      mtext(matlab::fileparts(save)$name)
      dev.off()
    }
    rgl::rgl.surface(mat_$x,mat_$y,mat_$z,color="green",alpha=c(0.5)) 
    
    te_ = try(cor.test(r[ij],r_pred))
    if(class(te_)=="try-error"){
      te_ = list(estimate=NA, p.value=NA)
    }
    te[m,] = c(te_$estimate, te_$p.value)
  }
  return(list(pcs=pcs, te=te, pct_explained=pct_explained,mat=mat, ml2=ipol_birth))
}


############################################
#### Estimate S-phase duration per clone ###
############################################
S_Phase_Duration_PerClone <- function(cellLineInfo, outF="Duration_PerCellCyclePhase_PerClone.txt"){
  
  THESTATES=c("G0G1","apoptotic","S")
  allCPC = list()
  for (cellLine in rownames(cellLineInfo)){
    print(paste("Processing", cellLine,"..."))
    closeAllConnections()
    mydb = .connect2DB()
    i2p   = cloneid::identity2perspectiveMap(cellLine,persp = "GenomePerspective")
    
    ## Gather cycling cells and G0G1 cells for each identity
    membership = c()
    for (x in names(i2p)){
      rs = dbSendQuery(mydb, paste0("select cloneid from Perspective where parent=",x,";"))
      kids = unlist(fetch(rs, n=-1))
      m = rep(x, length(kids))
      names(m) = kids
      membership =c(membership,m)
    }
    states=sapply(as.numeric(names(membership)),getState, "GenomePerspective");
    names(states) = names(membership)
    
    a=plyr::count(c(unique(membership),membership[states==THESTATES[1]]));   rownames(a)=as.character(a$x)      
    b=plyr::count(c(unique(membership),membership[states==THESTATES[2]]));   rownames(b)=as.character(b$x); 
    c=plyr::count(c(unique(membership),membership[states==THESTATES[3]]));   rownames(c)=as.character(c$x); 
    cpc=cbind(a$freq, b$freq, c$freq, (b$freq+a$freq) );
    colnames(cpc)=c(THESTATES,"noncycling"); 
    rownames(cpc)=rownames(a)
    cpc[is.na(cpc)]=0;    
    cpc=as.data.frame(cpc);
    cpc=sweep(cpc,MARGIN = 2, FUN="/", apply(cpc,2,sum));
    cpc$cellLine = cellLine
    cpc$Identity = i2p[rownames(cpc)]
    cpc$GenomePerspective = rownames(cpc)
    cpc$S_Duration_hours = cpc$S/(cpc$S+cpc$noncycling)*cellLineInfo[cellLine,]$doublingTime_hours
    allCPC[[cellLine]] = cpc
    print(cpc)
  }
  write.table(do.call(rbind, allCPC),  outF, sep="\t", quote=F)
  return(allCPC);
}


#########################################
#### Polymerase expression vs. ploidy ###
#########################################
PolymeraseExpression_Ploidy <- function(cellLines, indir, outF =  "Ploidy_PolyE_PerCellCyclePhase_PerClone.txt"){
  ## DNA polymerase epsilon is a member of the DNA polymerase family of enzymes found in eukaryotes. 
  ## It is composed of the following four subunits: POLE (central catalytic unit), POLE2 (subunit 2), POLE3 (subunit 3), and POLE4
  ## DNA polymerase delta (DNA Pol δ) is an enzyme complex found in eukaryotes that is involved in DNA replication and repair. 
  ## The DNA polymerase delta complex consists of 4 subunits: POLD1, POLD2, POLD3, and POLD4.
  ## DNA Pol δ is an enzyme used for both leading and lagging strand synthesis
  goi = c("POLD1", "POLD2", "POLD3", "POLD4", "POLE", "POLE2", "POLE3", "POLE4")
  O = list()
  for (sName in cellLines){
    for(state in c("S","G0G1")){
      # Gather genome- and transcriptome perspective per clone first if data not there already:
      if(!file.exists(paste0(indir,sName,"_",state,"_",x,".txt"))){
        gatherGenome_TranscriptomePerspective_PerClone(sName, outdir = indir)
      }
      fr = lapply(c(DNA,RNA), function(x) read.table(paste0(indir,sName,"_",state,"_",x,".txt"),sep="\t", check.names = F, stringsAsFactors = F))
      names(fr) = c(DNA,RNA)
      loci = cloneid::parseLOCUS(rownames(fr$GenomePerspective))
      ploidy = sapply(colnames(fr$GenomePerspective), function(x) calcPloidy(cbind(loci,fr$GenomePerspective), cnColumn = x))
      polyE = apply(fr$TranscriptomePerspective[goi,],2,mean, na.rm=T)
      sName_ = paste0(state, ".", sName)
      O[[sName_]] = as.data.frame(cbind(ploidy,polyE[names(ploidy)]))
      colnames(O[[sName_]]) = c("ploidy","polymerase_epsilon_delta")
      plot(O[[sName_]]$ploidy, O[[sName_]]$polymerase_epsilon_delta, pch=20, main=sName_)
    }
  }
  O = do.call(rbind,O)
  O$state = sapply(strsplit(rownames(O),".", fixed = T),"[[",1)
  O$cellLine = sapply(strsplit(rownames(O),".", fixed = T),"[[",2)
  write.table(O, outF, sep="\t", quote=F)
  return(O)
}

################################################################
#### Gather genome- and transcriptome perspectives per clone ####
################################################################
gatherGenome_TranscriptomePerspective_PerClone <- function(sName, outdir){
  print(paste("Reading genome- and transcriptome perspectives for", sName, "..."))
  
  rc=cloneid::getSubclones(sName,RNA)
  r=gatherSCprofiles(rc, ccState=state,whichP = RNA); 
  commonLoci=grep(":",rownames(r$X),value = T);
  
  dc=cloneid::getSubclones(sName,DNA)
  d=gatherSCprofiles(dc, ccState=state,whichP = DNA); 
  
  scseq = list(d,r)
  names(scseq) = c(DNA,RNA)
  identities=cloneid::getSubclones(sName,"Identity")
  pOnId=list()
  for(assay in c(DNA,RNA)){
    ##Perspective on Identity
    pOnId[[assay]]=sapply(identities,function(x) x$getPerspective(J("core.utils.Perspectives")$valueOf(assay))$toString()); 
    pOnId[[assay]] = pOnId[[assay]][extractID(pOnId[[assay]]) %in% scseq[[assay]]$clones]
  }
  
  ## From single cell to Clone resolution
  ids = intersect(names(pOnId[[1]]), names(pOnId[[2]]))
  fr = lapply(scseq, function(o) grpstats(t(o$X), o$clones,'mean')$mean)
  names(fr) = names(scseq)
  fr = sapply(names(fr), function(assay) fr[[assay]][extractID(pOnId[[assay]][ids]),])
  rownames(fr[[2]]) <- rownames(fr[[1]]) <- extractID(ids)
  sapply(names(fr), function(x) write.table(t(fr[[x]]), file=paste0(outdir, filesep, sName,"_",state,"_",x,".txt"), sep="\t", quote = F))
}

