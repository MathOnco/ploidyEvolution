EnergyParameterization_dNTP_O2 <- function(cancer = "STAD"){
  options(java.parameters = "-Xmx7g")
  library(cloneid)
  library(xlsx)
  library("RMySQL")
  devtools::source_url("https://github.com/noemiandor/Utils/blob/master/grpstats.R?raw=TRUE")
  setwd("~/Projects/PMO/HighPloidy_DoubleEdgedSword/code/MathModel/GrowOrGo_PloidyEnergyRobustness/R")
  source("Utils.R")
  indir = "../data/"
  cellLines <- c("SNU-16","MKN-45", "NUGC-4","SNU-638", "KATOIII", "NCI-N87", "HGC-27", "SNU-668", "SNU-601")
  F_SDURATION = paste0(indir,"S-phase_duration_Params.xlsx")
  F_CCDURATION = paste0(indir,"Duration_PerCellCyclePhase_PerClone.txt")
  F_PLOIDYPOLYE= paste0(indir,"Ploidy_PolyE_PerCellCyclePhase_PerClone.txt")
  ########################
  ### Model parameters ###
  ########################
  params = read.xlsx(F_SDURATION, sheetIndex = "Params",as.data.frame = T)
  rownames(params)=params$Parameter
  P = as.data.frame(t(params[,"Value", drop=F]))
  UNITS = list(dNTP="mMol dNTP", o2=paste(params$Unit[params$Parameter=="K_M_RNR"],"O2"), po4="mmol_per_liter")
  if(cancer=="CL"){
    FLAG_SAMPLECOHORT = "Stomach Cell Lines"
    FLAG_SUBSTRATE = names(UNITS)[1]
  }else if(cancer=="STAD"){
    FLAG_SAMPLECOHORT = "Stomach Cancer"
    FLAG_SUBSTRATE = names(UNITS)[3]
  }else if(cancer=="GBM"){
    FLAG_SAMPLECOHORT = "Glioblastoma"
    FLAG_SUBSTRATE = names(UNITS)[2]
  }else{
    print(paste("cancer must be set to either 'STAD' or 'GBM'")) 
    return()
  }
  
  ############################
  ### Load cell lines info ###
  ############################
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  cellLineInfo = plotCellLineHistory()
  rownames(cellLineInfo) = cellLineInfo$name
  cellLineInfo = cellLineInfo[sort(cellLineInfo$doublingTime_hours,index.return=T)$ix,]
  cellLineInfo = cellLineInfo[cellLineInfo$name!="SNU-601", ]
  
  
  ########################
  ### Model parameters ###
  ########################
  cc_params = read.xlsx(F_SDURATION, sheetIndex = "O2_PO4_dNTP_CellCycle",as.data.frame = T)
  rownames(cc_params) = cc_params$CellCycle
  ##nondividing cells consume oxygen at a rate of:
  o2_nd =  P$O2_consumption_WT * 0.000001;  ## umol/cell/h
  # dividing cells consume 1.6 times more oxygen and display 3.2 times increased ATP levels compared to non-dividing
  o2_d = 1.6 * o2_nd; ## umol/cell/h
  ntp_nd = cc_params["G1",]$dNTP_mMol; ## dNTP consumed by non-dividing cell
  
  if(!file.exists(F_CCDURATION) || !file.exists(F_CCDURATION)){
    ##Estimate S-phase duration per clone from scDNA-seq data
    O = S_Phase_Duration_PerClone(cellLineInfo, outF=F_CCDURATION)
    ## Polymerase expression vs. ploidy
    O = PolymeraseExpression_Ploidy(cellLines, indir = indir, outF =  F_PLOIDYPOLYE)
    # O = read.table(  F_PLOIDYPOLYE, sep="\t", check.names = F, stringsAsFactors = F)
    O$cellLine = paste(match(O$cellLine, rownames(cellLineInfo)), O$cellLine); ## Enforce in order of doubling time
    te=cor.test(O$ploidy,O$polymerase_epsilon_delta)
    p <- ggplot(O, aes(x=ploidy, y =polymerase_epsilon_delta, col=cellLine, shape=state)) + ggtitle(paste(state,": r =", round(te$estimate,3),"; p =", round(te$p.value,10))) +
      geom_point(size=3) + geom_smooth(method='lm', formula= y~x, se=FALSE) + labs(x="ploidy", y = "polymerase_epsilon_delta") 
    ggsave(p, filename = paste0("~/Downloads/",state,"_gastricCLs_ploidy_vs_polymerase_epsilon_delta.png"), width = 5, height = 4.5) 
  }
  
  ## Read ploidy and cell cycle duration info
  if(FLAG_SAMPLECOHORT == "Stomach Cell Lines"){
    # Cell lines
    tmp = read.table(  F_CCDURATION, sep="\t", check.names = F, stringsAsFactors = F)
    cpc = read.table(  F_PLOIDYPOLYE, sep="\t", check.names = F, stringsAsFactors = F)
    cpc$Identity = sapply(strsplit(rownames(cpc),".", fixed = T),"[[",3)
    cpc = cpc[cpc$state=="G0G1" & cpc$Identity %in% tmp$Identity,]
    cpc$S_Duration_hours = tmp$S_Duration_hours[match(cpc$Identity, tmp$Identity)]
    cpc$ploidy = round(cpc$ploidy,2)
    cpc = cpc[sort(cpc$ploidy, index.return=T)$ix,]
    cpc = cpc[cpc$cellLine %in% c("MKN-45","SNU-16"),]
    cpc = as.data.frame(grpstats(cpc[,c("ploidy", "S_Duration_hours")],cpc$cellLine, statscols = "mean")$mean)
    cpc$Sample = rownames(cpc)
    
    doi = cc_params["S",]$dNTP_mMol; #dNTP_sPhase_mMol; 
    col = gg_color_hue(nrow(cellLineInfo));
    names(col) = rownames(cellLineInfo);
    # x_dNTP=(50:(5000*P$K_M_poly*18))/5000
    x_dNTP=seq(5E-5,6.5E-3,by=1E-6)
  }else { ##if(FLAG_SAMPLECOHORT = "Glioblastoma")
    # P+R GBM
    PATIENT = "P06269"
    cpc = read.table(paste0(indir,PATIENT,"_ploidyPerClone.txt"), sep="\t", check.names = F, stringsAsFactors = F )
    colnames(cpc) = "ploidy"
    cpc$S_Duration_hours = 7
    cpc$Sample = paste("Clone",sapply(strsplit(rownames(cpc),"ID"),"[[", 2))
    cpc$Size = sapply(strsplit(rownames(cpc),"_"),"[[", 2)
    cpc = cpc[sort(cpc$Size, decreasing = T, index.return=T)$ix,]
    cpc$ploidy = c(4,3,2)
    
    # col = gg_color_hue(nrow(cpc));
    col = c("darkolivegreen3","yellow","darkred")
    names(col) = cpc$Sample
    cpc = cpc[sort(abs(cpc$ploidy-2), decreasing = T, index.return=T)$ix,]
    doi = P$O2_brain_min
    
    # x=(30:(10000*P$K_M_poly*3))/5000
    x_o2 = seq(0.75,P$K_M_RNR*3.5,by=0.01)
  }
  
  
  #############################################################################################
  #### Visualize dependency of DNA replication speed on substrate concentration and ploidy ####
  #############################################################################################
  activeRNRPerCell=1
  out = list()
  logF = "x"
  # pdf(paste0("~/Downloads/",FLAG_SUBSTRATE,".pdf"), width = 4, height = 4)
  for(i in 1:nrow(cpc)){
    ploidy = cpc$ploidy[i]
    s_Hours = mean(cpc$S_Duration_hours[cpc$Sample == cpc$Sample[i]])
    humanGenome = ploidy*P$humanGenome 
    numPolymerases = 1/(s_Hours*P$v_max_poly/(2000*P$humanGenome)); ## human cells have on the order of 100,000 origins of replication. So depending on how many origins are active simultaneously, there are likely many thousands of polymerase molecules acting at once to replicate the DNA of a single eukaryotic cell.
    
    if(FLAG_SUBSTRATE=="o2"){
      ## Replication time vs. O2 levels
      x_ = x_o2
      x_dNDP = activeRNRPerCell* P$v_max_RNR * x_o2 / (ploidy*P$K_M_RNR + x_o2)
      y = numPolymerases*P$v_max_poly * x_dNDP / (ploidy*P$K_M_poly + x_dNDP)
    }else{ 
      y = numPolymerases*P$v_max_poly * x_dNTP / (ploidy*P$K_M_poly + x_dNTP)
      if(FLAG_SUBSTRATE=="po4"){
        ## Replication time vs. PO4 levels
        dNTP = x_dNTP*numPolymerases 
        x_ = (cc_params["S",]$PO4_mMol_per_liter - max(cc_params[c("M","G1"),]$PO4_mMol_per_liter))* dNTP/800; ## conversion from XX
      }else  if(FLAG_SUBSTRATE=="dNTP"){
        ## Replication time vs. dNTP levels
        x_ = x_dNTP
      } 
    }
    out[[cpc$Sample[i]]] = humanGenome/y
  } 
  
  
  plot(x_, out[[cpc$Sample[1]]], ylim=quantile(do.call("c",out),c(0,1)),ylab="hours to replicate genome",log=logF,col="white",xlab=UNITS[[FLAG_SUBSTRATE]])
  # plot(x_, out[[cpc$Sample[i]]], ylim=quantile(c(hours,hours2),c(0,1)),ylab="hours to replicate genome",log=logF,col="white")
  for(i in 1:length(out)){
    lines(x_, out[[cpc$Sample[i]]], col=col[cpc$Sample[i]], lwd=5)
    points(x_, out[[cpc$Sample[i]]], cex=0.4, col=col[cpc$Sample[i]], pch=20)
    ## plot single data point
    idx = which.min(abs(x_-doi))
    if(idx>1 & idx<length(x_)){
      print(humanGenome/(numPolymerases* y[idx]))
      points(x_[idx], out[[cpc$Sample[i]]][idx], pch=3, cex=2)
    }
  }
  legend("topright",paste("ploidy", round(cpc$ploidy,2),"(",cpc$Sample,")"), fill=col[cpc$Sample],bty="n",cex=0.7)
  mtext(doi)
  
  
  
  
  
  # ## Map onto risk of death
  # z = (humanGenome/(P$N_poly_Hela* y))/100
  # z[z>1] = 1
  # plot(o2_mmhg, z, xlab="O2 mmHg", ylab="risk of death",log="x",pch=20,col=col[as.character(ploidy)])
  # plot(o2_mmhg, 1-z, xlab="O2 mmHg", ylab="chance of complete mitosis",log="x",pch=20,col=col[as.character(ploidy)])
  
  
}
