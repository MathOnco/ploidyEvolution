setwd("C:/Users/4473331/Documents/projects/ploidyEvolution-main/Julia/")

varyer <- function(par,val){
  x <- readLines("input_spatial.toml")
  
  id <- grepl(par,x)
  if(par=="difc") id <- grepl("\"Î"\"",x) 
  if(par=="difE") id <- grepl("\"Î"â,'\"",x)   
  if(par=="Ethresh") id <- grepl("\"Ï.\"" ,x)
  if(par=="taxis") id <- grepl("\"Îz\"" ,x)  
  if(par=="consumption") id <- grepl("\"Î´\"",x) 
  
  
  y <- strsplit(x[id],split="=")[[1]]
  
  new_y <- as.numeric(y[2])*val
  new_y <- paste(y[1],new_y,sep="= ")
  x[id] <- new_y
  writeLines(x,"input_tmp.toml")
  system("julia --project=. -L ploidyMovement.jl driver.jl -m -s -i input_tmp.toml")
  yname <- as.numeric(y[2])*val
  nm <- paste(par,yname,sep="_")
  print(paste("test_output/ABC/sweep_00/",nm,"_population.csv",sep=""))
  file.copy(from = "0300_population.csv",
            to   = paste("test_output/ABC/sweep_00/",nm,"_population.csv",sep=""))
  
  file.copy(from = "0300_pdestates.csv",
            to   = paste("test_output/ABC/sweep_00/",nm,"_pdestates.csv",sep=""))
  
  file.copy(from = "0300_stochstates.csv",
            to   = paste("test_output/ABC/sweep_00/",nm,"_stochstates.csv",sep=""))
  print(paste("test_output/ABC/sweep_00/",nm,"_stochstates.csv",sep=""))
}

wrapper <- function(par,vals=c(0.1,0.2,0.4,0.8,1.6,2.4,4.8)){
  lapply(vals,function(val) varyer(par,val))
}

dir.create("test_output/ABC/sweep_00")
pars <- c("deathRate","misRate","difc","difE","Ethresh","taxis","consumption","k")


lapply(pars,wrapper)
                            
          
