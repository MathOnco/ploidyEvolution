chrmod <- function(time,state,parms){
  with(as.list(parms),{
    
    ds <- q%*%state
    
    return(list(ds))
  })
}

## generate the transition matrix (q) for the model
generate_q <- function(eta,delta,B,D,maxchrom){
  
  lambda <- 1
  q <- matrix(0,maxchrom,maxchrom)
  for(i in 1:maxchrom){
    for(j in 1:maxchrom){
      q[i,j] <- lambda*choose(j,abs(i-j))*B(eta,j,k=maxchrom)^abs(i-j)*(1-B(eta,j,k=maxchrom))^(j-abs(i-j))
      if(i==j){
        qij <- choose(j,abs(i-j))*B(eta,j,k=maxchrom)^abs(i-j)*(1-B(eta,j,k=maxchrom))^(j-abs(i-j))
        ## equivalently
        qij <- 1*1*(1-B(eta,j,k=maxchrom))^(j)
        q[i,j] <- lambda*(2*qij-1)-D(delta,i,k=maxchrom)
      }
    }}
  
  return(q)
}


## determine whether given parameters lie on the critical curve (occurs when function returns 0)
c_crit <- function(delta,eta,B,maxchrom,D,very_small_flag=FALSE){
  
  if(very_small_flag==TRUE){
    delta <- 10^delta
  }
  ## generate transition matrix
  q <- generate_q(eta,delta,B,D,maxchrom)
  z <- tryCatch(eigen(q)$value, error = function(e) return("error"))
  if(z[1]=="error") return(10^5)
  res <- 10^5
  if(sum(Im(z)==0)>0){
    res <- abs(max(Re(z[Im(z)==0])))
  }
  
  return(res)
}

## find the value of delta lying closest to the critical curve
find_c_crit <- function(eta,B,maxchrom,D){
  id <- D(NaN,NaN,label=TRUE)
  interval <- c(0,1)
  very_small_flag<-FALSE
  if(id=="D(i)==1-1/(mu+e^sqrt( abs( i-3 ) ))") interval=c(0,1000)
  if(id=="log(i*mu)^2/K^2"){
    interval=c(-20,0)
    very_small_flag <- TRUE
  }
  if(B(NaN,NaN,TRUE)=="B(i)==1/(beta+e^sqrt( abs( i-3 ) ))" & id==" D(i)==mu"){
    interval=c(-10,0)
    very_small_flag <- TRUE
  }
  if(B(NaN,NaN,TRUE)==" B(i)==beta" & id==" D(i)==mu"){
    interval=c(-10,0)
    very_small_flag <- TRUE
  }
  out <- optimise(c_crit,interval=interval,tol=0.000001,eta=eta,B=B,maxchrom=maxchrom,D=D,very_small_flag=very_small_flag)
  res <- data.frame(v=out$objective,delta=out$minimum,eta=eta)
  if(very_small_flag==TRUE){
    res <- data.frame(v=out$objective,delta=10^out$minimum,eta=eta)
  }
  return(res)
}

find_c_crit_tran <- function(eta,delta,B,D,maxchrom){
  
  q <- generate_q(eta,delta,B,D,maxchrom)
  z <- eigen(q)$value
  v <- eigen(q)$vectors[,which(Re(z)==max(Re(z[Im(z)==0])))]
  v <- Re(v/sum(v))
  net_g <- v*(1-sapply(1:length(v), function(i) D(delta,i,k=maxchrom)))
  net_m <- q%*%v
  extra_D <- sum(net_g-net_m)
  pop_B <- sapply(1:length(v), function(j) B(eta=eta,j,k=maxchrom))
  pop_B <- 1-(1-pop_B)^(1:length(v)) 
  pop_D <- sapply(1:length(v), function(j) D(delta=delta,j,k=maxchrom))
  pop_B <- sum(v*pop_B)
  pop_D <- sum(v*pop_D)#+extra_D
  
  data.frame(pop_B=pop_B,pop_D=pop_D,eta=eta,delta=delta)
}


assess_model <- function(B,D,maxchrom,granularity=50){
  print(paste("evaluating model... ", B(NaN,NaN,TRUE), D(NaN,NaN,TRUE) ))
  
  starters <- c(0.00001,0.00002,0.0000225,0.0001,0.001,0.002,0.004,0.008)
  enders <- 1-starters
  
  eta <- c(starters,seq(1/granularity,1-1/granularity,1/granularity),enders)
  
  if(B(NaN,NaN,TRUE)=="B(i)==1/(beta+e^sqrt( abs( i-3 ) ))"){
    eta <- 2000/2^seq(0,20,20/granularity)
    eta <- c(8:1*10000,eta)
  }
  
  if(B(NaN,NaN,TRUE)=="B(i)==e^(-beta*i)"){
    eta <- c(eta,seq(1,3,0.25))
  }
  

  df <- do.call(rbind,lapply(eta,find_c_crit,B,maxchrom,D))
  df <- df[df$v<0.01,]
  x <- lapply(1:nrow(df), function(i){
    find_c_crit_tran(df$eta[i],df$delta[i],B,D,maxchrom)
  })
  x <- do.call(rbind,x)
  x <- x[x$pop_B<1&x$pop_D<1,]
  
  x$B <- B(NaN,NaN,TRUE)
  x$D <- D(NaN,NaN,TRUE)
  x$maxchrom <- str_pad(maxchrom,width=2)
  return(x)
}

assess_all_models <- function(model,maxchrom,granularity=100){
  x <- do.call(rbind,lapply(maxchrom, function(mc) assess_model(model$B,model$D,mc,granularity)))
  return(x)
}

run_model <- function(eta,delta,B,D,maxchrom,t_end){
  #model <- 2
  lambda <- 1
  #maxchrom <- 8
  state <- rep(0,maxchrom)
  state[2] <- 1
  times <- seq(0,t_end,t_end/1000)
  
  q <- generate_q(eta,delta,B,D,maxchrom)
  
  parms <- list(q=q)
  out <- ode.1D(y = state, times = times, func = chrmod, parms = parms,
                nspec = 1, names = c("CHR"),method="vode")
  return(out)
  
}

wrap_run_model <- function(eta,delta,maxchrom,model){
  B <- model$B
  D <- model$D
  t_end <- 100 
  out <- run_model(eta,delta,B,D,maxchrom,t_end)
  state <- out[nrow(out),2:ncol(out)]
  state <- state/sum(state)
  
  
  pop_B <- sapply(1:length(state), function(j) B(eta=eta,j,k=maxchrom))
  pop_B <- 1-(1-pop_B)^(1:length(state))
  pop_D <- sapply(1:length(state), function(j) D(delta=delta,j,k=maxchrom))
  pop_B <- sum(state*pop_B)
  pop_D <- sum(state*pop_D)
  
  data.frame(eta=eta,delta=delta,pop_B=pop_B,pop_D=pop_D,maxchrom=str_pad(maxchrom,width=2))
}

process_cloud <- function(out,B,D){
  eta_u <- unique(out$eta)
  delta_u <- unique(out$delta)
  
  if(B==" B(i)==beta"&D=="D(i)==1-1/(mu+e^sqrt( abs( i-3 ) ))"){
    out <- out[out$delta==head(delta_u,1),]
    return(out)
  }
  
  if(B=="B(i)==1/(beta+e^sqrt( abs( i-3 ) ))"&D=="D(i)==1-1/(mu+e^sqrt( abs( i-3 ) ))"){
    out <- out[out$eta==eta_u[1]|out$delta==delta_u[1],]
    return(out)
  }
  
  if(B=="B(i)==1/(beta+e^sqrt( abs( i-3 ) ))"&D==" D(i)==mu"){
    out <- out[out$eta==eta_u[1],]
    return(out)
  }
  
  if(D=="log(i*mu)^2/K^2") return(out)
  
  if(D=="D(i)==1-1/(mu+e^sqrt( abs( i-3 ) ))"){
    out <- out[out$delta==delta_u[1],]
    return(out)
  }
  
}

grid_search_models <- function(model,maxchrom=c(4,8,16,32)){
  B <- model$B(NaN,NaN,label=TRUE)
  D <- model$D(NaN,NaN,label=TRUE)
  print(B)
  print(D)
  
  
  
  pars <- expand.grid(eta=eta_seq(model$B),
                      delta=delta_seq(model$D))
  
  
  
  out <- do.call(rbind,lapply(maxchrom, function(k){
    if(D=="log(i*mu)^2/K^2"&k==4){
      pars <- expand.grid(eta=0.0001,
                          delta=c(0.013013,0.013014,0.013015,0.01304,0.0132,0.014,0.015,seq(0.02,1,0.01)))
    }
    if(D=="log(i*mu)^2/K^2"&k==8){
      pars <- expand.grid(eta=0.0001,
                          delta=c(0.0000658,0.0001,0.001,0.01,0.1))
    }
    
    if(D=="log(i*mu)^2/K^2"&k==16){
      pars <- expand.grid(eta=0.0001,
                          delta=c(c(1.022,1.024,1.026,1.028,1.03)*10^-8,10^(-7:0)))
    }
    
    if(D=="log(i*mu)^2/K^2"&k==32){
      pars <- expand.grid(eta=0.0001,
                          delta=c(c(5.71,5.72,5.74,5.78,6)*10^-16,10^(-15:0)))
    }
    out <- do.call(rbind,lapply(1:nrow(pars),function(i){
      wrap_run_model(pars$eta[i],pars$delta[i],maxchrom=k,model)
    }))
    out <- process_cloud(out,B,D)
  }))
  
  out$B <- B
  out$D <- D
  return(out)
  
}

eta_seq <- function(B){
  id <- B(NaN,NaN,label=TRUE)
  if(id==" B(i)==beta") return(c(0.0001,seq(0.1,0.9,0.1),0.9999))
  if(id=="B(i)==1/(beta+e^sqrt( abs( i-3 ) ))") return(c(0.001,0.01,0.1,0.5,2.5,10,40,100,250))
  if(id=="B(i)==e^(-beta*i)") return(c(0.0001,seq(0.1,0.9,0.1),0.9999))
  
  
}

delta_seq <- function(D){
  id <- D(NaN,NaN,label=TRUE)
  if(id==" D(i)==mu") return(c(0.001,seq(0.1,0.9,0.1),0.99,0.999,0.9999))
  if(id=="D(i)==1-1/(mu+e^sqrt( abs( i-3 ) ))") return(c(0.0000001,0.001,0.01,0.1,0.5,2.5,10,40,100,250,1000))
  if(id=="log(i*mu)^2/K^2") return(c(0.0000001,0.001,0.01,0.1,0.5,2.5,10,40,100,250,1000))
  
}

run_model_to_ss <- function(eta,delta,B,D,maxchrom){
  
  B <- list(B0,B1,B2)[[which(c(B0(NaN,NaN,TRUE),B1(NaN,NaN,TRUE),B2(NaN,NaN,TRUE))==B)]]
  D <- list(D0,D1,D2)[[which(c(D0(NaN,NaN,TRUE),D1(NaN,NaN,TRUE),D2(NaN,NaN,TRUE))==D)]]
  
  #model <- 2
  lambda <- 1
  state <- rep(0,maxchrom)
  state[2] <- 1
  
  q <- generate_q(eta,delta,B,D,maxchrom)
  
  parms <- list(q=q)
  steady_r <- FALSE
  t_end <- 100
  
  while(!is.na(steady_r)&!steady_r&t_end<1e6){
    out <- ode.1D(y=state,times=seq(0,t_end,t_end/1000),func=chrmod,parms=parms,nspec=1)
    tt <- out[,1]
    out <- out[,2:ncol(out)]
    out <- out/rowSums(out)
    d_out <- apply(out,2,diff)
    d_out <- sapply(1:nrow(d_out), function(i) max(abs(d_out[i,]/out[i,]))<0.001)
    steady_r <- as.logical(max(d_out))
    t_end <- t_end*10
  }
  
  
  tt <- head(tt[d_out],1)
  y <- head(out[d_out,],1)
  if(is.na(steady_r)) steady_r <- FALSE
  
  if(!steady_r){
    y <- NaN
    tt <- NaN
  }
  
  
  ploidy <- sum(y*(1:length(y)))
  data.frame(ploidy=ploidy,time=tt)
}

run_model_check_growth <- function(eta,delta,B,D,maxchrom,t_end=100){
  #model <- 2
  lambda <- 1
  state <- rep(0,maxchrom)
  state[2] <- 1
  
  q <- generate_q(eta,delta,B,D,maxchrom)
  
  parms <- list(q=q)
  
  out <- ode.1D(y=state,times=seq(0,t_end,t_end/1000),func=chrmod,parms=parms,nspec=1)
  tt <- out[,1]
  out <- out[,2:ncol(out)]
  n_out <- rowSums(out)
  out <- out/n_out
  d_out <- apply(out,2,diff)
  d_out <- sapply(1:nrow(d_out), function(i) max(abs(d_out[i,]/out[i,]))<0.001)
  steady_r <- as.logical(max(d_out))
  
  dNdt <- tail(diff(n_out),1)/n_out[length(n_out)-1]/(t_end/1000)
  N_end <- tail(n_out,1)
  grew <- dNdt>0
  
  data.frame(eta=eta,delta=delta,dNdt=dNdt,N_end=N_end,steady_r=steady_r,grew=grew)
}



B0 <- function(eta,i,label=FALSE,k=8) {
  if(label==TRUE){
    return(" B(i)==beta")
  }
  return(eta)
}
D0 <- function(delta,i,label=FALSE,k=8){
  if(label==TRUE){
    return(" D(i)==mu")
  }
  return(delta)
}

B1 <- function(eta,i,label=FALSE,k=8) {
  if(label==TRUE){
    return("B(i)==1/(beta+e^sqrt( abs( i-3 ) ))")
  }
  return(1/(eta+exp(sqrt(abs(i-3)))))
}

B2 <- function(eta,i,label=FALSE,k=8) {
  if(label==TRUE){
    return("B(i)==e^(-beta*i)")
  }
  return(exp(-eta*i))
}

D1 <- function(delta,i,label=FALSE,k=8){
  if(label==TRUE){
    return("D(i)==1-1/(mu+e^sqrt( abs( i-3 ) ))")
  }
  return(1-1/(delta+exp(sqrt(abs(i-3)))) )
} 

D2 <- function(delta,i,label=FALSE,k=8){
  if(label==TRUE){
    return("log(i*mu)^2/K^2")
  }
  return(log(i*delta)^2/k^2 )
} 