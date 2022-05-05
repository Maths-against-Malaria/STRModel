
#---------------------------------------

varsets <- function(l,n){   #calculate all var sets #n = number of loci #l = number of observations per locus
  B <- array(0,c(l^n,n))
  B[1:l,1] <- 0:(l-1)
  lkmo <- l
  if(n>1){
    for(k in 2:n){
      lk <- lkmo*l
      pick1 <- (lkmo+1):lk
      B[pick1,] <- B[rep(1:lkmo,l-1),]
      B[pick1,k] <- rep(1:(l-1),each=lkmo)
      lkmo <- lk                        
    }
  }
  B
}

nbialModel <- function(Nx, X){
  
  eps <- 10^-8 #Error
  N <- sum(Nx) #Sample size
  nn <- nrow(X)  #Number of different observations present in the dataset
  Ax <- list()
  n <- ncol(X)  #Number of loci
  
  for(u in 1:nn){
    xx <- array(X[u,],c(1,n)) #xx = observation
    sel <- (1:n)[xx==2] #Identifying the loci where the 2 alleles are observed
    l <- length(sel)    #Counting the number of loci where the 2 alleles are observed
    
    if(l==0){   # If the infection is a haplotype (only one allele per locus)
      yy <- xx
    }else{ 
      yy <- xx[rep(1,3^l),]
      yy[,sel] <- varsets(3,l) # Set of all possible observations which combinations can form xx $\mathscr{A}_{y}$
    }
    bin <- 2^(0:(n-1))
    iilist <- list()
    siglist <- list()
    for(i in 1:3^l){
      y1 <- array(yy[i,],c(1,n)) # Observation {\pmb y} in the set $\mathscr{A}_{y}$
      sel <- (1:n)[y1==2]
      l1 <- length(sel)
      if(l1==0){
        ii <- y1
      }else{
        ii <- y1[rep(1,2^l1),]
        ii[,sel] <- varsets(2,l1)
      }
      iilist[[i]] <- as.character(ii%*%bin+1)
      siglist[[i]] <- (-1)^(l-l1)
      
    }
    Ax[[u]] <- list(iilist,siglist,3^l)
  }
  
  # list of all occuring halotypes 
  hapl1 <- c()
  for(u in 1:nn){
    hapl1 <- c(unlist(Ax[[u]][[1]]),hapl1)
  }
  hapl1 <- unique(hapl1)
  H <- length(hapl1)
  
  pp <- array(rep(1/H,H),c(H,1))   #Initial frequency distribution for the EM algorithm
  rownames(pp) <- hapl1
  
  la <- 2                         #Initial value of lambda for the EM algo.
  
  num0 <- pp*0
  cond1 <- 1  ## Initializing the condition to stop EM alg! 
  Bcoeff <- num0
  num <- num0
  rownames(num) <- hapl1
  rownames(Bcoeff) <- hapl1
  
  while(cond1>eps){
    Ccoeff <- 0
    Bcoeff <- num0 #reset B coefficients to 0 in next iteration
    num <- num0  #reset numerator to 0 in next iteration
    for(u in 1:nn){
      denom <- 0
      num <- num0
      CC <- 0
      for(k in 1:Ax[[u]][[3]]){
        p <- sum(pp[Ax[[u]][[1]][[k]],])
        vz <- Ax[[u]][[2]][[k]]
        lap <- la*p
        exlap <- vz*exp(lap)
        denom <- denom + exlap-vz 
        num[Ax[[u]][[1]][[k]],] <- num[Ax[[u]][[1]][[k]],]+ exlap
        CC <- CC + exlap*p
      }
      num <- num*pp
      denom <- Nx[u]/denom
      denom <- la*denom
      Ccoeff <- Ccoeff + CC*denom
      
      Bcoeff <- Bcoeff + num*denom
      
    }
    
    Ccoeff <- Ccoeff/N
    ppn <- Bcoeff/(sum(Bcoeff))
    
    ##*************Newton-Raphson to estimate the lambda parameter
    cond2 <- 1
    xt <- Ccoeff   
    
    while(cond2 > eps){
      ex <- exp(-xt)
      xtn <- xt + (1-ex)*(xt + Ccoeff*ex - Ccoeff)/(ex*xt+ex-1)
      if(is.nan(xtn)){ #Replacing NA value by random values of lambda
        xtn <- sample(seq(0.1,3,0.1),1)
      }
      cond2 <- abs(xtn-xt)
      xt <- xtn
    }
    # Filling the frequency values = 0 in the frequency vector
    cond1 <- abs(xt-la)+sqrt(sum((pp-ppn)^2))
    la <- xt
    pp <- ppn
  }
  nhapl <- 2^n
  
  if(length(pp)<nhapl){
    out <- t(pp)
    name <- colnames(out)
    cnt <- 0
    
    for (i in 1:nhapl) {
      if (is.element(as.character(i), name)){
        cnt <- cnt + 1
      }else{
        pp <- append(pp, list(x = 0), i-1)
      }
    }
  }
  
  return(c(la, pp))
}