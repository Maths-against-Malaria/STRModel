# Title        : MLE method for SNPs data
# Objective    : Contains implementation of the model (EM-algorithm) and supporting functions
# Created by   : Christian Tsoungui Obama, Kristan. A. Schneider
# Created on   : 05.05.22
# Last modified: 05.05.22

#################################
# Function varsets(n,l) outputs all possible vectors of length n with entries 0,.., l-1
# in an l^n x n matrix
#################################
varsets1 <- function(l){   #calculate all var sets
  # n number of loci
  # l number of alleles per locus
  n <- length(l)
  B <- array(0,c(prod(l),n))
  B[1:l[1],1] <- 0:(l[1]-1)
  lkmo <- l[1]
  if(n>1){
    for(k in 2:n){
      if(l[k]>1){
        lk <- lkmo*l[k]
        pick1 <- (lkmo+1):lk
        B[pick1,] <- B[rep(1:lkmo,l[k]-1),]
        B[pick1,k] <- rep(1:(l[k]-1),each=lkmo)
        lkmo <- lk   
      }
    }
  }
  B
}

varsets2 <- function(l){   #calculate all var sets
  # n number of loci
  # l number of alleles per locus
  n <- length(l)
  B <- array(1,c(prod(l),n))
  B[1:l[1],1] <- 1:l[1]
  lkmo <- l[1]
  if(n>1){
    for(k in 2:n){
      if(l[k]>1){
        lk <- lkmo*l[k]
        pick1 <- (lkmo+1):lk
        B[pick1,] <- B[rep(1:lkmo,l[k]-1),]
        B[pick1,k] <- rep(2:l[k],each=lkmo)
        lkmo <- lk
      }
    }
  }
  B
}

gead <- function(x,l,n){   ## calculates geadic expression of each element of vectorx 
  l <- rep(l,n)
  out <- array(0,c(length(x),n))
  div <- c(1,cumprod(l[1:(n-1)]))
  for(k in n:1){
    r <- x%%div[k]
    out[,k] <- (x-r)/div[k]
    x <- r
  }
  out
}

#################################
# The function estsnpmodel(X,Nx) implements the EM algorithm and returns the MLEs, i.e., 
# estimates of haplotype frequencies and Poisson parameter.
#################################
eststrmodel <- function(dat, arch){
  X <- dat[[1]] + 1
  Nx <- dat[[2]]
  N <- sum(Nx)
  nn <- nrow(X)
  l <- arch
  n <- length(l)
  x <- X
  hapll <- c()

  ggead <- c(l[2], 1)
  Hx <- list() 
  Ax <- list() 
  alnum <- list()
  bin2num <- list()
  for(k in 1:n){
    alnum[[k]] <- 1:l[k]
    bin2num[[k]] <- 2^(0:(l[k]-1))
  }
  alcnt <- array(,n)
  for(u in 1:nrow(x)){
    Hx[[u]] <- list()
    Ax[[u]] <- list()
    Hx[[u]][[1]] <- array(0,n)
    Hx[[u]][[2]] <- list()
    Hx[[u]][[3]] <- list()
    Hx[[u]][[4]] <- list()
    Hx[[u]][[5]] <- list()
    Hx[[u]][[6]] <- array(,n)
    Ax[[u]][[1]] <- list()
    Ax[[u]][[2]] <- list()
    Ax[[u]][[3]] <- list()
    Ax[[u]][[4]] <- list()
    for(k in 1:n){
      temp <- gead(x[u,k],2,l[[k]])
      temp1 <- varsets1(temp+1)[-1,]
      Hx[[u]][[1]][k] <- sum(temp)
      Hx[[u]][[2]][[k]] <- (alnum[[k]])[temp*alnum[[k]]]
      Hx[[u]][[3]][[k]] <- temp1
      Hx[[u]][[4]][[k]] <- temp1%*%bin2num[[k]] 
      Hx[[u]][[5]][[k]] <- Hx[[u]][[3]][[k]]%*%rep(1,l[k])
      Hx[[u]][[6]][k] <- length(temp1%*%bin2num[[k]])
    }
    vz1 <- sum(Hx[[u]][[1]])
    temp2 <- prod(Hx[[u]][[6]])
    Ax[[u]][[1]] <- temp2
    Ax[[u]][[2]] <- varsets2(Hx[[u]][[6]])
    for(k in 1:n){
      Ax[[u]][[2]][,k] <- Hx[[u]][[4]][[k]][Ax[[u]][[2]][,k]]
    }
    for(j in 1:temp2){
      Ax[[u]][[3]][[j]] <- list() 
      for(k in 1:n){
        temp <- gead(Ax[[u]][[2]][j,k],2,l[[k]])
        temp1 <- (alnum[[k]])[temp*alnum[[k]]]
        alcnt[k] <- length(temp1)
        Ax[[u]][[3]][[j]][[k]] <- temp1
      }
      Ax[[u]][[4]][[j]] <- varsets2(alcnt)
      for(k in 1:n){
        Ax[[u]][[4]][[j]][,k] <- Ax[[u]][[3]][[j]][[k]][Ax[[u]][[4]][[j]][,k]]
      }
      Ax[[u]][[4]][[j]] <- as.character((Ax[[u]][[4]][[j]]-1)%*%ggead+1)
      Ax[[u]][[3]][[j]] <- (-1)^(vz1+sum(alcnt))
    }
    hapll[[u]] <- Ax[[u]][[4]][[temp2]]
  }  

  hapl1 <- unique(unlist(hapll))
  H <- length(hapl1)
  pp <- array(rep(1/H,H),c(H,1))
  rownames(pp) <- hapl1

  num0 <- pp*0
  num <- num0
  rownames(num) <- hapl1
  Bcoeff <- num0
  rownames(Bcoeff) <- hapl1
  la <- 2
  eps <- 10^-8
  cond1 <- 1 
  t <- 0
  while(cond1>eps && t<500){
    Ccoeff <- 0
    Bcoeff <- num0    # reset B coefficients to 0 in next iteration
    num <- num0       # reset numerator to 0 in next iteration
    for(u in 1:nn){
      denom <- 0
      num <- num0
      CC <- 0
      for(k in 1:Ax[[u]][[1]]){
        p <- sum(pp[Ax[[u]][[4]][[k]],])
        vz <- Ax[[u]][[3]][[k]]
        lap <- la*p
        exlap <- vz*exp(lap)
        denom <- denom + exlap-vz
        num[Ax[[u]][[4]][[k]],] <- num[Ax[[u]][[4]][[k]],] + exlap
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
    
    ### Newton step
    cond2 <- 1
    xt <- Ccoeff   ### good initial condition
    tau <- 0
    while(cond2 > eps  &&  tau<300){
      ex <- exp(-xt)
      xtn <- xt + (1-ex)*(xt + Ccoeff*ex - Ccoeff)/(ex*xt+ex-1)
      if(is.nan(xtn) || (tau == 299) || xtn < 0){
          xtn <- runif(1, 0.1, 2.5)
      }
      cond2 <- abs(xtn-xt)
      xt <- xtn
    }
    cond1 <- abs(xt-la)+sqrt(sum((pp-ppn)^2))
    la <- xt
    pp <- ppn
  }
  list(la, pp)
}

strmodel <- function(dat, arch){
  X <- dat[[1]] + 1
  Nx <- dat[[2]]
  N <- sum(Nx)
  nn <- nrow(X)
  l <- arch
  n <- length(l)
  x <- X
  hapll <- c()

  ggead <- c(l[2], 1)
  Hx <- list() 
  Ax <- list() 
  alnum <- list()
  bin2num <- list()
  for(k in 1:n){
    alnum[[k]] <- 1:l[k]
    bin2num[[k]] <- 2^(0:(l[k]-1))
  }
  alcnt <- array(,n)
  for(u in 1:nrow(x)){
    Hx[[u]] <- list()
    Ax[[u]] <- list()
    Hx[[u]][[1]] <- array(0,n)
    Hx[[u]][[2]] <- list()
    Hx[[u]][[3]] <- list()
    Hx[[u]][[4]] <- list()
    Hx[[u]][[5]] <- list()
    Hx[[u]][[6]] <- array(,n)
    Ax[[u]][[1]] <- list()
    Ax[[u]][[2]] <- list()
    Ax[[u]][[3]] <- list()
    Ax[[u]][[4]] <- list()
    for(k in 1:n){
      temp <- gead(x[u,k],2,l[[k]])
      temp1 <- varsets1(temp+1)[-1,]
      Hx[[u]][[1]][k] <- sum(temp)
      Hx[[u]][[2]][[k]] <- (alnum[[k]])[temp*alnum[[k]]]
      Hx[[u]][[3]][[k]] <- temp1
      Hx[[u]][[4]][[k]] <- temp1%*%bin2num[[k]] 
      Hx[[u]][[5]][[k]] <- Hx[[u]][[3]][[k]]%*%rep(1,l[k])
      Hx[[u]][[6]][k] <- length(temp1%*%bin2num[[k]])
    }
    vz1 <- sum(Hx[[u]][[1]])
    temp2 <- prod(Hx[[u]][[6]])
    Ax[[u]][[1]] <- temp2
    Ax[[u]][[2]] <- varsets2(Hx[[u]][[6]])
    for(k in 1:n){
      Ax[[u]][[2]][,k] <- Hx[[u]][[4]][[k]][Ax[[u]][[2]][,k]]
    }
    for(j in 1:temp2){
      Ax[[u]][[3]][[j]] <- list() 
      for(k in 1:n){
        temp <- gead(Ax[[u]][[2]][j,k],2,l[[k]])
        temp1 <- (alnum[[k]])[temp*alnum[[k]]]
        alcnt[k] <- length(temp1)
        Ax[[u]][[3]][[j]][[k]] <- temp1
      }
      Ax[[u]][[4]][[j]] <- varsets2(alcnt)
      for(k in 1:n){
        Ax[[u]][[4]][[j]][,k] <- Ax[[u]][[3]][[j]][[k]][Ax[[u]][[4]][[j]][,k]]
      }
      Ax[[u]][[4]][[j]] <- as.character((Ax[[u]][[4]][[j]]-1)%*%ggead+1)
      Ax[[u]][[3]][[j]] <- (-1)^(vz1+sum(alcnt))
    }
    hapll[[u]] <- Ax[[u]][[4]][[temp2]]
  }  

  hapl1 <- unique(unlist(hapll))
  H <- length(hapl1)
  pp <- array(rep(1/H,H),c(H,1))
  rownames(pp) <- hapl1

  num0 <- pp*0
  num <- num0
  rownames(num) <- hapl1
  Bcoeff <- num0
  rownames(Bcoeff) <- hapl1
  la <- 2
  eps <- 10^-8
  cond1 <- 1 
  t <- 0
  while(cond1>eps && t<500){
    Ccoeff <- 0
    Bcoeff <- num0    # reset B coefficients to 0 in next iteration
    num <- num0       # reset numerator to 0 in next iteration
    for(u in 1:nn){
      denom <- 0
      num <- num0
      CC <- 0
      for(k in 1:Ax[[u]][[1]]){
        p <- sum(pp[Ax[[u]][[4]][[k]],])
        vz <- Ax[[u]][[3]][[k]]
        lap <- la*p
        exlap <- vz*exp(lap)
        denom <- denom + exlap-vz
        num[Ax[[u]][[4]][[k]],] <- num[Ax[[u]][[4]][[k]],] + exlap
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
    
    ### Newton step
    cond2 <- 1
    xt <- Ccoeff   ### good initial condition
    tau <- 0
    while(cond2 > eps  &&  tau<300){
      ex <- exp(-xt)
      xtn <- xt + (1-ex)*(xt + Ccoeff*ex - Ccoeff)/(ex*xt+ex-1)
      if(is.nan(xtn) || (tau == 299) || xtn < 0){
          xtn <- runif(1, 0.1, 2.5)
      }
      cond2 <- abs(xtn-xt)
      xt <- xtn
    }
    cond1 <- abs(xt-la)+sqrt(sum((pp-ppn)^2))
    la <- xt
    pp <- ppn
  }
  ## Ordering the frequencies
  pp <- pp[order(as.numeric(rownames(pp))), ]

  ## Setting the frequencies of the unobserved haplotypes to 0.0
  nhapl <- prod(arch)

  if(length(pp)<nhapl){
    out <- t(pp)
    name <- colnames(out)
    cnt <- 0
    for (i in 1:nhapl) {
      if (is.element(as.character(i), name)){
        cnt <- cnt + 1
      }else{
        pp <- append(pp, list(x = 0.0), i-1)
      }
    }
  }
  list(la, pp)
}

dat <- readRDS('/Volumes/GoogleDrive-117934057836063832284/My Drive/Maths against Malaria/Christian/Models/STRModel/dataset/data.rds')
