# Title        : MLE method for SNPs data
# Objective    : Contains implementation of the model (EM-algorithm) and supporting functions
# Created by   : Christian Tsoungui Obama, Kristan. A. Schneider
# Created on   : 05.05.22
# Last modified: 05.05.22

#################################
# Function varsets(n,l) outputs all possible vectors of length n with entries 0,.., l-1
# in an l^n x n matrix
#################################
varsets <- function(l,n) {
  # l = 2^2-1 (case of biallelic loci) # n = number of heterozygote loci
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

#################################
# Function hapl(arch) takes as input a vector of the number of alleles per locu and outputs
#  a nm x 2 matrix of all possible haplotypes in geadic representation
#################################
hapl <- function(arch){
  H <- array(0,c(prod(arch),2))
  H[,1] <- rep(0:(arch[1]-1), each=arch[2])
  H[,2] <- rep(0:(arch[2]-1), arch[1])
  H
}

#################################
# Function cpoiss(lambda,n) outputs n randomly drawn integer from a condidtional Poisson distribution 
# with parameter lambda
#################################
cpoiss <- function(lambda,n){
  m <- 100 # to accelerate computation it is assumed that m<100 is generically drawn
  out <- c()
  x <- runif(n,min=0,max=1)
  p0 <- ppois(0,lambda)
  nc <- 1/(1-exp(-lambda))
  pvec <- (ppois(1:m,lambda)-p0)*nc
  pvec <- c(pvec,1) 
  for (i in 1:n){
    k <- 1
    while(x[i] > pvec[k]){
      k <- k+1
    }
    if(k >= m){ # if a m>=100 is drawn this is executed
      k <- k+1
      a <- dpois(k,lambda)*nc
      b <- pvec[m]+a
      while(x[i]>b){
        k <- k+1
        a <- a*lambda/k
        b <- b+a
      }
    }
    out <- c(out, k) 
  }
  out
}

#################################
# The function obs(M) gives a representation for an infection with haplotypes given by the matrix M. 
# The input M is  k x n matrix with entries 0 and 1, where each row  is a haplotype corresponding to 
# a 0-1 vector. Obs returns the corresponding vector representation of the observ
#################################
obs <- function(M, arch){
  M <- matrix(M, ncol = 2)
  out <- array(0,2)
  arch1 <- arch - 1
  x <- list()
  x[[1]] <- seq(0,arch1[1])
  x[[2]] <- seq(0,arch1[2])

  for(i in 1:2){
    bin <- 2^x[[i]]
    binx <- is.element(x[[i]], M[,i]) 
    out[i] <- bin %*% binx - 1
  }
  out
}

#################################
# Function datasetgen(P,lambda,N,n) is used to simulate data. It generates N observations assuming n biallelic loci with haplotype distribution P 
# wich must be a vector of length 2n a N x n matrix of observations sampled using the multinomial and Poisson distribution 
# of parameters (m, P) and lambda respectively. m is the MOI for the corresponding sample.
#################################
datasetgen <- function(P,lambda,N,arch){ 
  H <- hapl(arch)                    # Set of possible haplotypes
  out <- matrix(0,nrow=N, ncol=2)
  m <- cpoiss(lambda,N)              # MOI values for each sample following CPoiss(lambda)
  for(j in 1:N){                     
    s <- rmultinom(1, m[j], P) #multinomially select M[j] haplotypes from the haplotype pool
    out[j,] <- obs(H[s!=0,], arch)-1 #Summing up the trianary representation of a number representing the infection
  } #vector of infections
  out
}

datagen <- function(P,lambda,N,arch){ 
  n <- 2
  H <- hapl(arch)       # Set of possible haplotypes
  vec <- rev((2^arch-1)^c(0,1)) #c((2^arch[2]-1),1)    # Vector for geadic representation
  out <- array(0,N)
  m <- cpoiss(lambda,N) # MOI values for each sample following CPoiss(lambda)
  for(j in 1:N){        # for each sample
    s <- rmultinom(1, m[j], P) #multinomially select M[j] haplotypes from the haplotype pool
    out[j] <- vec%*%obs(H[s!=0,], arch)+1 #Summing up the trianary representation of a number representing the infection
  } #vector of infections
  out <- t(as.data.frame(summary.factor(out))) #vector of how many times each infection that is effectively present appears in the dataset
  vals <- as.integer(colnames(out))-1 #Infections present in the dataset
  dat <- array(0,c(length(vals),n))
  for(k in 1:n){ #for each locus
    re <- vals%%vec[k]
    dat[,k] <- (vals-re)/vec[k]
    vals <- re
  } #Trianary representation of each infection present in the dataset
  list(dat, c(out))
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
# Function gen_func(x,lambd) calculates the value of the generating function of the conditional Poisson distribution for x.
#################################
gen_func <- function(x, lambd){
  (exp(x*lambd)-1)/(exp(lambd) - 1)
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
  cond1 <- 1 
  eps <- 10^-8
  while(cond1>eps){
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
    while(cond2 > eps){
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
  list(pp, la)
}

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

#################################
# The function reform(X1,id) takes as input the dataset in the 0-1-2-notation and returns a matrix of the observations,
# and a vector of the counts counts of those observations, i.e., number of times each observation is made in the dataset.
#################################
reform <- function(X1, id = TRUE){
    # This function formats the data for the MLE function
    # Remove the id column
    if(id){
        X1 <- X1[,-1]
    }
    
    # Deriving the number of time Nx each observation is made in the dataset
    Nx <- sampl(X1)

    # Matrix of  observed observations
    nloci <- ncol(X1)
    nHapl <- 2^nloci
    H <- hapl(nloci)
    trin <- 3^((nloci-1):0)

    X1 <- as.matrix(X1, ncol=nloci)%*%trin + 1
    X1 <- t(as.data.frame(summary.factor(X1)))    # Observations present in the dataset
    vals <- as.integer(colnames(X1))-1           
    dat <- array(0,c(length(vals),nloci))
    for(k in 0:(nloci-1)){ 
        re <- vals%%(3^(nloci-k-1))
        dat[,nloci-k] <- (vals-re)/(3^(nloci-k-1))  
        vals <- re
    }
    list(dat, Nx)
}

#################################
# The function mle(df,id) wraps the reform(X1,id) and estsnpmodel(X, Nx) to find the MLEs and outputs the estimates
# for haplotype frequencies, Poisson parameters, and a matrix of detected haplotypes.
#################################
mle <- function(df, id = TRUE){
    # This function removes the ID column if there is one,
    # then it derives the number of time each observation is made in the dataset,
    # finally, the MLE are obtained and return in a list.

    dat1 <- reform(df, id=TRUE)
    X <- dat1[[1]]
    Nx <- dat1[[2]]
    nloci <- ncol(X)

    # MLEs
    out <- estsnpmodel(X, Nx)
    out2 <- out[[2]]
    rnames <- as.integer(rownames(out2)) - 1
    nh <- length(rnames)
    dat <- array(0,c(nh,nloci))
    for(k in 0:(nloci-1)){ #for each locus
        re <- rnames%%(2^(nloci-k-1))
        dat[,nloci-k] <- (rnames-re)/(2^(nloci-k-1))
        rnames <- re
    } 
    for(i in 1:nh){
        rnames[i] <- paste(dat[i,], collapse = '')
    }
    rownames(out2) <- rnames
    out <- list(unlist(out[1]), t(out2), dat)
    names(out) <- c(expression(lambda), 'p', 'haplotypes')
    out
}

#################################
# The function adhocmodel(X,Nx) calculates the relative prevalence of the haplotypes
# conditionned on unambiguous observations.
#################################
adhocmodel <- function(X, Nx){

  X <- cbind(X, Nx)
  n <- ncol(X)
  nloci <- n-1
  # estimate haplotype frequencies
  nhpl <- 2^nloci
  
  # extract unambiguous observations
  tmp1 <- matrix(X[,1:nloci]==2, ncol = nloci)
  X1 <- X[rowSums(tmp1)<2,]

  if(!all(is.na(X1))){  # if there are unambiguous infections
    n1 <- nrow(X1)
    if(is.null(n1)){ # if there is only one unambiguous infection
      n1 <- 1
    }
    X <- matrix(X1, nrow = n1)
    # find indexes of multiple infections
    tmp2 <- matrix(X[,1:nloci]==2, ncol = nloci)
    idx1 <- which(rowSums(tmp2)==1)

    if(length(idx1)>0){
      # single infections
      s <- X[-idx1,]
    }else {
      s <- X
    }
    
    # find all the haplotypes in X
    for(i in idx1){
      y <- X[i,]
      idx2 <- which(y[1:nloci]==2)
      h <- array(rep(y,2), c(n, 2))
      h[idx2,] <- c(0,1)
      # add haplotypes in s
      s <- rbind(s,t(h))
    }
    
    # binary representation
    bin <- 2^((nloci-1):0)
    pp  <- s[,1:nloci]%*%bin+1
    pp  <- cbind(pp,s[,n])

    # observed haplotypes
    idx3 <- as.integer(colnames(t(as.data.frame(summary.factor(pp[,1])))))
    
    # Frequencies estimates
    p <- matrix(0, ncol=length(idx3))
    tot <- sum(pp[,2])
    for (i in idx3){
      idx4 <- which(pp[,1]==i)
      p[,which(i==idx3)] <- sum(pp[idx4,2])/tot
    }
      colnames(p) <-  idx3
  }else{
    # Frequencies estimates
    p <- matrix(0, ncol=nhpl)
  }
  p
}

#################################
# The function sampl(dat) finds the number of occurences of each observation in the dataset dat.
# The output is a vector of those numbers.
#################################
sampl <- function(dat){
    nloci <- ncol(dat)
    trin <- 3^((nloci-1):0)
    out <- table(as.matrix(dat, ncol=nloci)%*%trin + 1)
    out
}

#################################
# The function estunobsprev(estim) calculates the unobservable prevalence of the haplotypes.
# The input is the MLEs obtained by the function estsnpmodel(X, Nx).
#################################
estunobsprev <- function(estim){
  # This function estimates the unobservable prevalence as defined in the manuscript of tsoungui et.al, titled
  # "A maximum-likelihood method to estimate haplotype frequencies and prevalence alongside multiplicity of infection from SNPs data"

    ## For each set of estimates, compute prevalence
    prev <- (exp(estim[[1]]) - exp(1-estim[[2]])^estim[[1]])/(exp(estim[[1]])-1)
    
    prev
}

#################################
# The function estcondprev(estim) calculates the prevalence of the haplotypes conditioned on.
# unambiguous observations. The input is the MLEs obtained by the function estsnpmodel(X, Nx).
#################################
estcondprev <- function(estim){
  # This function estimates the unambiguous prevalence as defined in the manuscript of tsoungui et.al, titled
  # "A maximum-likelihood method to estimate haplotype frequencies and prevalence alongside multiplicity of 
  # infection from SNPs data"

    # Number of loci
    numb_Loci <- ncol(estim[[3]])

    # Table of all possible haplotypes
    Hapl <- hapl(numb_Loci)

    ## For each haplotype in the table, build the set of observation Uh
    numb_Hapl_Uh <- numb_Loci + 1

    cnames <- colnames(estim[[2]])

    ## Access the estimates
    tmp2 <- estim[[2]]
    pickhap <- estim[[3]]%*%2^(0:(numb_Loci-1))+1

    nHapl <- length(pickhap)

    prev <- matrix(0, ncol = nHapl)
    colnames(prev) <- cnames
    numh <- rep(0, nHapl)
    denh <- rep(0, nHapl)

    # Find ambiguous prevalence for each observed haplotype
    for (idx in pickhap[,1]){ 
        # Build uh
        uh                  <- t(array(rep(Hapl[idx,], numb_Loci), dim=c(numb_Loci, numb_Hapl_Uh)))
        uh[2:numb_Hapl_Uh,] <- (uh[2:numb_Hapl_Uh,]+diag(numb_Loci))%%2

        # Indexes of haplotypes forming unambiguous observations with h (Uh)
        pick1 <- uh%*%2^((numb_Loci-1):0)+1

        ## Pick the right frequencies estimates 
        pickh <- which(pickhap == pick1[1])
        GPh   <- gen_func(tmp2[pickh], estim[[1]])

        # Picking haplotypes in Uh with non zero frequencies
        pick2 <- which(pickhap%in%pick1) # indices of observed haplotypes in uh
        tmp3 <- rep(0, numb_Hapl_Uh)
        tmp3[which(pick1%in%pickhap)] <- tmp2[pick2]
        rem <- which(tmp3==tmp2[pickh])[1]
        tmp3 <- tmp3[-rem]
        i <- which(idx==pickhap)
        for(j in 1:(numb_Hapl_Uh-1)){ 
            GPartFreq  <- gen_func(tmp3[j], estim[[1]])
            GFreq      <- gen_func(sum(c(tmp2[pickh],tmp3[j])), estim[[1]])
            tmp        <- GFreq - GPartFreq
            numh[i]    <- numh[i] + tmp
            denh[i]    <- denh[i] + tmp/2
        }
        numh[i] <- numh[i] - (numb_Loci - 1)*GPh
        denh[i] <- denh[i] - (numb_Loci/2 - 1)*GPh
    }
    for(q in 1:nHapl){
        prev[,q] <- numh[q]/sum(denh)
    }
    colnames(prev) <- cnames
    prev
}

#################################
# The function estrelprev(df,id) wraps the functions reform(df,id) and adhocmodel(X, Nx) to calculate 
#  the relative prevalence of the haplotypes conditionned on unambiguous observations.
#################################
estrelprev <- function(df, id = TRUE){
    # This function removes the ID column if there is one,
    # then it derives the number of time each observation is made in the dataset,
    # finally, the MLE are obtained and return in a list.

    dat1 <- reform(df, id=TRUE)
    X <- dat1[[1]]
    Nx <- dat1[[2]]
    nloci <- ncol(X)

    # MLEs
    out <- adhocmodel(X, Nx)
    cnames <- as.integer(colnames(out)) - 1
    nh <- length(cnames)
    dat <- array(0,c(nh,nloci))
    for(k in 0:(nloci-1)){ #for each locus
        re <- cnames%%(2^(nloci-k-1))
        dat[,nloci-k] <- (cnames-re)/(2^(nloci-k-1))
        cnames <- re
    } 
    for(i in 1:nh){
        cnames[i] <- paste(dat[i,], collapse = '')
    }
    colnames(out) <- cnames
    out
}