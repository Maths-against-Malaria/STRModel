# Title        : MLE method for SNPs data
# Objective    : Contains implementation of the model (EM-algorithm) and supporting functions
# Created by   : Christian Tsoungui Obama, Kristan. A. Schneider
# Created on   : 05.05.22
# Last modified: 26.10.22

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

setUh <- function(haplo, cardUh, arch){
  uh  <- t(array(rep(haplo, cardUh), dim = c(2, cardUh)))
  idx <- which(arch==max(arch))
  if(length(idx) > 1){
    idx <- idx[1]
  }
  v1 <- 0:(arch[idx]-1)
  v  <- v1[-which(v1==uh[1,idx])]
  uh[2:(length(v)+1), idx] <- v
  idx  <- which(1:2 !=idx)
  idx2 <- length(v)+2
  v2   <- 0:(arch[idx]-1)
  v    <- v2[-which(v2==uh[1,idx])]
  uh[idx2:cardUh, idx] <- v
  uh
}

#################################
# The function estsnpmodel(X,Nx) implements the EM algorithm and returns the MLEs, i.e., 
# estimates of haplotype frequencies and Poisson parameter.
#################################
strmodel0 <- function(dat, arch){
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
    t <- t+1
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

    # Replacing NaN's in Ak by 0
    cnt <- sum(is.nan(Bcoeff))
    if(cnt > 0){
        break
      }else{
        ppn <- Bcoeff/(sum(Bcoeff))
    }

    ### Newton step
    cond2 <- 1
    xt    <- Ccoeff   ### good initial condition
    tau   <- 0
    while(cond2 > eps  &&  tau<300){
      tau <- tau + 1
      ex  <- exp(-xt)
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

eststrmodel <- function(dat, arch){
  out <- strmodel0(dat,arch)
  la <- out[[1]]
  pp <- out[[2]]
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
  list(la, unlist(pp))
}

#################################
# The function strmodel_plugin(dat, arch, lam) implements the EM algorithm with the Poisson parameter as plugin estimate
# and returns the plugin Poisson parameter and the MLEs for haplotype frequencies.
#################################
strmodel_plugin <- function(dat, arch, lam){
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
  #la <- lam
  eps <- 10^-8
  cond1 <- 1 
  t <- 0
  while(cond1>eps && t<500){
    t <- t+1
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
        lap <- lam*p
        exlap <- vz*exp(lap)
        denom <- denom + exlap-vz
        num[Ax[[u]][[4]][[k]],] <- num[Ax[[u]][[4]][[k]],] + exlap
        CC <- CC + exlap*p
      }
      num <- num*pp
      denom <- Nx[u]/denom
      denom <- lam*denom
      Ccoeff <- Ccoeff + CC*denom
      Bcoeff <- Bcoeff + num*denom
    }
    
    Ccoeff <- Ccoeff/N

    # Replacing NaN's in Ak by 0
    cnt <- sum(is.nan(Bcoeff))
    if(cnt > 0){
        break
      }else{
        ppn <- Bcoeff/(sum(Bcoeff))
    }

    ### Newton step
    # cond2 <- 1
    # xt    <- Ccoeff   ### good initial condition
    # tau   <- 0
    # while(cond2 > eps  &&  tau<300){
    #   tau <- tau + 1
    #   ex  <- exp(-xt)
    #   xtn <- xt + (1-ex)*(xt + Ccoeff*ex - Ccoeff)/(ex*xt+ex-1)
    #   if(is.nan(xtn) || (tau == 299) || xtn < 0){
    #       xtn <- runif(1, 0.1, 2.5)
    #   }
    #   cond2 <- abs(xtn-xt)
    #   xt <- xtn
    # }
    # cond1 <- abs(xt-la)+sqrt(sum((pp-ppn)^2))
    # la <- xt
    cond1 <- sqrt(sum((pp-ppn)^2))
    pp <- ppn
  }
  ## Ordering the frequencies
  # pp <- pp[order(as.numeric(rownames(pp))), ]

  # ## Setting the frequencies of the unobserved haplotypes to 0.0
  # nhapl <- prod(arch)

  # if(length(pp)<nhapl){
  #   out <- t(pp)
  #   name <- colnames(out)
  #   cnt <- 0
  #   for (i in 1:nhapl) {
  #     if (is.element(as.character(i), name)){
  #       cnt <- cnt + 1
  #     }else{
  #       pp <- append(pp, list(x = 0.0), i-1)
  #     }
  #   }
  # }
  list(lam, pp)
}

#################################
# The function strmodel1(X,Nx,plugin) implements the EM algorithm with the plugin argument defining the value of the Poisson parameter 
# if plugin=NULL, the Poisson parameter is estimated from the data in which case the MLEs are obtained using the function estsnpmodel0(X, Nx),
# otherwise, the value of plugin is used as pugin-estimate of the Poisson parameter and the MLEs are obtained using the function estsnpmodel_plugin(X, Nx, plugin).
#################################
strmodel1 <- function(dat, arch, plugin=NULL){
  if(is.null(plugin)){
    out <- strmodel0(dat, arch)               # calculates the uncorrected estimate
  }else{
    out <- strmodel_plugin(dat, arch, plugin) # calculates the corrected estimate
  }
}

#################################
# The function strmodel(dat, arch) implements the EM algorithm with the option for bias correction (BC) using either
# a "bootstrap", or a "Jacknife" method, and returns the MLEs, i.e., estimates of haplotype frequencies and Poisson parameter.
#################################
strmodel <- function(dat, arch, BC=FALSE, method='bootstrap', Bbias=10000, plugin=NULL){
  out.temp <- strmodel1(dat, arch, plugin=plugin)
  rnames1 <- as.integer(rownames(out.temp[[2]])) - 1
  rnames <- rnames1
  nhap <- length(out.temp[[2]])
  X <- dat[[1]]
  Nx <- dat[[2]]
  if(BC){
    N <- sum(Nx)
    if(method == 'bootstrap'){
      infct <- vector(mode = "list", length = 2)
      prob <- Nx/N
      Estim <- array(0, dim = c((nhap+1), Bbias))
      rownames(Estim) <- c('l',(rnames1+1))
      for (l in 1:Bbias){
        samp  <- rmultinom(N, 1, prob)
        tmp   <- rowSums(samp)
        pick  <- tmp == 0
        infct[[1]]  <- X[!pick,]
        infct[[2]]  <- tmp[!pick]
        tmp1        <- strmodel1(infct, arch, plugin=plugin)
        rnames      <- as.integer(rownames(tmp1[[2]]))
        Estim[1,l]  <- unlist(tmp1[[1]])  
        Estim[as.character(rnames),l] <- unlist(tmp1[[2]])
      }
      bias  <- rowSums(Estim)/Bbias
      lamBC <- 2*out.temp[[1]][1]  - bias[1]
      ppBC  <- 2*out.temp[[2]] - bias[-1]
    }else{
      if(method=="jackknife"){
        J = length(Nx)
        Estim <- array(0, dim = c((nhap+1), J))
        rownames(Estim) <- c('l',(rnames1+1))
        for(j in 1:J){
          NxJ         <- Nx 
          NxJ[j]      <- NxJ[j]-1 
          pick        <- NxJ !=0
          infct[[1]]  <- X[pick,]
          infct[[2]]  <- NxJ[pick]
          tmp1        <- strmodel1(infct, arch, plugin=plugin)
          rnames      <- as.integer(rownames(tmp1[[2]]))
          Estim[1,j]  <- unlist(tmp1[[1]])  
          Estim[as.character(rnames),j] <- unlist(tmp1[[2]]) 
        }
        bias  <- Estim %*% Nx/N
        lamBC <- out.temp[[1]][1] - (N-1)*( bias[1] - out.temp[[1]][1]) 
        ppBC  <- out.temp[[2]] - (N-1)*( bias[-1] - out.temp[[2]])
      }else{
        warning("method needs to be either bootstrap or jackknife")
      }
    }
    out <- list(lamBC, ppBC)
  }else{
    out <- out.temp
  }
  out
}

#################################
# The function reform(X1,id) takes as input the dataset in the 0-1-2-notation and returns a matrix of the observations,
# and a vector of the counts of those observations, i.e., number of times each observation is made in the dataset.
#################################
reform <- function(DATA, arch, id = TRUE){
    # This function formats the data for the MLE function
    # Remove the id column
    X1 <- DATA
    if(id){
        X1 <- X1[,-1]
    }
    
    # Deriving the number of time Nx each observation is made in the dataset
    Nx <- sampl(X1, arch)

    # Matrix of  observed observations
    trin <- rev((2^arch-1)^c(0,1)) # geadic representation

    X1 <- as.matrix(X1, ncol=2)%*%trin + 1
    X1 <- t(as.data.frame(summary.factor(X1)))    # Observations present in the dataset
    vals <- as.numeric(colnames(X1))-1        
    dat <- array(0,c(length(vals),2))
    for(k in 1:2){ #for each locus
      re <- vals%%trin[k]
      dat[,k] <- (vals-re)/trin[k]
      vals <- re
    }
    list(dat, Nx)
}

#################################
# The function mle(df,id) wraps the reform(X1,id) and either strmodel(dat, arch) or strmodel_plugin(dat, arch, plugin) to find the MLEs 
# with or without the Poisson parameter as a plug-in estimate, respectively. Moreover, the option to ouput t bias corrected (BC) estimates with 
# confidence intervals (CI) is available. The function outputs the estimates for haplotype frequencies, Poisson parameters, and a matrix of detected haplotypes.
#################################
mle <-function(Data, arch, id=TRUE, plugin=NULL, CI=FALSE, BC=FALSE, method="bootstrap", Bbias=10000, B=10000, alpha=0.05){
  
  dat1  <- reform(Data, arch, id=id)
  X     <- dat1[[1]]
  Nx    <- dat1[[2]]
  nloci <- 2 #ncol(X)

  # MLEs
  out <- strmodel(dat1, arch, BC=BC, method=method, Bbias=Bbias, plugin=plugin)
  trin <- rev((arch)^c(0,1)) # geadic representation
  out2 <- out[[2]]
  rnames1 <- as.integer(rownames(out2)) - 1
  rnames <- rnames1
  nh <- length(rnames)
  dat <- array(0,c(nh,nloci))
  for(k in 1:2){ #for each locus
    re <- rnames%%trin[k]
    dat[,k] <- (rnames-re)/trin[k]
    rnames <- re
  }
  for(i in 1:nh){
      rnames[i] <- paste(dat[i,], collapse = '')
  }
  rownames(out2) <- rnames

  # Bootstrap CIs
  if(CI){
    nhap  <- length(out2)
    N     <- sum(Nx)
    prob  <- Nx/N
    Estim <- array(0, dim = c((nhap+1), B))
    rownames(Estim) <- c('l',(rnames1+1))
    for (l in 1:B){
      infct <- vector(mode = "list", length = 2)
      samp  <- rmultinom(N, 1, prob)
      tmp   <- rowSums(samp)
      pick  <- tmp == 0
      infct[[1]]  <- X[!pick,]
      infct[[2]]  <- tmp[!pick]
      tmp1        <- strmodel(infct, arch, BC=BC, method=method, Bbias=Bbias, plugin=plugin)
      rnames      <- as.integer(rownames(tmp1[[2]]))
      Estim[1,l]  <- unlist(tmp1[[1]])  
      Estim[as.character(rnames),l] <- unlist(tmp1[[2]])                             ## Evaluating and saving the Estimates
    }
    perc <- t(apply(Estim, 1, quantile, c(alpha/2, (1-alpha/2))))
    if(is.null(plugin)){
      out3 <- c(unlist(out[[1]]), perc[1,])
      names(out3) <- c('', paste0(as.character((alpha/2)*100), '%'), paste0(as.character((1-alpha/2)*100), '%')) 
    }else{
      out3 <- out[[1]]
      names(out3) <- c('')
    }
    out4 <- cbind(out2,perc[2:(nhap+1),])

    out <- list(out3, out4, dat)
  }else{
    out1 <- out[[1]]
    names(out1) <- c('')
    out <- list(out1, t(out2), dat)
  }
  names(out) <- c(expression(lambda), 'p', 'haplotypes')
  out
}

mle2 <-function(dat1, arch, id=TRUE, plugin=NULL, CI=FALSE, BC=FALSE, method="bootstrap", Bbias=10000, B=10000, alpha=0.05){
  
  X     <- dat1[[1]]
  Nx    <- dat1[[2]]
  nloci <- 2 #ncol(X)

  # MLEs
  out <- strmodel(dat1, arch, BC=BC, method=method, Bbias=Bbias, plugin=plugin)
  trin <- rev((arch)^c(0,1)) # geadic representation
  out2 <- out[[2]]
  rnames1 <- as.integer(rownames(out2)) - 1
  rnames <- rnames1
  nh <- length(rnames)
  dat <- array(0,c(nh,nloci))
  for(k in 1:2){ #for each locus
    re <- rnames%%trin[k]
    dat[,k] <- (rnames-re)/trin[k]
    rnames <- re
  }
  for(i in 1:nh){
      rnames[i] <- paste(dat[i,], collapse = '')
  }
  rownames(out2) <- rnames

  # Bootstrap CIs
  if(CI){
    nhap  <- length(out2)
    N     <- sum(Nx)
    prob  <- Nx/N
    Estim <- array(0, dim = c((nhap+1), B))
    rownames(Estim) <- c('l',(rnames1+1))
    for (l in 1:B){
      infct <- vector(mode = "list", length = 2)
      samp  <- rmultinom(N, 1, prob)
      tmp   <- rowSums(samp)
      pick  <- tmp == 0
      infct[[1]]  <- X[!pick,]
      infct[[2]]  <- tmp[!pick]
      tmp1        <- strmodel(infct, arch, BC=BC, method=method, Bbias=Bbias, plugin=plugin)
      rnames      <- as.integer(rownames(tmp1[[2]]))
      Estim[1,l]  <- unlist(tmp1[[1]])  
      Estim[as.character(rnames),l] <- unlist(tmp1[[2]])                             ## Evaluating and saving the Estimates
    }
    perc <- t(apply(Estim, 1, quantile, c(alpha/2, (1-alpha/2))))
    if(is.null(plugin)){
      out3 <- c(unlist(out[[1]]), perc[1,])
      names(out3) <- c('', paste0(as.character((alpha/2)*100), '%'), paste0(as.character((1-alpha/2)*100), '%')) 
    }else{
      out3 <- out[[1]]
      names(out3) <- c('')
    }
    out4 <- cbind(out2,perc[2:(nhap+1),])

    out <- list(out3, out4, dat)
  }else{
    out1 <- out[[1]]
    names(out1) <- c('')
    out <- list(out1, t(out2), dat)
  }
  names(out) <- c(expression(lambda), 'p', 'haplotypes')
  out
}
#################################
# The function adhocfreqmodel(X,Nx) calculates the relative prevalence of the haplotypes
# conditionned on unambiguous observations.
#################################
adhocfreqmodel <- function(X, Nx, arch){
  X     <- cbind(X, Nx)
  n     <- ncol(X)
  nloci <- n-1

  # number of haplotypes
  nhpl <- prod(arch)
  
  # extract unambiguous observations
  bin <- list(2^(0:arch[1]), 2^(0:arch[2]))
  Y <-  X[,1:2]+1  
  pick1 <- Y
  for(i in 1:2){
    Y[,i] <- sapply((Y[,i]), function(x) sum(as.integer(intToBits(x))[-seq((arch[i]+1),32)]))
    pick1[,i] <- Y[,i] == 1
  }
  pick <- rowSums(pick1) > 0
  X1   <- matrix(X[pick,], ncol = (nloci+1))

  if(!all(is.na(X1))){  # if there are unambiguous infections
    n1 <- nrow(X1)
    if(is.null(n1)){ # if there is only one unambiguous infection
      n1 <- 1
    }
    X <- matrix(X1, nrow = n1)

    # find indexes of multiple infections
    pick <- rowSums(matrix(pick1[pick,]==0, ncol = 2)) == 1 
    tmp2 <- matrix(X[pick,1:nloci], ncol = nloci)
    idx1 <- which(pick)

    if(length(idx1)>0){
      # single infections
      tmp <- matrix(X[-idx1,], ncol = 3)
      s <- cbind(tmp, tmp[,3])
    }else {
      s <- cbind(X, X[,3])
    }

    # Single infections to corresponding haplotypes
    for(i in 1:2){
      s[,i] <- sapply((s[,i]+1), function(x) which(as.integer(intToBits(x))[-seq((arch[i]+1),32)]==1)-1)
    }

    # Saving single infections
    s1 <- matrix(s[,-4],ncol=3)

    # find all the haplotypes in X
    for(i in idx1){
      y <- matrix(X[i,], ncol = 3, byrow = TRUE)
      idx2 <- which(!(y[c(1,2)]+1)%in%2^(0:max(arch)))
      y2 <- y[idx2]+1
      all <- which(as.integer(intToBits(y2))[-seq((arch[idx2]+1),32)]==1)-1
      nl <- length(all)
      h <- matrix(rep(c(y,y[,3]/nl),length(all)), ncol=4, byrow = TRUE)
      h[,idx2] <- all
      idx3 <- c(1,2)[-idx2] 
      h[,idx3] <- sapply((h[,idx3]+1), function(x) which(as.integer(intToBits(x))[-seq((arch[idx3]+1),32)]==1)-1)
      # add haplotypes in s
      s <- rbind(s,h)
    }
    
    # binary representation
    bin <- c(arch[2],1)
    pp  <- s[,1:2]%*%bin+1
    pp  <- cbind(pp,s[,3:4])
    pp1  <- s1[,1:2]%*%bin+1
    pp1  <- cbind(pp1,s1[,3])

    # observed haplotypes
    idx3 <- as.integer(colnames(t(as.data.frame(summary.factor(pp[,1])))))
    idx5 <- as.integer(colnames(t(as.data.frame(summary.factor(pp1[,1])))))
    
    # Frequencies estimates Method 1
    p <- matrix(0, ncol=length(idx3), nrow = 3)
    tot <- colSums(pp[,2:3])
    for (i in idx3){
      idx4 <- which(pp[,1]==i)
      p[1:2,which(i==idx3)] <- colSums(matrix(pp[idx4,2:3],ncol=2))/tot
    }
    for (i in idx5){
      idx6 <- which(pp1[,1]==i)
      p[3,which(i==idx3)] <- sum(pp1[idx6,2])/sum(pp1[,2])
    }
      colnames(p) <-  idx3
  }else{
    # Frequencies estimates
    p <- matrix(0, ncol=nhpl)
  }

  list(p[1,], p[2,], p[3,])
}

adhocfreqmodelsim <- function(X, Nx, arch){
  
  XX     <- cbind(X, Nx)
  n     <- ncol(X)
  nloci <- n

  # number of haplotypes
  nhpl <- prod(arch)
  
  # extract unambiguous observations
  bin <- list(2^(0:arch[1]), 2^(0:arch[2]))
  Y <-  XX[,1:2]+1  
  pick1 <- Y
  for(i in 1:2){
    Y[,i] <- sapply((Y[,i]), function(x) sum(as.integer(intToBits(x))[-seq((arch[i]+1),32)]))
    pick1[,i] <- Y[,i] == 1
  }
  pick <- rowSums(pick1) > 0
  X1   <- matrix(XX[pick,], ncol = (nloci+1))

  if(!all(is.na(X1))){  # if there are unambiguous infections
    n1 <- nrow(X1)
    if(is.null(n1)){ # if there is only one unambiguous infection
      n1 <- 1
    }
    XX <- matrix(X1, nrow = n1)

    # find indexes of multiple infections
    pick <- rowSums(matrix(pick1[pick,]==0, ncol = 2)) == 1 
    #tmp2 <- matrix(XX[pick,1:nloci], ncol = nloci)
    idx1 <- which(pick)

    if(length(idx1)>0){
      # single infections
      tmp <- matrix(XX[-idx1,], ncol = 3)
      s <- cbind(tmp, tmp[,3])
    }else {
      s <- cbind(XX, XX[,3])
    }

    # Single infections to corresponding haplotypes
    for(i in 1:2){
      s[,i] <- sapply((s[,i]+1), function(x) which(as.integer(intToBits(x))[-seq((arch[i]+1),32)]==1)-1)
    }

    # Saving single infections
    s1 <- matrix(s[,-4],ncol=3)

    # find all the haplotypes in XX
    for(i in idx1){
      y <- matrix(XX[i,], ncol = 3, byrow = TRUE)
      idx2 <- which(!(y[c(1,2)]+1)%in%2^(0:max(arch)))
      y2 <- y[idx2]+1
      all <- which(as.integer(intToBits(y2))[-seq((arch[idx2]+1),32)]==1)-1
      nl <- length(all)
      h <- matrix(rep(c(y,y[,3]/nl),length(all)), ncol=4, byrow = TRUE)
      h[,idx2] <- all
      idx3 <- c(1,2)[-idx2] 
      h[,idx3] <- sapply((h[,idx3]+1), function(x) which(as.integer(intToBits(x))[-seq((arch[idx3]+1),32)]==1)-1)
      # add haplotypes in s
      s <- rbind(s,h)
    }
    
    # binary representation
    bin <- c(arch[2],1)
    pp  <- s[,1:2]%*%bin+1
    pp  <- cbind(pp, matrix(s[,3:4], ncol = 2))  #cbind(pp,s[,3:4])   # frequencies from uambiguous observations, assuming (i) non-weighted, (ii) weighted
    pp1  <- s1[,1:2]%*%bin+1
    pp1  <- cbind(pp1,s1[,3])  # frequencies from single infections

    # observed haplotypes
    idx3 <- as.integer(colnames(t(as.data.frame(summary.factor(pp[,1])))))
    idx5 <- as.integer(colnames(t(as.data.frame(summary.factor(pp1[,1])))))
    
    # Frequencies estimates Method 1
    p <- matrix(0, ncol=length(idx3), nrow = 3)
    tot <- colSums(matrix(pp[,2:3], ncol = 2))
    for (i in idx3){
      idx4 <- which(pp[,1]==i)
      p[1:2,which(i==idx3)] <- colSums(matrix(pp[idx4,2:3],ncol=2))/tot
    }
    for (i in idx5){
      idx6 <- which(pp1[,1]==i)
      p[3,which(i==idx3)] <- sum(pp1[idx6,2])/sum(pp1[,2])
    }
      colnames(p) <-  idx3
  }else{
    # Frequencies estimates
    p <- matrix(0, ncol=nhpl)
  }

  ## Ordering the frequencies 
  p <- list(p[1,], p[2,], p[3,])

  for(i in 1:3){
    pp <- p[[i]]
    pp <- matrix(pp[order(as.integer(idx3))], nrow = 1)
    colnames(pp) <- idx3

    if(length(pp)<nhpl){
      name <- colnames(pp)
      cnt <- 0
      for (j in 1:nhpl) {
        if (!(j%in%name)){
          pp <- append(pp,0.0,after=j-1)
        }
      }
    }
    p[[i]] <- pp
  }
  p
}

#################################
# The function sampl(dat) finds the number of occurences of each observation in the dataset dat.
# The output is a vector of those numbers.
#################################
sampl <- function(dat, arch){
    nloci <- ncol(dat)
    trin <- rev((2^arch-1)^c(0,1)) # vector of geadic representaion
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
estcondprev <- function(estim, arch){
  # This function estimates the unambiguous prevalence as defined in the manuscript of tsoungui et.al, titled
  # "A maximum-likelihood method to estimate haplotype frequencies and prevalence alongside multiplicity of 
  # infection from SNPs data"

    # Number of loci
    numb_Loci <- 2

    # Table of all possible haplotypes
    Hapl <- hapl(arch)

    # Cardinality of Uh
    cardUh <- sum(arch)-1

    cnames <- colnames(t(estim[[2]]))

    ## Access the estimates
    tmp2 <- estim[[2]]
    pickhap <- estim[[3]]%*%c(arch[2],1)+1

    nHapl <- length(pickhap)

    prev <- matrix(0, ncol = nHapl)
    colnames(prev) <- cnames
    numh <- rep(0, nHapl)
    denh <- rep(0, nHapl)

    # Find ambiguous prevalence for each observed haplotype
    for (idx in pickhap[,1]){ 
        # Build uh
        uh <- setUh(Hapl[idx,], cardUh, arch)

        # Indexes of haplotypes forming unambiguous observations with h (Uh)
        pick1 <- uh%*%c(arch[2],1)+1

        ## Pick the right frequencies estimates 
        pickh <- which(pickhap == pick1[1])
        GPh   <- gen_func(tmp2[pickh], estim[[1]])

        # Picking haplotypes in Uh with non zero frequencies
        pick2 <- which(pickhap%in%pick1) # indices of observed haplotypes in uh
        tmp3 <- rep(0, cardUh)
        tmp3[which(pick1%in%pickhap)] <- tmp2[pick2]
        rem <- which(tmp3==tmp2[pickh])[1]
        tmp3 <- tmp3[-rem]
        i <- which(idx==pickhap)
        for(j in 1:(cardUh-1)){ 
            GPartFreq  <- gen_func(tmp3[j], estim[[1]])
            GFreq      <- gen_func(sum(c(tmp2[pickh],tmp3[j])), estim[[1]])
            tmp        <- GFreq - GPartFreq
            numh[i]    <- numh[i] + tmp
            denh[i]    <- denh[i] + tmp/2
        }
        numh[i] <- numh[i] - ((cardUh-1) - 1)*GPh
        denh[i] <- denh[i] - ((cardUh-1)/2 - 1)*GPh
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

#################################
# The function adhocLDmodel(X,Nx) calculates the strength of pairwise Linkage-disequilibrium using the D' r^2, and Q*
# measures for multi-allelic loci.
#################################
adhocLDsim <- function(freq, gen){
  # This function implements the true linkage disequilibrium measures (D', r^2, Q*) as defined in the manuscript of Tsoungui & Schneider, titled
  # "A maximum-likelihood method to estimate haplotype frequencies and prevalence alongside multiplicity of infection from STR data"

  pmat  <- matrix(freq, gen, byrow = TRUE)
  Afreq <- rowSums(pmat)
  Bfreq <- colSums(pmat)

  idxA <- Afreq != 0
  idxB <- Bfreq != 0
  gen_new <- c(sum(idxA), sum(idxB))

  pmatnew  <- pmat[idxA, idxB]
  freqnew  <- matrix(freq, nrow = gen[1], byrow = TRUE)[idxA, idxB]
  Afreqnew <- Afreq[idxA]
  Bfreqnew <- Bfreq[idxB]

  pij   <- Afreqnew%*%t(Bfreqnew)        # pipj
  picjc <- (1-Afreqnew)%*%t(1-Bfreqnew)  # (1-pi)(1-pj)
  picj  <- (1-Afreqnew)%*%t(Bfreqnew)    # (1-pi)pj
  pijc  <- Afreqnew%*%t(1-Bfreqnew)      # pi(1-pj)

  D     <- round(freqnew - pij, 5)       # Dij
  pick  <- D > 0

  D_pos       <- array(0, gen_new)
  D_pos[pick] <- D[pick]

  D_neg        <- array(0, gen_new)
  D_neg[!pick] <- D[!pick]

  Dmax_pos  <- D
  Dmax_neg  <- D

  Dmax_pos  <- do.call(pmin, list(picj, pijc))  # if Dij > 0
  Dmax_neg  <- do.call(pmin, list(pij, picjc))  # if Dij < 0

  pij_pos       <- array(0, gen_new)
  pij_pos[pick] <- pij[pick]

  pij_neg        <- array(0, gen_new)
  pij_neg[!pick] <- pij[!pick]

  Dp <- sum(pij_pos*D_pos/Dmax_pos) + sum(pij_neg*abs(D_neg)/Dmax_neg) # D'
  tmp_sum <- sum(D^2/pij)
  r  <- sum(D^2)/((1-sum(Afreqnew**2))*(1-sum(Bfreqnew**2))) # D* # tmp_sum/min(gen-1)     # r^2
  Q  <- tmp_sum/prod(gen-1)    # Q*
  
  list(Dp, r, Q)
}

ldestim0 <- function(est, gen){
  # Ordering frequencies estimates from 1 to H=n1*n2. The ordering is necessary to obtain the alleles marginal frequencies.
  freq <- array(0, c(1,prod(gen)))
  frequ <- est$p
  trin <- rev((gen)^c(0,1)) #gen^c(1,0)
  hap <- est$haplotypes%*%trin + 1
  colnames(frequ) <- hap
  for(i in hap){
    idx <- which(colnames(frequ) == i)
    freq[,i] <- frequ[,idx]
  }

  pmat  <- matrix(freq, gen, byrow = TRUE)
  Afreq <- rowSums(pmat) # Alleles marginal frequencies locus 1
  Bfreq <- colSums(pmat) # Alleles marginal frequencies locus 2

  idxA <- Afreq != 0
  idxB <- Bfreq != 0
  gen_new <- c(sum(idxA), sum(idxB))

  #pmatnew  <- pmat[idxA, idxB]
  freqnew  <- matrix(freq, nrow = gen[1], byrow = TRUE)[idxA, idxB]
  Afreqnew <- Afreq[idxA]
  Bfreqnew <- Bfreq[idxB]

  pij   <- Afreqnew%*%t(Bfreqnew)        # pipj
  picjc <- (1-Afreqnew)%*%t(1-Bfreqnew)  # (1-pi)(1-pj)
  picj  <- (1-Afreqnew)%*%t(Bfreqnew)    # (1-pi)pj
  pijc  <- Afreqnew%*%t(1-Bfreqnew)      # pi(1-pj)

  D     <- freqnew - pij #round(freqnew - pij, 5)       # Dij
  pick  <- D > 0

  D_pos       <- array(0, gen_new)
  D_pos[pick] <- D[pick]

  D_neg        <- array(0, gen_new)
  D_neg[!pick] <- D[!pick]

  Dmax_pos  <- D
  Dmax_neg  <- D

  Dmax_pos  <- do.call(pmin, list(picj, pijc))  # if Dij > 0
  Dmax_neg  <- do.call(pmin, list(pij, picjc))  # if Dij < 0

  pij_pos       <- array(0, gen_new)
  pij_pos[pick] <- pij[pick]

  pij_neg        <- array(0, gen_new)
  pij_neg[!pick] <- pij[!pick]

  Dp <- sum(pij_pos*D_pos/Dmax_pos) + sum(pij_neg*abs(D_neg)/Dmax_neg) # D'
  tmp_sum <- sum(D^2/pij)
  r  <- sum(D^2)/((1-sum(Afreqnew**2))*(1-sum(Bfreqnew**2))) # D* # tmp_sum/min(gen-1)     # r^2
  Q  <- tmp_sum/prod(gen-1)    # Q*
  
  list(Dp, r, Q)
}

ldestim <- function(Data, arch, plugin=NULL, id=TRUE, CI=FALSE, B=10000, alpha=0.05){
  dat1  <- reform(Data, arch, id=id)
  X     <- dat1[[1]]
  Nx    <- dat1[[2]]
  #nloci <- 2 #ncol(X)

  # MLEs
  est <- mle(Data, arch, id=FALSE, plugin=NULL, CI=FALSE, BC=FALSE, method="bootstrap", Bbias=100, B=100, alpha=0.05) #strmodel(dat1, arch, BC=FALSE, method=method, Bbias=FALSE, plugin=plugin)
<<<<<<< HEAD
  ldvals <- ldestim0(est, arch)
=======
  ld <- ldestim0(est, arch)
>>>>>>> a3368880b8bfa8163092b26e9f16092c7cf7e2b7
  # Bootstrap CIs
  if(CI){
    N     <- sum(Nx)
    prob  <- Nx/N
    Estim <- array(0, dim = c(3, B))
    for (l in 1:B){
      infct <- vector(mode = "list", length = 2)
      samp  <- rmultinom(N, 1, prob)
      tmp   <- rowSums(samp)
      pick  <- tmp == 0
      infct[[1]]  <- X[!pick,]
      infct[[2]]  <- tmp[!pick]

      #infct <- Data[samp,]
      esttmp1     <- mle2(infct, arch, id=FALSE, plugin=NULL, CI=FALSE, BC=FALSE, method="bootstrap", Bbias=100, B=100, alpha=0.05)
      tmp1        <- ldestim0(esttmp1, arch)
      Estim[,l]  <- unlist(tmp1) 
    }
<<<<<<< HEAD
    Estim <- Estim[ , colSums(is.na(Estim))==0]
=======
>>>>>>> a3368880b8bfa8163092b26e9f16092c7cf7e2b7
    perc <- t(apply(Estim, 1, quantile, c(alpha/2, (1-alpha/2))))
    out <- cbind(ldvals,perc)
    rownames(out) <- c("D'", expression(r^2), "Q")
  }else{
<<<<<<< HEAD
    out <- ldvals
  }
=======
    out <- ld
    #names(out) <- c("D'", expression(r^2), "Q")
  }
 # names(out) <- c("D'", expression(r^2), "Q")
>>>>>>> a3368880b8bfa8163092b26e9f16092c7cf7e2b7
  out
}
