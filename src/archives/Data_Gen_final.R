hapl <- function(n){
  H <- array(0,c(2^n,n))
  H[1:2,1] <- c(0,1)
  for(k in 2:n){
    H[(2^(k-1)+1):2^k,1:(k-1)] <- H[1:2^(k-1),1:(k-1)]
    H[(2^(k-1)+1):2^k,k] <- 1
  }
  H <- H[,n:1]
  H
}

obs <- function(M){
  n <- ncol(as.matrix(M)) #number of loci
  if(n==1){
    M <- t(M)
    n <- ncol(M)
  }
  out <- array(0,n)
  for(k in 1:n){
    out[k] <- sum(unique(M[,k])+1)
  }
  out
}

cpoissc<-function(lambda,n){
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

# generate example data
sampNew <- function(P,lambda,k, n){ #P = haplotype distro., lambda = Poisson parameter, K = Sample size, n = N° of loci
  H <- hapl(n) #Set of possible haplotypes
  vec <- 3^(0:(n-1)) #Vector for trianary representation
  out <- array(0,k)
  m <- cpoiss(lambda,k) #MOI values for each sample following CPoiss(lambda)
  for(j in 1:k){# for each sample
    s <- rmultinom(1, m[j], P) #multinomially select M[j] haplotypes from the haplotype pool
    out[j] <- sum(vec*(obs(H[s!=0,])-1))+1 #Summing up the trianary representation of a number representing the infection
  } #vector of infections
  out <- t(as.data.frame(summary.factor(out))) #vector of how many times each infection that is effectively present appears in the dataset
  vals <- as.integer(colnames(out))-1 #Infections present in the dataset
  dat <- array(0,c(length(vals),n))
  for(k in 0:(n-1)){ #for each locus
    re <- vals%%(3^(n-k-1))
    dat[,n-k] <- (vals-re)/(3^(n-k-1))
    vals <- re
  } #Trianary representation of each infection present in the dataset
  list(dat,c(out)) 
}