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

cpoiss<-function(lambda,n){
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

## Haplotype frequency distribution
pp1 <- abs(rnorm(16,0,1))
pp1 <- pp1/sum(pp1)

P <- pp1
lambda <- 1
k <- 250
n <- 4

ddd <- sampNew(P,lambda,k, n) #Generated data (vector containing the number of time the infections are present in the dataset,
                                               #Dataset of the present infections in their trinary representation)


#---------------------------------------

varsets <- function(l,n){   #calculate all var sets
  # n number of loci
  # l number of observations per locus
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

#---------------------------------------
# this calculates the list to pick proper haplotype freuencies, i.e. sets Ay
X <- array(c(0,2,2,2,1,2,1,2),c(2,4)) ### data

X <- ddd[[1]]
la <- 2
Nx <- c(2,3)  ## Nx
Nx <- ddd[[2]]
eps <- 10^-8
N <- sum(Nx)
nn <- nrow(X)  # of different observations
Ax <- list()
n <- ncol(X)  # loci
#cardAx <- array(0,nn)
for(u in 1:nn){
  xx <- array(X[u,],c(1,n))
  sel <- (1:n)[xx==2]
  l <- length(sel)
  #cardAx[u] <-l 
  if(l==0){
    yy <- xx
  }else{
    yy <- xx[rep(1,3^l),]
    yy[,sel] <- varsets(3,l)
  }
  bin <- 2^(0:(n-1))
  iilist <- list()
  siglist <- list()
  for(i in 1:3^l){
    y1 <- array(yy[i,],c(1,n))
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
pp <- array(rep(1/H,H),c(H,1))
rownames(pp) <- hapl1

#t <- 0
#initial list#

num0 <- pp*0
cond1 <- 1  ## condition to stop EM alg! 

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
      denom <- denom + exlap-vz  ##   = (1-)^(Nx-Ny)*(Exp(lambda*sum p)-1) = (1-)^(Nx-Ny)*G(sum p)
      num[Ax[[u]][[1]][[k]],] <- num[Ax[[u]][[1]][[k]],]+ exlap#*pp[Ax[[u]][[1]][[k]],]
      ## exlap =  (1-)^(Nx-Ny) G'(sum p)   --- denominator of generating functions cancels out!
      CC <- CC + exlap*p
    }
    num <- num*pp
    #print("____________")
    denom <- Nx[u]/denom
    denom <- la*denom
    Ccoeff <- Ccoeff + CC*denom
   
    Bcoeff <- Bcoeff + num*denom
    
  }
  
  Ccoeff <- Ccoeff/N
  #print(Bcoeff/sum(Bcoeff))
  #print(Ccoeff)
  ppn <- Bcoeff/(sum(Bcoeff))
  
  ### Newton step
  cond2 <- 1
  xt <- Ccoeff   ### good initial condition
  while(cond2 > eps){
    ex <- exp(-xt)
    xtn <- xt + (1-ex)*(xt + Ccoeff*ex - Ccoeff)/(ex*xt+ex-1)
    cond2 <- abs(xtn-xt)
    xt <- xtn
  }
  cond1 <- abs(xt-la)+sqrt(sum((pp-ppn)^2))
  print(c(cond1,abs(xt-la),sqrt(sum((pp-ppn)^2))))
  la <- xt
  pp <- ppn
  #print(c(pp))
  print(la)
  print(like(pp,la,Nx,N,Ax))
  
}


like <- function(pp,la,Nx,N,Ax){
  logli <- 0
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
      denom <- denom + (exlap-vz) 
    }
    
    logli <-  logli+ Nx[u]*log(denom/(exp(la)-1))
  }
 logli 
  
    
}  

#------------------------------


xx <- array(c(1,2,2,2),c(1,4))
n <- length(xx)
sel <- (1:n)[xx==2]
l <- length(sel)
if(l==0){
  yy <- xx
}else{
  yy <- xx[rep(1,3^l),]
  yy[,sel] <- varsets(3,l)
}
yy

bin <- 2^(0:(n-1))
iilist <- list()
siglist <- list()
for(i in 1:3^l){
  y1 <- yy[i,]
  sel <- (1:n)[y1==2]
  l1 <- length(sel)
  if(l1==0){
    ii <- y1
  }else{
    ii <- yy[rep(1,2^l1),]
    ii[,sel] <- varsets(2,l1)
  }
  iilist[[i]] <- ii%*%bin+1
  siglist[[i]] <- l1
  
}

psum <- iilist

psum <- function(inlist,pp,len){
  p <- 0
  for(i in len){
    p <- p+ sum(pp[inlist[[i]]])
  }
}


 

l <- length(A1)
B <- subsets(l)  # all subsetes where 2 has to be replaced
pick <- rowSums(B)
out <- array(,c(3^l-1,n))
t <- 1


subsets <- function(l){
  B <- array(FALSE,c(2^l,l))
  B[1:2,1] <- c(FALSE,TRUE)
  for(k in 2:l){
    B[(2^(k-1)+1):2^k,1:(k-1)] <- B[1:2^(k-1),1:(k-1)]
    B[(2^(k-1)+1):2^k,k] <- TRUE
  }
  B
}
