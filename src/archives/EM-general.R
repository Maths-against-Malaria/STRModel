pp1 <- abs(rnorm(16,0,1))
pp1 <- pp1/sum(pp1)
#lambda <- 2

sampNew <- function(P,lambda,k,n){
  H <- hapl(n) 
  vec <- 3^(0:(n-1))
  out <- array(0,k)
  m <- cpoiss(lambda,k)
  for(j in 1:k){
    s <- rmultinom(1, m[j], P)
    out[j] <- sum(vec*(obs(H[s!=0,])-1))+1
  }
  out <- t(as.data.frame(summary.factor(out)))
  vals <- as.integer(colnames(out))-1
  dat <- array(0,c(length(vals),n))
  for(k in 0:(n-1)){
    re <- vals%%(3^(n-k-1))
    dat[,n-k] <- (vals-re)/(3^(n-k-1))
    vals <- re
  }
  list(dat,c(out))  
}

ddd <- sampNew(pp1,1,250,4)
#---------------------------------------
# help functions 

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

varsets <- function(l,n){   #calculate all var sets
  # n number of loci
  # l number of alleles per locus
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

gead <- function(x,l,n){   ## calculates geadic expression of each element of vectorx 
  l <- rep(l,n)
  out <- array(0,c(length(x),n))
  div <- c(1,cumprod(l[1:(n-1)]))
  for(k in n:1){
    r <- x%%div[k]
    #print(c(r,div[k]))
    out[,k] <- (x-r)/div[k]
    x <- r
  }
  out
}

gead1 <- function(x,l){   ## calculates general geadic expression of each element of vectorx
  n <- length(l)
  out <- array(0,c(length(x),n))
  div <- c(1,cumprod(l[1:(n-1)]))
  for(k in n:1){
    r <- x%%div[k]
    print(c(r,div[k]))
    out[,k] <- (x-r)/div[k]
    x <- r
  }
  out
}

#_____________________________
# example data
X <- dat[[1]] + 1 # array(c(3,0,7,2,2,1),c(2,3))
Nx <- dat[[2]]
N <- sum(Nx)
nn <- nrow(X)
l <- arch #c(2,3)   # number of alleles per locus
n <- length(l)
x <- X

#_____________________________
# this calculates the list to pick proper haplotype freuencies, i.e. sets Ay

#allconf <- function(x,l,n){  # x array
hapll <- c()
if(length(l)==1){
  l <- rep(l,n)
}else{
  n <- length(l)
}
ggead <- c(1,cumprod(l[1:(n-1)]))
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
    Hx[[u]][[4]][[k]] <- temp1%*%bin2num[[k]] # subsets of alleles
    Hx[[u]][[5]][[k]] <- Hx[[u]][[3]][[k]]%*%rep(1,l[k]) # number of alleles in subset at locus k
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




#}
#---------------------------------------
# initialize parameters

H <- length(hapl1)
pp <- array(rep(1/H,H),c(H,1))
rownames(pp) <- hapl1

#t <- 0
#initial list#
eps <- 10^(-8)
num0 <- pp*0
cond1 <- 1  ## condition to stop EM alg! 
la <- 2
num <- num0
Bcoeff <- num0
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
    for(k in 1:Ax[[u]][[1]]){
      p <- sum(pp[Ax[[u]][[4]][[k]],])
      vz <- Ax[[u]][[3]][[k]]
      lap <- la*p
      exlap <- vz*exp(lap)
      denom <- denom + exlap-vz  ##   = (1-)^(Nx-Ny)*(Exp(lambda*sum p)-1) = (1-)^(Nx-Ny)*G(sum p)
      num[Ax[[u]][[4]][[k]],] <- num[Ax[[u]][[4]][[k]],]+ exlap#*pp[Ax[[u]][[1]][[k]],]
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
  #print(like(pp,la,Nx,N,Ax))
  
}

likeGen <- function(pp,la,Nx,N,Ax){
  logli <- 0
  Bcoeff <- num0 #reset B coefficients to 0 in next iteration
  num <- num0  #reset numerator to 0 in next iteration
  for(u in 1:nn){
    denom <- 0
    num <- num0
    CC <- 0
    for(k in 1:Ax[[u]][[1]]){
      p <- sum(pp[Ax[[u]][[4]][[k]],])
      vz <- Ax[[u]][[3]][[k]]
      lap <- la*p
      exlap <- vz*exp(lap)
      denom <- denom + (exlap-vz) 
    }
    
    logli <-  logli+ Nx[u]*log(denom/(exp(la)-1))
  }
 logli 
  
    
}