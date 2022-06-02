# Title        : Simulation study
# Objective    : Compute the bias and coefficient of variation of the estimates,
#                and the estimates of prevalence
# Created by   : Christian Tsoungui Obama
# Created on   : 05.05.22
# Last modified: 23.05.22

# Importing libraries
library(wordspace)

# Functions
psi <- function (inp){
  inp / (1 - exp(-inp))
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

bias     <- function(estim_Param, sim_Param, name){
  # This function implements bias of mean MOI and haplotype frequencies as defined in the manuscript
  # "A maximum-likelihood method to estimate haplotype frequencies and prevalence 
  # alongside multiplicity of infection from SNPs data"

  moibiasloc <- vector(mode = "list", length = n_Sim_Loci)
  freqbiasloc <- vector(mode = "list", length = n_Sim_Loci)

  for (l in 1:n_Sim_Loci){   # For each number of locus
    freqbiasSamp <- vector(mode = "list", length = n_Sampl)
    moibiasSamp  <- vector(mode = "list", length = n_Sampl)

    # Number of haplotypes
    num_Hapl   <- n_Hapl[l]
    num_Hapl_PlusOne <- num_Hapl + 1

    for (k in 1:n_Sampl){  # For each sample size
      freqbias_lamb <- vector(mode = "list", length = n_Freq_Distr)
      moibias_lamb  <- vector(mode = "list", length = n_Freq_Distr)

      for (j in 1:n_Lbda){ # For each true Lambda
        freq_bias <- vector(mode = "list", length = n_Freq_Distr)
        moi_bias  <- vector(mode = "list", length = n_Freq_Distr)

        for (i in 1:n_Freq_Distr){ # For each choice of true frequency distribution
          # Access each of the 10000 estimates
          tmp1 <- estim_Param[[l]][[k]][[j]][, , i]

          # True frequencies
          tmp2 <- t(sim_Param[[1]][[l]])
          tmp3 <- tmp2[, i]

          # Relative bias of haplotype freq. in percent
          freqbias <- rowMeans((tmp1[2:num_Hapl_PlusOne,]/tmp3 - 1), na.rm = TRUE) * 100

          # Remove the estimates of Lambda that yield an "Infinite mean MOI"
          mean_moi_estim <- psi(tmp1[1, ])
          
          # Relative bias of MOI in percent
          bias_moi <- (mean_moi_estim/true_Mean_MOI[j] - 1)
          bias_moi <- mean(bias_moi[!is.infinite(bias_moi)], na.rm = TRUE) * 100

          # Save frequency and MOI bias in lists
          freq_bias[[i]] <- freqbias
          moi_bias[[i]]  <- bias_moi
        }
        freqbias_lamb[[j]] <- freq_bias
        moibias_lamb[[j]]  <- moi_bias
      }
      freqbiasSamp[[k]] <- freqbias_lamb
      moibiasSamp[[k]]  <- moibias_lamb
    }
    freqbiasloc[[l]] <- freqbiasSamp
    moibiasloc[[l]]  <- moibiasSamp
  }

  # Saving the estimates
  saveRDS(freqbiasloc, file = paste0(path, "dataset/freqbias",name ,".rds"))
  saveRDS(moibiasloc, file = paste0(path, "dataset/moibias",name ,".rds"))
}

coefvar  <- function(estim_Param, sim_Param, name){
  # This function implements coefficient of variation of mean MOI as defined in the manuscript
  # "A maximum-likelihood method to estimate haplotype frequencies and prevalence 
  # alongside multiplicity of infection from SNPs data"

  moicvloc <- vector(mode = "list", length = n_Sim_Loci)

  for (l in 1:n_Sim_Loci){   # For each number of locus
    moicvSamp  <- vector(mode = "list", length = n_Sampl)

    # Number of haplotypes
    num_Hapl   <- n_Hapl[l]
    num_Hapl_PlusOne <- num_Hapl + 1

    for (k in 1:n_Sampl){  # For each sample size
      moicv_lamb  <- vector(mode = "list", length = n_Freq_Distr)

      for (j in 1:n_Lbda){ # For each true Lambda
        moi_cv  <- vector(mode = "list", length = n_Freq_Distr)

        for (i in 1:n_Freq_Distr){ # For each choice of true frequency distribution
          # Access each of the 10000 estimates
          tmp1 <- estim_Param[[l]][[k]][[j]][, , i]

          # Remove the estimates of Lambda that yield an "Infinite mean MOI"
          mean_moi_estim <- psi(tmp1[1, ])
          mean_moi_estim <- mean_moi_estim[!is.infinite(mean_moi_estim)]
          
          # Dimensionless coefficient of variation of MOI
          moicv <- sd(mean_moi_estim, na.rm = TRUE) / true_Mean_MOI[j]

          # Save moi coefficient of variation in list
          moi_cv[[i]]  <- moicv
        }
        moicv_lamb[[j]]  <- moi_cv
      }
      moicvSamp[[k]]  <- moicv_lamb
    }
    moicvloc[[l]]  <- moicvSamp
  }

  # Saving the estimates
  saveRDS(moicvloc, file = paste0(path, "dataset/moicv", name, ".rds"))
}

true_amb_prevalence  <- function(reshap_Sim_Param, name){
  # This function implements the true ambiguous prevalence as defined in the manuscript of tsoungui et.al, titled
  # "A maximum-likelihood method to estimate haplotype frequencies and prevalence alongside multiplicity of infection from SNPs data"

  qh_loc <- vector(mode = "list", length = n_Sim_Loci)

  for (l in 1:n_Sim_Loci){
    qh_loc[[l]] <- vector(mode = "list", length = n_Freq_Distr)

    # Number of haplotypes
    num_Hapl         <- n_Hapl[l]
    num_Hapl_PlusOne <- n_Hapl[l] + 1

    for (i in 1:n_Freq_Distr){ # For each choice of true frequency distribution
      amb_prevalence <- array(0, dim = c(num_Hapl, n_Lbda))

      # Access each of the 10.000 parameters estimates
      tmp1 <- reshap_Sim_Param[[l]][[i]]

      # For each set of parameters estimates
      for (j in 1:n_Lbda){
        amb_prevalence[,j] <- (exp(tmp1[1,j]) - exp(1-tmp1[2:num_Hapl_PlusOne,j])^tmp1[1,j])/(exp(tmp1[1,j])-1)
      }

      qh_loc[[l]][[i]] <- amb_prevalence
    }
  }

  # Save the ambiguous prevalence estimates
  saveRDS(qh_loc, file = paste0(path, "dataset/true_Amb_Prevalence", name, ".rds"))
  qh_loc
}

true_conditional_prevalence <- function(reshap_Sim_Param, sim_Param, name, gen){
  # This function implements the true conditional prevalence as defined in the manuscript of tsoungui et.al, titled
  # "A maximum-likelihood method to estimate haplotype frequencies and prevalence alongside multiplicity of infection from SNPs data"
 
  qh_loc <- vector(mode = "list", length = n_Sim_Loci)

  for (l in 1:n_Sim_Loci){   # For each number of locus
    qh_loc[[l]] <- vector(mode = "list", length = n_Freq_Distr)
    true_freq <- sim_Param[[1]][[l]]

    # Number of haplotypes
    num_Hapl <- n_Hapl[l]

    # Number of loci
    numb_Loci <- 2

    # Table of all possible haplotypes
    arch <- gen[l,]
    Hapl <- hapl(arch)

    ## Cardinality of the set Uh
    cardUh <- sum(arch)-1

    for (i in 1:n_Freq_Distr){ # For each choice of true frequency distribution
      numh <- array(0, dim = c(num_Hapl, n_Lbda))
      denh <- array(0, dim = c(num_Hapl, n_Lbda))
      rh   <- array(0, dim = c(num_Hapl, n_Lbda))

      ## For each combination of true lanmbda and haplotype frequencies values
      tmp1 <- reshap_Sim_Param[[l]][[i]]
      tmp2 <- tmp1[2:(num_Hapl+1),]

      ## For eachobserved haplotype build the set Uh for all l
      trufreq_vec <- true_freq[i,]
      pickhap     <- which(trufreq_vec != 0)

      for (idx in pickhap){ 
        uh <-  setUh(Hapl[idx,], cardUh, arch)

        ## Pick the right frequencies estimates 
        pickh <- which(colSums(uh[1,] == t(Hapl))==numb_Loci)
        GPh   <- gen_func(tmp2[pickh,], tmp1[1,])

        GPartFreq <- rep(0, n_Lbda)
        GFreq     <- rep(0, n_Lbda)

        pick1 <- rep(0, cardUh)
        pick2 <- rep(0, cardUh)

        for(idxUh in 1:cardUh){ 
          pick1[idxUh] <- which(colSums(uh[idxUh,] == t(Hapl))==numb_Loci)
        }

        pick2    <- pick1
        pick2[1] <- 0

        for(idxUh in 2:cardUh){ 
          GPartFreq  <- gen_func(tmp2[pick2[idxUh],], tmp1[1,])
          GFreq      <- gen_func(colSums(tmp2[pick1[c(1,idxUh)],]), tmp1[1,])
          tmp        <- GFreq - GPartFreq
          numh[idx,] <- numh[idx,] + tmp
          denh[idx,] <- denh[idx,] + tmp/2
        }
          numh[idx,] <- numh[idx,] - ((cardUh-1) - 1)*GPh
          denh[idx,] <- denh[idx,] - ((cardUh-1)/2 - 1)*GPh
      }

      den <- colSums(denh, na.rm = TRUE)

      for(q in 1:n_Lbda){
        rh[,q] <- numh[,q]/den[q]
      }
      qh_loc[[l]][[i]] <- rh
    }
  }  
  # Saving the cunditional prevalence estimates
  saveRDS(qh_loc, file = paste0(path,  "dataset/true_Cond_Prevalence", name, ".rds"))
  qh_loc
}

true_relative_prevalence    <- function(reshap_Sim_Param, sim_Param, name, gen){
  # This function implements the true relative prevalence as defined in the manuscript of tsoungui et.al, titled
  # "A maximum-likelihood method to estimate haplotype frequencies and prevalence alongside multiplicity of infection from SNPs data"
 
  qh_loc <- vector(mode = "list", length = n_Sim_Loci)

  for (l in 1:n_Sim_Loci){   # For each number of locus
    qh_loc[[l]] <- vector(mode = "list", length = n_Freq_Distr)
    true_freq   <- sim_Param[[1]][[l]]

    # Number of haplotypes
    num_Hapl <- n_Hapl[l]

    # Number of loci
    numb_Loci <- 2

    # Table of all possible haplotypes
    arch <- gen[l,]
    Hapl <- hapl(arch)

    ## Cardinality of the set Uh
    cardUh <- sum(arch)-1

    for (i in 1:n_Freq_Distr){ # For each choice of true frequency distribution
      rh <- array(0, dim = c(num_Hapl, n_Lbda))

      ## For each combination of true lanmbda and haplotype frequencies values
      tmp1 <- reshap_Sim_Param[[l]][[i]]
      tmp2 <- tmp1[2:(num_Hapl+1),]

      ## For each haplotype build the set Uh for all l
      trufreq_vec <- true_freq[i,]
      pickhap <- which(trufreq_vec != 0)

      for (idx in pickhap){ 
        uh <-  setUh(Hapl[idx,], cardUh, arch)

        ## Pick the right frequencies estimates 
        pickh <- which(colSums(uh[1,] == t(Hapl))==numb_Loci)
        GPh   <- gen_func(tmp2[pickh,], tmp1[1,])

        GPartFreq <- rep(0, n_Lbda)
        GFreq     <- rep(0, n_Lbda)

        pick1 <- rep(0, cardUh)
        pick2 <- rep(0, cardUh)

        for(idxUh in 1:cardUh){ 
          pick1[idxUh] <- which(colSums(uh[idxUh,] == t(Hapl))==numb_Loci)
        }

        pick2    <- pick1
        pick2[1] <- 0

        for(idxUh in 2:cardUh){ 
          GPartFreq <- gen_func(tmp2[pick2[idxUh],], tmp1[1,])
          GFreq     <- gen_func(colSums(tmp2[pick1[c(1,idxUh)],]), tmp1[1,])
          rh[idx,]  <- rh[idx,] + GFreq - GPartFreq
        }
          rh[idx,]  <- rh[idx,] - ((cardUh-1) - 1)*GPh
      }
    for (q in 1:n_Lbda){
      rh[,q] <- rh[,q]/sum(rh[,q], na.rm = TRUE)
    }
    qh_loc[[l]][[i]] <- rh
    }
  }
  # Saving the relative prevalence estimates
  saveRDS(qh_loc, file = paste0(path,  "dataset/true_Rel_Prevalence", name, ".rds"))
  qh_loc
}

estim_amb_prevalence <- function(estim_Param, true_prev, name){
  # This function estimates the ambiguous prevalence as defined in the manuscript of tsoungui et.al, titled
  # "A maximum-likelihood method to estimate haplotype frequencies and prevalence alongside multiplicity of infection from SNPs data"

  qh_loc <- vector(mode = "list", length = n_Sim_Loci)

  for (l in 1:n_Sim_Loci){   # For each number of locus
    qh_Samp      <- vector(mode = "list", length = n_Sampl)

    # Number of haplotypes
    num_Hapl         <- n_Hapl[l]
    num_Hapl_PlusOne <- num_Hapl + 1

    for (k in 1:n_Sampl){  # For each true sample size
      qh_lamb <- vector(mode = "list", length = n_Freq_Distr)

      for (j in 1:n_Lbda){ # For each true Lambda
        qh_freq <- vector(mode = "list", length = n_Freq_Distr)

        for (i in 1:n_Freq_Distr){ # For each choice of true frequency distribution
          amb_prevalence <- array(0, dim = c(num_Hapl, n_Sampl_Gen))

          ## Access each of the 10000 estimates
          tmp1 <- estim_Param[[l]][[k]][[j]][, , i]

          ## For each set of estimates, compute prevalence
          amb_prevalence <- (exp(tmp1[1,]) - exp(1-tmp1[2:num_Hapl_PlusOne,])^tmp1[1,])/(exp(tmp1[1,])-1)

          ## Replace entries with NAN values by 0
          amb_prevalence[is.na(amb_prevalence)] <- 0.0

          qh <- rowMeans(amb_prevalence, na.rm = TRUE)
          
          ## Save the prevalence in a list
          qh_freq[[i]] <- qh
        }
        qh_lamb[[j]] <- qh_freq
      }
      qh_Samp[[k]] <- qh_lamb
    }
    qh_loc[[l]] <- qh_Samp
  }
  # Saving the estimates
  saveRDS(qh_loc, file = paste0(path, "dataset/estim_Amb_Prevalence", name, ".rds"))
  qh_loc
}

estim_conditional_prevalence <- function(estim_Param, sim_Param, true_prev, name, gen){
  # This function estimates the unambiguous prevalence as defined in the manuscript of tsoungui et.al, titled
  # "A maximum-likelihood method to estimate haplotype frequencies and prevalence alongside multiplicity of infection from SNPs data"
 
  qh_loc <- vector(mode = "list", length = n_Sim_Loci)

  for (l in 1:n_Sim_Loci){   # For each number of locus
    qh_Samp <- vector(mode = "list", length = n_Sampl)

    # True haplotype frequencies
    true_freq <- sim_Param[[1]][[l]]

    for (k in 1:n_Sampl){  # For each true sample size
      qh_lamb <- vector(mode = "list", length = n_Freq_Distr)

      # Number of haplotypes
      num_Hapl <- n_Hapl[l]

      # Number of loci
      numb_Loci <- 2

      # Table of all possible haplotypes
      arch <- gen[l,]
      Hapl <- hapl(arch)

      ## Cardinality of the set Uh
      cardUh <- sum(arch)-1

      for (j in 1:n_Lbda){ # For each true Lambda
        qh_freq      <- vector(mode = "list", length = n_Freq_Distr)

        for (i in 1:n_Freq_Distr){ # For each choice of true frequency distribution
          rh   <- array(0, dim = c(num_Hapl, n_Sampl_Gen))
          numh <- array(0, dim = c(num_Hapl, n_Sampl_Gen))
          denh <- array(0, dim = c(num_Hapl, n_Sampl_Gen))
          prev <- array(0, dim = c(num_Hapl, n_Sampl_Gen))

          prevalence <- rep(0, num_Hapl)

          ## Access each of the 10000 estimates
          tmp1 <- estim_Param[[l]][[k]][[j]][, , i]
          tmp2 <- tmp1[2:(num_Hapl+1),]

          ## For each haplotype build the set Uh for all l
          trufreq_vec <- true_freq[i,]
          pickhap <- which(trufreq_vec != 0)

          # Find ambiguous prevalence for each haplotype
          for (idx in pickhap){ 
            uh <-  setUh(Hapl[idx,], cardUh, arch)

            ## Pick the right frequencies estimates 
            pickh <- which(colSums(uh[1,] == t(Hapl))==numb_Loci)
            GPh   <- gen_func(tmp2[pickh,], tmp1[1,])

            GPartFreq <- rep(0, n_Lbda)
            GFreq     <- rep(0, n_Lbda)

            pick1 <- rep(0, cardUh)
            pick2 <- rep(0, cardUh)

            for(idxUh in 1:cardUh){ 
              pick1[idxUh] <- which(colSums(uh[idxUh,] == t(Hapl))==numb_Loci)
            }

            pick2 <- pick1
            pick2[1] <- 0

            for(idxUh in 2:cardUh){ 
              GPartFreq  <- gen_func(tmp2[pick2[idxUh],], tmp1[1,])
              GFreq      <- gen_func(colSums(tmp2[pick1[c(1,idxUh)],]), tmp1[1,])
              tmp        <- GFreq - GPartFreq
              numh[idx,] <- numh[idx,] + tmp
              denh[idx,] <- denh[idx,] + tmp/2
            }
              numh[idx,] <- numh[idx,] - ((cardUh -1) - 1)*GPh
              denh[idx,] <- denh[idx,] - ((cardUh -1)/2 - 1)*GPh
          }

          den <- colSums(denh, na.rm = TRUE)
          for (q in 1:n_Sampl_Gen){
            prev[,q] <- numh[,q]/den[q]
          }
          prevalence <- rowMeans(prev, na.rm = TRUE)

          ## Save the prevalence in a list
          qh_freq[[i]]  <- prevalence

        }
        qh_lamb[[j]] <- qh_freq
      }
      qh_Samp[[k]]   <- qh_lamb
    }
    qh_loc[[l]]      <- qh_Samp
  }
  # Saving the estimates
  saveRDS(qh_loc, file = paste0(path, "dataset/estim_Cond_Prevalence", name, ".rds"))
  qh_loc
}

estim_relative_prevalence <- function(estim, name, indx){
  # This function estimates the relative prevalence as defined in the manuscript of tsoungui et.al, titled
  # "A maximum-likelihood method to estimate haplotype frequencies and prevalence alongside multiplicity of infection from SNPs data"

  qh_loc <- vector(mode = "list", length = n_Sim_Loci)

  for (l in 1:n_Sim_Loci){   # For each number of locus
    qh_Samp <- vector(mode = "list", length = n_Sampl)

    # True haplotype frequencies
    true_freq <- sim_Param[[1]][[l]]

    for (k in 1:n_Sampl){  # For each true sample size
      qh_lamb <- vector(mode = "list", length = n_Freq_Distr)

      # Number of haplotypes
      num_Hapl <- n_Hapl[l]

      for (j in 1:n_Lbda){ # For each true Lambda
        qh_freq <- vector(mode = "list", length = n_Freq_Distr)

        for (i in 1:n_Freq_Distr){ # For each choice of true frequency distribution
          prev       <- array(0, dim = c(num_Hapl, n_Sampl_Gen))
          prevalence <- rep(0, num_Hapl)

          ## Access each of the 10000 estimates
          prev <- estim[[l]][[k]][[j]][, , i]
          prevalence <- rowMeans(prev, na.rm = TRUE)
  
          ## Save the prevalence in a list
          qh_freq[[i]]      <- prevalence
        }
        qh_lamb[[j]]      <- qh_freq
      }
      qh_Samp[[k]]        <- qh_lamb
    }
    qh_loc[[l]]      <- qh_Samp
  }

  # Saving the estimates
  saveRDS(qh_loc, file = paste0(path,  "dataset/estim_Rel_Prevalence",indx, name, ".rds"))
  qh_loc
}

main <- function(sim_Param, reshap_Sim_Param, name, gen){
  # Loading estimated haplotype frequencies and MOI
  estim_Param       <- readRDS(paste0(path, "dataset/modelEstimates", name, ".rds"))
  adhoc_estim_Param1 <- readRDS(paste0(path, "dataset/adhocModelEstimates1", name, ".rds"))
  adhoc_estim_Param2 <- readRDS(paste0(path, "dataset/adhocModelEstimates2", name, ".rds"))
  adhoc_estim_Param3 <- readRDS(paste0(path, "dataset/adhocModelEstimates3", name, ".rds"))

  # Bias of frequencies and MOI
  bias(estim_Param, sim_Param, name)

  # Coefficient of variation of MOI
  coefvar(estim_Param, sim_Param, name)

  # True ambiguous prevalence
  true_Amb_Prev          <- true_amb_prevalence(reshap_Sim_Param, name)

  # True conditional prevalence
  true_Conditional_Prev  <- true_conditional_prevalence(reshap_Sim_Param, sim_Param, name, gen)

  # True relative prevalence
  true_Relative_Prev     <- true_relative_prevalence(reshap_Sim_Param, sim_Param, name, gen)


  # Estimated ambiguous prevalence
  estim_amb_prevalence(estim_Param, true_Amb_Prev, name)

  # Estimated conditional prevalence
  estim_conditional_prevalence(estim_Param, sim_Param, true_Conditional_Prev, name, gen)

  # Estimated relative prevalence (adhoc Model)
  estim_relative_prevalence(adhoc_estim_Param1, name, 1)
  estim_relative_prevalence(adhoc_estim_Param2, name, 2)
  estim_relative_prevalence(adhoc_estim_Param3, name, 3)
}
 
# Relative path
path <- "/Volumes/GoogleDrive-117934057836063832284/My Drive/Maths against Malaria/Christian/Models/STRModel/"

# Define data origin ('' <- simulated data, 'Kenya' <- kenyan data)
namelist <- c('') #c('', 'Kenya')

for (name in namelist){
  print(paste0('Ongoing simulation for: ', name, ' data!'))

  # Loading true haplotype parameters for the simulation
  sim_Param <- readRDS(paste0(path, "dataset/true_Parameters", name, ".rds"))

  # Loading extra parameters
  extra_Sim_Param  <- readRDS(paste0(path, "dataset/extra_Parameters", name, ".rds"))

  # Simulation parameters
  n_Lbda        <- extra_Sim_Param[[1]]
  n_Sim_Loci    <- extra_Sim_Param[[2]]
  n_Hapl        <- extra_Sim_Param[[3]]
  n_Sampl       <- extra_Sim_Param[[4]]
  n_Sampl_Gen   <- extra_Sim_Param[[5]]
  n_Freq_Distr  <- extra_Sim_Param[[6]]
  genArch       <- extra_Sim_Param[[7]]
  true_Mean_MOI <- psi(sim_Param[[2]])

  # Reformatting true parameters to compute true prevalence
  reshap_Sim_Param <- vector(mode='list', length=n_Sim_Loci)

  for (i in 1:n_Sim_Loci){
    numb_Row <- n_Hapl[i]+1
    reshap_Sim_Param[[i]] <- vector(mode='list', length=n_Freq_Distr)

    for (j in 1:n_Freq_Distr){
      reshap_Sim_Param[[i]][[j]]             <- array(0, c(numb_Row, n_Lbda))
      reshap_Sim_Param[[i]][[j]][1,]         <- sim_Param[[2]]
      reshap_Sim_Param[[i]][[j]][2:numb_Row,] <- sim_Param[[1]][[i]][j,]
    }
  }

  main(sim_Param, reshap_Sim_Param, name, genArch)
}
