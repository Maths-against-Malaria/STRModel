# Title        : Simulation study using simulated data
# Objectives   : Implement the EM-algorithm on simulated data and save the estimates
# Created by   : christian Tsoungui Obama, Kristan A. Schneider
# Created on   : 03.04.21
# Last modified: 15.04.22

# Relative path
path <- "/Volumes/GoogleDrive-117934057836063832284/My Drive/Maths against Malaria/Christian/Models/MultiLociBiallelicModel/"

# Loading external ressources
source(paste0(path, "src/nbiallelicModel.R"))          ## Loading Model
source(paste0(path, "src/dataGenerator.R"))            ## Loading the data generaor for model

# True Poisson parameter
lbdavec <- c(0.1, 0.25, 0.5, 1, 1.5, 2, 2.5)
NLbd <- length(lbdavec)

# Number of loci considered
NumbLoci <- c(2, 5)
Nn <- length(NumbLoci)
nvec <- matrix(NumbLoci, nrow = 1, ncol = Nn)

# Number of possible haplotypes
Hvec <- 2^nvec
NH <- length(Hvec)
Hvecpo <- Hvec + 1

# Sample sizes considered
Nvec <- c(50, 100, 150, 200, 500)
NN <- length(Nvec)

# Number of estimates generated in the simulation
NEst <- 10000

# Number of frequencies set per (number of loci) case
NFreq <- length(NumbLoci)

# Extra parameters
ParExtra <- list(NLbd, Nn, Hvec, NN, NEst, NFreq)

# True haplotype frequencies
Pvec <- vector(mode="list", length=Nn)
for (i in 1:Nn){
  Pvec[[i]] <- array(1/Hvec[i], c(1, Hvec[i]))
  Pvec[[i]] <- rbind(Pvec[[i]], c(0.7, rep(0.3/(Hvec[i]-1), (Hvec[i]-1))))
}

# True parameter
True_param <- list(Pvec, lbdavec, Nvec)

# Simulation
out <- vector(mode = "list", length = Nn)
out2 <- vector(mode = "list", length = Nn)

for (i in 1:Nn){
  print(paste0("processing frequency distributions of ", nvec[,i], " loci."))
  sizelist <- vector(mode = "list", length = NN)
  sizelist2 <- vector(mode = "list", length = NN)
  for (j in 1:NN){                                                                                ## For each value of the sample size
    lbdalist <- vector(mode = "list", length = NLbd)
    lbdalist2 <- vector(mode = "list", length = NLbd)
    for (k in 1:NLbd){                                                                            ## For each value of the lambda parameter
      Estim <- array(0, dim = c(Hvecpo[i], NEst, NFreq))
      adhocEstim <- array(0, dim = c(Hvec[i], NEst, NFreq))
      for (cnt in 1:NFreq){
        for (l in 1:NEst){
          infct <- sampNew(unlist(Pvec[[i]][cnt,]) ,unlist(lbdavec[k]) ,Nvec[j], nvec[,i])        ## Generating data for the simulation
          Estim[,l,cnt] <- unlist(nbialModel(infct[[2]], infct[[1]]))                             ## Evaluating and saving the Estimates
          adhocEstim[,l,cnt] <- unlist(adhocModel(infct))                                    ## Ad hoc estimates for frequencies
        }
      }
      lbdalist[[k]] <- Estim
      lbdalist2[[k]] <- adhocEstim
    }
    sizelist[[j]] <- lbdalist
    sizelist2[[j]] <- lbdalist2
  }
  out[[i]] <- sizelist
  out2[[i]] <- sizelist2

  # Saving the MOI and frequencies estimates for post-processing
  saveRDS(out, file = paste0(path, "dataset/modelEstimates.rds"))
  saveRDS(out2, file = paste0(path, "dataset/adhocModelEstimates.rds"))

  # End of simulation warning
  print(paste0("Simulation finished for ", nvec[,i], " loci ;)"))
}

# Saving the true parameters
saveRDS(True_param, file = paste0(path, "dataset/true_Parameters.rds"))

# Saving the extra parameters
saveRDS(ParExtra, file = paste0(path, "dataset/extra_Parameters.rds"))
