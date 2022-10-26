# Title        : Simulation study using simulated data
# Objectives   : Implement the EM-algorithm on simulated data and save the estimates
# Created by   : christian Tsoungui Obama, Kristan A. Schneider
# Created on   : 05.05.22
# Last modified: 24.05.22

# Relative path
path <- "/Volumes/GoogleDrive-117934057836063832284/My Drive/Maths against Malaria/Christian/Models/STRModel/"

# Loading external ressources
source(paste0(path, "src/STRModel.R"))          ## Loading Model
#source(paste0(path, "src/dataGenerator.R"))            ## Loading the data generaor for model

# True Poisson parameter
lbdavec <- c(0.1, 0.25, 0.5, 1, 1.5, 2, 2.5)
NLbd <- length(lbdavec)

# Genetic architecture
genArch <- matrix(c(2,2,2,3,4,7), ncol = 2, byrow = TRUE)
Nn <- nrow(genArch)
#nvec <- matrix(NumbLoci, nrow = 1, ncol = Nn)

# Number of possible haplotypes
Hvec <-  apply(genArch,1, prod) 
NH <- length(Hvec)
Hvecpo <- Hvec + 1

# Sample sizes considered
Nvec <- c(50, 100, 150, 200, 500)
NN <- length(Nvec)

# Number of estimates generated in the simulation
NEst <- 10000

# Number of frequencies set per (number of loci) case
NFreq <- 2

# Extra parameters
ParExtra <- list(NLbd, Nn, Hvec, NN, NEst, NFreq, genArch)

# True haplotype frequencies
Pvec <- vector(mode="list", length=Nn)
for (i in 1:Nn){
  Pvec[[i]] <- array(1/Hvec[i], c(1, Hvec[i]))
  Pvec[[i]] <- rbind(Pvec[[i]], c(0.7, rep(0.3/(Hvec[i]-1), (Hvec[i]-1))))
}

# True parameter
True_param <- list(Pvec, lbdavec, Nvec)

# Simulation
out  <- vector(mode = "list", length = Nn)

outfreq1 <- vector(mode = "list", length = Nn)
outfreq2 <- vector(mode = "list", length = Nn)
outfreq3 <- vector(mode = "list", length = Nn)

outDp1 <- vector(mode = "list", length = Nn)
outDp2 <- vector(mode = "list", length = Nn)

outr1 <- vector(mode = "list", length = Nn)
outr2 <- vector(mode = "list", length = Nn)

outQ1 <- vector(mode = "list", length = Nn)
outQ2 <- vector(mode = "list", length = Nn)

for (i in 1:Nn){
  print(paste0("processing frequency distributions for m=", genArch[i,1], " and n=",  genArch[i,2], " alleles, respectively."))
  sizelist <- vector(mode = "list", length = NN)

  sizefreqlist1 <- vector(mode = "list", length = NN)
  sizefreqlist2 <- vector(mode = "list", length = NN)
  sizefreqlist3 <- vector(mode = "list", length = NN)

  sizeDplist1 <- vector(mode = "list", length = NN)
  sizeDplist2 <- vector(mode = "list", length = NN)

  sizerlist1 <- vector(mode = "list", length = NN)
  sizerlist2 <- vector(mode = "list", length = NN)

  sizeQlist1 <- vector(mode = "list", length = NN)
  sizeQlist2 <- vector(mode = "list", length = NN)
  for (j in 1:NN){                                                                                ## For each value of the sample size
    lbdalist <- vector(mode = "list", length = NLbd)

    lbdalist1 <- vector(mode = "list", length = NLbd)
    lbdalist2 <- vector(mode = "list", length = NLbd)
    lbdalist3 <- vector(mode = "list", length = NLbd)

    lbdaDplist1 <- vector(mode = "list", length = NLbd)
    lbdaDplist2 <- vector(mode = "list", length = NLbd)

    lbdarlist1 <- vector(mode = "list", length = NLbd)
    lbdarlist2 <- vector(mode = "list", length = NLbd)

    lbdaQlist1 <- vector(mode = "list", length = NLbd)
    lbdaQlist2 <- vector(mode = "list", length = NLbd)
    for (k in 1:NLbd){                                                                            ## For each value of the lambda parameter
      Estim <- array(0, dim = c(Hvecpo[i], NEst, NFreq))

      adhocFreqEstim1 <- array(0, dim = c(Hvec[i], NEst, NFreq))
      adhocFreqEstim2 <- array(0, dim = c(Hvec[i], NEst, NFreq))
      adhocFreqEstim3 <- array(0, dim = c(Hvec[i], NEst, NFreq))

      adhocDpEstim1 <- array(0, dim = c(1, NEst, NFreq))
      adhocDpEstim2 <- array(0, dim = c(1, NEst, NFreq))

      adhocrEstim1 <- array(0, dim = c(1, NEst, NFreq))
      adhocrEstim2 <- array(0, dim = c(1, NEst, NFreq))

      adhocQEstim1 <- array(0, dim = c(1, NEst, NFreq))
      adhocQEstim2 <- array(0, dim = c(1, NEst, NFreq))
      for (cnt in 1:NFreq){
        for (l in 1:NEst){
          infct              <- datagen(unlist(Pvec[[i]][cnt,]) ,unlist(lbdavec[k]) ,Nvec[j], genArch[i,])      ## Generating data for the simulation
          Estim[,l,cnt]      <- unlist(strmodel(infct, genArch[i,]))                                            ## Evaluating and saving the Estimates
          tmp <- adhocfreqmodelsim(infct[[1]], infct[[2]], genArch[i,])        ## Ad hoc estimates for frequencies

          adhocFreqEstim1[,l,cnt] <- tmp[[1]]                                  ## Ad hoc estimates for frequencies
          adhocFreqEstim2[,l,cnt] <- tmp[[2]]                                  ## Ad hoc estimates for frequencies assuming unambiguous infections
          adhocFreqEstim3[,l,cnt] <- tmp[[3]]                                  ## Ad hoc estimates for frequencies assuming single infections

          tmp2 <- adhocLDsim(tmp[[2]], genArch[i,])
          adhocDpEstim1[,l,cnt] <- tmp2[[1]]                                    ## Ad hoc estimates for frequencies assuming unambiguous infections
          adhocrEstim1[,l,cnt]  <- tmp2[[2]]                                    ## Ad hoc estimates for frequencies assuming single infections
          adhocQEstim1[,l,cnt]  <- tmp2[[3]]                                    ## Ad hoc estimates for frequencies assuming single infections

          tmp3 <- adhocLDsim(tmp[[3]], genArch[i,])
          adhocDpEstim2[,l,cnt] <- tmp3[[1]]                                    ## Ad hoc estimates for frequencies assuming unambiguous infections
          adhocrEstim2[,l,cnt]  <- tmp3[[2]]                                    ## Ad hoc estimates for frequencies assuming single infections
          adhocQEstim2[,l,cnt]  <- tmp3[[3]]                                    ## Ad hoc estimates for frequencies assuming single infections
        }
      }
      lbdalist[[k]]  <- Estim
      lbdalist1[[k]] <- adhocFreqEstim1
      lbdalist2[[k]] <- adhocFreqEstim2
      lbdalist3[[k]] <- adhocFreqEstim3

      lbdaDplist1[[k]] <- adhocDpEstim1
      lbdaDplist2[[k]] <- adhocDpEstim2

      lbdarlist1[[k]] <- adhocrEstim1
      lbdarlist2[[k]] <- adhocrEstim2

      lbdaQlist1[[k]] <- adhocQEstim1
      lbdaQlist2[[k]] <- adhocQEstim2
    }
    sizelist[[j]]  <- lbdalist
    sizefreqlist1[[j]] <- lbdalist1
    sizefreqlist2[[j]] <- lbdalist2
    sizefreqlist3[[j]] <- lbdalist3

    sizeDplist1[[j]] <- lbdaDplist1
    sizeDplist2[[j]] <- lbdaDplist2

    sizerlist1[[j]] <- lbdarlist1
    sizerlist2[[j]] <- lbdarlist2

    sizeQlist1[[j]] <- lbdaQlist1
    sizeQlist2[[j]] <- lbdaQlist2
  }
  #out[[i]]  <- sizelist
  #outfreq1[[i]] <- sizefreqlist1
  #outfreq2[[i]] <- sizefreqlist2
  #outfreq3[[i]] <- sizefreqlist3

  #outDp1[[i]] <- sizeDplist1
  #outDp2[[i]] <- sizeDplist2
  #
  #outr1[[i]] <- sizerlist1
  #outr2[[i]] <- sizerlist2

  #outQ1[[i]] <- sizeQlist1
  #outQ2[[i]] <- sizeQlist2

  # Saving the MOI and frequencies estimates for post-processing
  saveRDS(sizelist,      file = paste0(path, "dataset/full_modelEstimates_", i,".rds"))
  saveRDS(sizefreqlist1, file = paste0(path, "dataset/full_adhocModelEstimates1_", i,".rds"))
  saveRDS(sizefreqlist2, file = paste0(path, "dataset/full_adhocModelEstimates2_", i,".rds"))
  saveRDS(sizefreqlist3, file = paste0(path, "dataset/full_adhocModelEstimates3_", i,".rds"))

  saveRDS(sizeDplist1, file = paste0(path, "dataset/full_adhocModelDpEstimates1_", i,".rds"))
  saveRDS(sizeDplist2, file = paste0(path, "dataset/full_adhocModelDpEstimates2_", i,".rds"))

  saveRDS(sizerlist1, file = paste0(path, "dataset/full_adhocModelREstimates1_", i,".rds"))
  saveRDS(sizerlist2, file = paste0(path, "dataset/full_adhocModelREstimates2_", i,".rds"))

  saveRDS(sizeQlist1, file = paste0(path, "dataset/full_adhocModelQEstimates1_", i,".rds"))
  saveRDS(sizeQlist2, file = paste0(path, "dataset/full_adhocModelQEstimates2_", i,".rds"))

  # End of simulation warning
  print(paste0("Simulation finished ;)"))
}

# Saving the true parameters
saveRDS(True_param, file = paste0(path, "dataset/full_true_Parameters.rds"))

# Saving the extra parameters
saveRDS(ParExtra, file = paste0(path, "dataset/full_extra_Parameters.rds"))
