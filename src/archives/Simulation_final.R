#************************************************************************************
############################### Simulation ##########################################
#************************************************************************************

simul <- function(lambda ,P ,Nsz ,np, nn){
  
N <- c()
freq <- c()
cv_psi <- c()
real_psi <- c()
Rbias_p <- c()
Rbias_m <- c()

for (i in 1:length(P)){               ## For each value of the frequency distribution
 
  BB <- c()
  CC <- c()
  fr <- i
  
  for (j in 1:length(Nsz)){           ## For each value of the sample size
    Ns <- Nsz[j]
    
    for (k in 1:length(lambda)){         ## For each value of the lambda parameter
    
      Estim <- matrix(0, nrow = np, ncol = ss)
      
      for (ii in 1:ss){              ## Generating ss datasets
        infct <-  sampNew(unlist(P[i]) ,unlist(lambda[k]) ,Ns  ,nn)#as.data.frame(list(samp(P[,i], lambda[k], Ns, nl))) ## Generating data for the simulation 
        Estim[,ii] <- unlist(nbialModel(infct[[2]], infct[[1]]))  ## Evaluating and saving the Estimates
        } 
      
      est_moi  <- Estim[1,]/(1 - exp(-Estim[1,]))                     ## Estimated values of MOI
      Real_moi <- lambda[k]/(1 - exp(-lambda[k]))                     ## Real value of MOI
      bias_m   <- ((est_moi/Real_moi) - 1)*100                          ## Bias of MOI estimates
      Rbias_m  <- c(Rbias_m, mean(bias_m))                            ## Mean relative bias for MOI estimates
      
      p <- P[,i]                                                      ## Real values of the frequencies
      est_p    <- Estim[2:np,]                                        ## Estimated values of the frequencies  
      
      bias_p   <- ((est_p/p) - 1)*100                                   ## Relative bias for frequencies estimates
      Rbias_p  <- c(Rbias_p, mean(bias_p))                            ## Mean relative bias for frequencies estimates
      
      cv_psi   <- c(cv_psi, sd( (est_moi/lambda[k])*100 ))            ## Coefficient of variation for MOI

      N        <- c(N, Nsz[j])
      freq     <- c(freq, fr)
      real_psi <- c(real_psi, lambda[k])
    }
    
  }

}

FDat <- data.frame(lambda = real_psi, frequency = freq, Sample = N, Bias_moi = Rbias_m, 
                   Variation_coef_moi = cv_psi, Bias_freq = Rbias_p)                   ## Storing the estimates in a dataframe

write.csv(FDat, "Dataset/Model2_simulated_data.csv")       ## Saving the estimates for further analysis

}

#************************************************************************************
################################# Simulation ########################################

source("n_biallelic_model_final.R")          ## Loading Model2
source("Data_Gen_final.R")  ## Loading the data generaor for model2

######################################## MAIN #########################################


################################ Generating the Data #################################

# Defining simulation values for the poisson parameter
lambda <- c(0.1, 0.25, 0.5, 1, 1.5, 2, 2.5)

# Number of loci
n <- 2

# Number of haplotypes
nh <- 2^n

# Generating a set of 4 p vectors
P <- as.data.frame(cbind(rep(1/nh, nh), c(0.6, rep(0.4/(nh - 1), nh - 1)) ))
Np <- nrow(P)+1     

# Set of sample sizes
Nsz <- c(50, 100, 150, 200, 500)

# Defining the number of datasets generated in the simulation
ss <- 100

simul(lambda, P, Nsz, Np, n)
