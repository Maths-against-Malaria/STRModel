# Title        : Plotting the hapl. freq. and MOI Bias & CV, and also prevalence
# Objective    : Plot the bias and coefficient of variation of the estimates
# Created by   : christian Tsoungui Obama
# Created on   : 03.04.21
# Last modified: 10.08.22

# Importing libraries
library(dplyr)
library(ggplot2)

# Functions
beautify <- function (p, legende1, legende2, pos, colpal, linety, lgdetitle, form){
  p <- p + theme(panel.grid.minor = element_blank(),panel.grid=element_blank())
  p <- p + theme(panel.grid.major = element_blank(),panel.grid=element_blank())
  p <- p + theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank())
  p <- p + theme(panel.grid.minor.x=element_blank(),panel.grid.major.x=element_blank())
  p <- p + theme(panel.background = element_rect(colour='black',fill="transparent"))
  p <- p + theme(plot.background = element_rect(fill = "transparent", colour = NA))

  # Title
  p <- p + theme(plot.title =element_text(face ="italic",size=rel(2.2),hjust = 0.5))

  # Legend
  p <- p + theme(legend.key=element_blank(), legend.text.align = 0, legend.background = element_blank())
  p <- p + theme(legend.key=element_blank(), legend.text.align = 0, legend.position = pos, legend.background = element_blank())
  p <- p + scale_linetype_manual(values=linety, labels=Ilab <-  legende2, guide = guide_legend(title = NULL))
  p <- p + scale_colour_manual(values=colpal, labels=Ilab <-  legende1, guide = guide_legend(title = lgdetitle))
  p <- p + theme(legend.text = element_text(size = rel(2.2)))
  p <- p + theme(legend.title = element_text(size = rel(1.5),face=form), legend.margin = margin(t = 1, b = 0.00001))
  p <- p + theme(legend.key.width = unit(9.5,"mm"))

  # Axis
  p <- p + theme(axis.text = element_text(colour='black'))
  p <- p + theme(axis.text = element_text(size = rel(2.1)))
  p <- p + theme(axis.title = element_text(size = rel(2.2)))
  p <- p + theme(axis.title.x = element_text(face="italic"))
  p <- p + theme(axis.title.y = element_text(face="plain"))
  p <- p + theme(axis.ticks = element_line(color = "black"))

  p
}

psi <- function (inp){
  inp / (1 - exp(-inp))
}

dataframe_builder_prev <- function(prev_estim, type_prev, locNumb, true_prev){
  NRow <- length(samp_Vec)*length(lbda_Vec)*n_Freq_Distr*n_Hapl [locNumb] 

  cnames <- c('prev', 'sample', 'freq', 'shape')
  df1 <- array(0, dim = c(NRow, length(cnames)))
  df1 <- as.data.frame(df1)
  colnames(df1) <- cnames
  samp_vec <- rep(samp_Vec, each=NRow/(length(samp_Vec)*n_Freq_Distr))
  df1[,'sample'] <- as.factor(rep(samp_vec, n_Freq_Distr))
  df1[,'freq']   <- as.factor(rep(1:n_Hapl [locNumb], NRow/(length(lbda_Vec)*n_Hapl [locNumb])))
  df1[,'shape']  <- as.factor(rep(c("sym", "asym"), each=NRow/n_Freq_Distr))

  exp_prev <- c()
  for (l in 1:n_Freq_Distr){
    for (k in 1:length(samp_Vec)) {
      for (j in 1:n_Lbda){
        for (i in 1:n_Hapl [locNumb]){
          exp_prev <- c(exp_prev, prev_estim[[locNumb]][[k]][[j]][[l]][i])
        }
      }
    }
  }
  df1[,'prev'] <- exp_prev
  df1$type <- as.factor(type_prev)
  df1$vers <- as.factor('estimate')

  df2 <- array(0, dim = c(NRow, length(cnames)))
  df2 <- as.data.frame(df2)
  colnames(df2) <- cnames
  samp_vec <- rep(samp_Vec, each=NRow/(length(samp_Vec)*n_Freq_Distr))
  df2[,'sample'] <- as.factor(rep(samp_vec, n_Freq_Distr))
  df2[,'freq']   <- as.factor(rep(1:n_Hapl [locNumb], NRow/(length(lbda_Vec)*n_Hapl [locNumb])))
  df2[,'shape']  <- as.factor(rep(c("sym", "asym"), each=NRow/n_Freq_Distr))

  exp_prev <- c()
  for (l in 1:n_Freq_Distr){
    for (k in 1:length(samp_Vec)) {
      for (j in 1:n_Lbda){
        exp_prev <- c(exp_prev, true_prev[[locNumb]][[l]][,j])
      }
    }
  }

  df2[,'prev'] <- exp_prev
  df2$type <- as.factor(type_prev)
  df2$vers <- as.factor('true')

  df <- rbind(df1, df2)
  df
}

dataframe_builder_LD <- function(ld_estim, type_ld, locNumb, true_ld){
  NRow <- length(samp_Vec)*length(lbda_Vec)*n_Freq_Distr*n_Hapl[locNumb] 

  cnames <- c('ld', 'sample', 'freq', 'shape')
  df1 <- array(0, dim = c((NRow/n_Hapl[locNumb]), length(cnames)))
  df1 <- as.data.frame(df1)
  colnames(df1) <- cnames
  samp_vec <- rep(samp_Vec, each=NRow/(length(samp_Vec)*n_Freq_Distr*n_Hapl[locNumb]))
  df1[,'sample'] <- as.factor(rep(samp_vec, n_Freq_Distr))
  df1[,'freq']   <- 1 #as.factor(rep(1:n_Hapl [locNumb], NRow/(length(lbda_Vec)*n_Hapl [locNumb])))
  df1[,'shape']  <- as.factor(rep(c("sym", "asym"), each=NRow/(n_Freq_Distr*n_Hapl[locNumb])))

  exp_prev <- c()
  for (l in 1:n_Freq_Distr){
    for (k in 1:length(samp_Vec)) {
      for (j in 1:n_Lbda){
        exp_prev <- c(exp_prev, ld_estim[[locNumb]][[k]][[j]][[l]][1])
      }
    }
  }

  df1[,'ld'] <- exp_prev
  df1$type <- as.factor(type_ld)
  df1$vers <- as.factor('estimate')

  df2 <- array(0, dim = c(NRow/n_Hapl[locNumb], length(cnames)))
  df2 <- as.data.frame(df2)
  colnames(df2) <- cnames
  samp_vec <- rep(samp_Vec, each=NRow/(length(samp_Vec)*n_Freq_Distr*n_Hapl[locNumb]))
  df2[,'sample'] <- as.factor(rep(samp_vec, n_Freq_Distr))
  df2[,'freq']   <- 1
  df2[,'shape']  <- as.factor(rep(c("sym", "asym"), each=NRow/(n_Freq_Distr*n_Hapl[locNumb])))

  exp_prev <- c()
  for (l in 1:n_Freq_Distr){
    for (k in 1:length(samp_Vec)) {
      for (j in 1:n_Lbda){
        exp_prev <- c(exp_prev, true_ld[[locNumb]][[l]][1])
      }
    }
  }

  df2[,'ld'] <- exp_prev
  df2$type <- as.factor(type_ld)
  df2$vers <- as.factor('true')

  df <- rbind(df1, df2)
  df
}

dataframe_builder_Freqperf <- function(perf_estim, locNumb){
  NRow <- length(samp_Vec)*length(lbda_Vec)*n_Freq_Distr*n_Hapl [locNumb]

  # frequencies
  cnames <- c('bias', 'sample', 'freq', 'shape')
  df1 <- array(0, dim = c(NRow, length(cnames)))
  df1 <- as.data.frame(df1)
  colnames(df1) <- cnames

  samp_vec <- rep(samp_Vec, each=NRow/(length(samp_Vec)*n_Freq_Distr))
  df1[,'sample'] <- as.factor(rep(samp_vec, n_Freq_Distr))
  df1[,'freq']   <- as.factor(rep(1:n_Hapl [locNumb], NRow/(length(lbda_Vec)*n_Hapl [locNumb])))
  df1[,'shape']  <- as.factor(rep(c("sym", "asym"), each=NRow/n_Freq_Distr))

  exp_perf <- c()
  for (l in 1:n_Freq_Distr){
    for (k in 1:length(samp_Vec)) {
      for (j in 1:n_Lbda){
          exp_perf <- c(exp_perf, perf_estim[[locNumb]][[k]][[j]][[l]])
      }
    }
  }
  
  df1[,'bias'] <- exp_perf
  df1$type <- as.factor('freq')
  df1
}

dataframe_builder_Moiperf <- function(perf_estim, locNumb){
  NRow <- length(samp_Vec)*length(lbda_Vec)*n_Freq_Distr

  # frequencies
  cnames <- c('bias', 'sample', 'shape')
  df1 <- array(0, dim = c(NRow, length(cnames)))
  df1 <- as.data.frame(df1)
  colnames(df1) <- cnames

  samp_vec <- rep(samp_Vec, each=NRow/(length(samp_Vec)*n_Freq_Distr))
  df1[,'sample'] <- as.factor(rep(samp_vec, n_Freq_Distr))
  df1[,'shape']  <- as.factor(rep(c("sym", "asym"), each=NRow/n_Freq_Distr))

  exp_perf <- c()
  for (l in 1:n_Freq_Distr){
    for (k in 1:length(samp_Vec)) {
      for (j in 1:n_Lbda){
          exp_perf <- c(exp_perf, perf_estim[[locNumb]][[k]][[j]][[l]])
      }
    }
  }
  
  df1[,'bias'] <- exp_perf
  df1$type <- as.factor('moi')
  df1
}

main <- function(sim_Param, name, gen){
  # Plots parameters
  shape_typ <- c('sym', 'asym')

  if(name=="Kenya"){
    dir   <- 'DD'
  }else{
    dir   <- 'SD'
  }

  # Color palette (color-blind friendly) for the plots
  cbPalette <- c(rgb(0,0,0), rgb(.35, .70, .90), rgb(.90,.60,0), rgb(0,.60,.50), rgb(0,.45,.70), rgb(.80,.40,0), rgb(.5, .5, .5))
  lty       <- c("dashed", "solid")
  lty2       <- c("solid", "dashed", "twodash")
  legende2  <- c('estimate', 'true')

  if(1==1){ # Plotting prevalence
    # Importing the data to plot
    amb_prev           <- readRDS(paste0(path, "dataset/estim_Amb_Prevalence",   name, ".rds"))
    relative_prev1     <- readRDS(paste0(path, "dataset/estim_Rel_Prevalence1",   name, ".rds"))
    conditional_prev   <- readRDS(paste0(path, "dataset/estim_Cond_Prevalence",  name, ".rds"))

    true_amb_prev         <- readRDS(paste0(path, "dataset/true_Amb_Prevalence",   name, ".rds"))
    true_relative_prev    <- readRDS(paste0(path, "dataset/true_Rel_Prevalence",   name, ".rds"))
    true_conditional_prev <- readRDS(paste0(path, "dataset/true_Cond_Prevalence",  name, ".rds"))

    # Plots parameters
    legende1 <- c('unobservable', 'conditional', 'relative')

    # Position of legend
    pos <- c(0.20, 0.70)

    for(l in 1:n_Sim_Loci){ # 2 or 5 loci
      # Building the prevalence dataframe
      df_ambprev    <- dataframe_builder_prev(amb_prev,         'amb_prev',         l, true_amb_prev)
      df_relprev1   <- dataframe_builder_prev(relative_prev1,    'relative_prev1',    l, true_relative_prev)
      df_condprev   <- dataframe_builder_prev(conditional_prev, 'conditional_prev', l, true_conditional_prev)

      df <- rbind(df_ambprev, df_condprev, df_relprev1)
      tru_freq <- sim_Param[[1]][[l]]

      for(k in 1:n_Freq_Distr){  # sym or asym
        trufreq_vec <- tru_freq[k,]

        for(i in 1:n_Hapl[l]){
          if(trufreq_vec[i] != 0){
              for(j in samp_Vec){
                df1 <- df %>%
                      filter(sample == j, freq == i, shape == shape_typ[k]) %>%
                      droplevels()

                df1$psiVec <- psi(lbda_Vec)
                
                p <- ggplot(data = df1, aes(x=psiVec))
                p <- p + geom_line(aes(y = df1[,'prev'], color = type, linetype = vers), size=1.)
                p <- p + geom_hline(yintercept = round(trufreq_vec[i], 3), linetype="dashed", color = "grey")
                p <- beautify(p, legende1, legende2, pos, cbPalette, lty, NULL, NULL)
                if(name == 'Kenya'){
                    p <- p + labs(x=expression(frac(lambda, 1 - e^-lambda)), y="", title=paste0("p = ", round(trufreq_vec[i], 3), ", N = ", j, ", year = ", estim_Years[k]))
                    outfile <- paste0(path,"plots/Prev_plots_", dir, "/prev_freq_", i, "_SSize_", j, "_year_", estim_Years[k], "_", name, ".pdf")
                }else{
                    p <- p + labs(x=expression(frac(lambda, 1 - e^-lambda)), y="", title=paste0("p = ", round(trufreq_vec[i], 3), ", N = ", j, ", m = ", gen[l,1], ", n = ", gen[l,2]))
                    outfile <- paste0(path,"plots/Prev_plots_", dir, "/prev_", shape_typ[k], "_freq_", i, "_SSize_", j, "_nloci_", gen[l,1],"_",gen[l,2], "_", name, ".pdf")
                }

                pdf(outfile, height=5, width=8)
                print(p)
                dev.off()
              }
          }
        }
      }
    }
  }

  if(1==1){ # Plotting relative vs weighted vs single
    # Importing the data to plot
    relative_prev1     <- readRDS(paste0(path, "dataset/estim_Rel_Prevalence1",   name, ".rds"))
    relative_prev2     <- readRDS(paste0(path, "dataset/estim_Rel_Prevalence2",   name, ".rds"))
    relative_prev3     <- readRDS(paste0(path, "dataset/estim_Rel_Prevalence3",   name, ".rds"))

    true_relative_prev    <- readRDS(paste0(path, "dataset/true_Rel_Prevalence",   name, ".rds"))

    # Plots parameters
    legende1 <- c('Relative', 'Weighted', 'Single')

    # Position of legend
    pos <- c(0.18, 0.80)

    for(l in 1:n_Sim_Loci){ # 2 or 5 loci
      # Building the prevalence dataframe
      df_relprev1   <- dataframe_builder_prev(relative_prev1, 'relative_prev1', l, true_relative_prev) # frequency estimates from framework
      df_relprev2   <- dataframe_builder_prev(relative_prev2, 'relative_prev2', l, true_relative_prev) # frequency estimates weighted
      df_relprev3   <- dataframe_builder_prev(relative_prev3, 'relative_prev3', l, true_relative_prev) # frequency estimates assuming only single observations

      df <- rbind(df_relprev1, df_relprev2, df_relprev3)
      tru_freq <- sim_Param[[1]][[l]]

      for(k in 1:n_Freq_Distr){  # sym or asym
        trufreq_vec <- tru_freq[k,]

        for(i in 1:n_Hapl[l]){
          if(trufreq_vec[i] != 0){
              for(j in samp_Vec){
                df1 <- df %>%
                      filter(sample == j, freq == i, shape == shape_typ[k], vers == 'estimate') %>%
                      droplevels()

                df1$psiVec <- psi(lbda_Vec)
                
                p <- ggplot(data = df1, aes(x=psiVec))
                p <- p + geom_line(aes(y = df1[,'prev'], color = type), size=1.)
                p <- p + geom_hline(yintercept = round(trufreq_vec[i], 3), linetype="dashed", color = "grey")
                p <- beautify(p, legende1, legende2, pos, cbPalette, lty, NULL, NULL)
                if(name == 'Kenya'){
                    p <- p + labs(x=expression(frac(lambda, 1 - e^-lambda)), y="", title=paste0("p = ", round(trufreq_vec[i], 3), ", N = ", j, ", year = ", estim_Years[k]))
                    outfile <- paste0(path,"plots/Prev_plots_", dir, "/prev_freq_", i, "_SSize_", j, "_year_", estim_Years[k], "_", name, ".pdf")
                }else{
                    p <- p + labs(x=expression(frac(lambda, 1 - e^-lambda)), y="", title=paste0("p = ", round(trufreq_vec[i], 3), ", N = ", j, ", m = ", gen[l,1], ", n = ", gen[l,2]))
                    outfile <- paste0(path,"plots/Prev_plots_", dir, "/prev_adhocVS_", shape_typ[k], "_freq_", i, "_SSize_", j, "_nloci_", gen[l,1],"_",gen[l,2], "_", name, ".pdf")
                }

                pdf(outfile, height=5, width=8)
                print(p)
                dev.off()
              }
          }
        }
      }
    }
  }

  if(1==1){ # Plotting LD
    # Importing the data to plot
    D_ld      <- readRDS(paste0(path, "dataset/estim_LD_D",   name, ".rds"))
    r_ld      <- readRDS(paste0(path, "dataset/estim_LD_r",   name, ".rds"))
    Q_ld      <- readRDS(paste0(path, "dataset/estim_LD_Q",   name, ".rds"))

    true_D_ld    <- readRDS(paste0(path, "dataset/true_LD_D",   name, ".rds"))
    true_r_ld    <- readRDS(paste0(path, "dataset/true_LD_r",   name, ".rds"))
    true_Q_ld    <- readRDS(paste0(path, "dataset/true_LD_Q",   name, ".rds"))

    # Plots parameters
    legende1 <- c("D'", expression(r^2), expression(Q^'*'))

    # Position of legend
    pos <- NULL #c(0.20, 0.70)

    for(l in 1:n_Sim_Loci){ # 2 or 5 loci
      # Building the prevalence dataframe
      df_D_ld   <- dataframe_builder_LD(D_ld, 'D', l, true_D_ld)
      df_r_ld   <- dataframe_builder_LD(r_ld, 'r',l, true_r_ld)
      df_Q_ld   <- dataframe_builder_LD(Q_ld, 'Q',l, true_Q_ld)

      df <- rbind(df_D_ld, df_r_ld, df_Q_ld)

      for(k in 1:n_Freq_Distr){  # sym or asym
        for(j in samp_Vec){
          df1 <- df %>%
                filter(sample == j, shape == shape_typ[k]) %>%
                droplevels()

          df1$psiVec <- psi(lbda_Vec)
          
          p <- ggplot(data = df1, aes(x=psiVec))
          p <- p + geom_line(aes(y = df1[,'ld'], color = type, linetype = vers), size=1.)
          p <- p + expand_limits(y = 0)
          p <- p + scale_y_continuous(limits = c(0, 1))
          p <- beautify(p, legende1, legende2, pos, cbPalette, lty, NULL, NULL)
          if(name == 'Kenya'){
              p <- p + labs(x=expression(frac(lambda, 1 - e^-lambda)), y="", title=paste0("N = ", j, ", year = ", estim_Years[k]))
              outfile <- paste0(path,"plots/LD_plots_", dir, "/ld_SSize_", j, "_year_", estim_Years[k], "_", name, ".pdf")
          }else{
              p <- p + labs(x=expression(frac(lambda, 1 - e^-lambda)), y="", title=paste0("N = ", j, ", m = ", gen[l,1], ", n = ", gen[l,2]))
              outfile <- paste0(path,"plots/LD_plots_", dir, "/ld_", shape_typ[k], "_SSize_", j, "_nloci_", gen[l,1],"_",gen[l,2], "_", name, ".pdf")
          }

          pdf(outfile, height=5, width=8)
          print(p)
          dev.off()
        }
      }
    }
  }

  if(1==1){ # Plotting LD vs Adhoc LD
    # Importing the data to plot
    D_ld      <- readRDS(paste0(path, "dataset/estim_LD_D",   name, ".rds"))
    r_ld      <- readRDS(paste0(path, "dataset/estim_LD_r",   name, ".rds"))
    Q_ld      <- readRDS(paste0(path, "dataset/estim_LD_Q",   name, ".rds"))

    adhoc_D_ld1  <- list(readRDS(paste0(path, "dataset/adhocModelDpEstimates1_1",  name, ".rds")), readRDS(paste0(path, "dataset/adhocModelDpEstimates1_2",  name, ".rds")), readRDS(paste0(path, "dataset/adhocModelDpEstimates1_3",  name, ".rds")))
    adhoc_r_ld1  <- list(readRDS(paste0(path, "dataset/adhocModelrEstimates1_1",   name, ".rds")), readRDS(paste0(path, "dataset/adhocModelrEstimates1_2",   name, ".rds")), readRDS(paste0(path, "dataset/adhocModelrEstimates1_3",   name, ".rds")))
    adhoc_Q_ld1  <- list(readRDS(paste0(path, "dataset/adhocModelQEstimates1_1",   name, ".rds")), readRDS(paste0(path, "dataset/adhocModelQEstimates1_2",   name, ".rds")), readRDS(paste0(path, "dataset/adhocModelQEstimates1_3",   name, ".rds")))

    adhoc_D_ld2  <- list(readRDS(paste0(path, "dataset/adhocModelDpEstimates2_1",  name, ".rds")), readRDS(paste0(path, "dataset/adhocModelDpEstimates2_2",  name, ".rds")), readRDS(paste0(path, "dataset/adhocModelDpEstimates2_3",  name, ".rds")))
    adhoc_r_ld2  <- list(readRDS(paste0(path, "dataset/adhocModelrEstimates2_1",   name, ".rds")), readRDS(paste0(path, "dataset/adhocModelrEstimates2_2",   name, ".rds")), readRDS(paste0(path, "dataset/adhocModelrEstimates2_3",   name, ".rds")))
    adhoc_Q_ld2  <- list(readRDS(paste0(path, "dataset/adhocModelQEstimates2_1",   name, ".rds")), readRDS(paste0(path, "dataset/adhocModelQEstimates2_2",   name, ".rds")), readRDS(paste0(path, "dataset/adhocModelQEstimates2_3",   name, ".rds")))

    true_D_ld     <- readRDS(paste0(path, "dataset/true_LD_D",  name, ".rds"))
    true_r_ld    <- readRDS(paste0(path, "dataset/true_LD_r",   name, ".rds"))
    true_Q_ld    <- readRDS(paste0(path, "dataset/true_LD_Q",   name, ".rds"))

    # Plots parameters
    legendeD <- c("D'", expression(D[w]^"'"), expression(D[s]^"'"))
    legendeR <- c( expression(paste(r^2)), expression(r[w]^2), expression(r[s]^2))
    legendeQ <- c( expression(Q^"*"), expression(Q[w]^"*"), expression(Q[s]^"*"))

    # Position of legend
    pos <- NULL # c(0.15, 0.55)

    for(l in 1:n_Sim_Loci){ # 2 or 5 loci
      # Building the prevalence dataframe
      df_D_ld   <- dataframe_builder_LD(D_ld, 'D', l, true_D_ld)
      df_r_ld   <- dataframe_builder_LD(r_ld, 'r', l, true_r_ld)
      df_Q_ld   <- dataframe_builder_LD(Q_ld, 'Q', l, true_Q_ld)

      df_D_ld1   <- dataframe_builder_LD(adhoc_D_ld1, 'adhD1', l, true_D_ld)
      df_D_ld2   <- dataframe_builder_LD(adhoc_D_ld2, 'adhD2', l, true_D_ld)

      df_r_ld1   <- dataframe_builder_LD(adhoc_r_ld1, 'adhr1', l, true_r_ld)
      df_r_ld2   <- dataframe_builder_LD(adhoc_r_ld2, 'adhr2', l, true_r_ld)

      df_Q_ld1   <- dataframe_builder_LD(adhoc_Q_ld1, 'adhQ1', l, true_Q_ld)
      df_Q_ld2   <- dataframe_builder_LD(adhoc_Q_ld2, 'adhQ2', l, true_Q_ld)

      df <- rbind(df_D_ld, df_D_ld1, df_D_ld2, df_r_ld, df_r_ld1, df_r_ld2, df_Q_ld, df_Q_ld1, df_Q_ld2)

      for(k in 1:n_Freq_Distr){  # sym or asym
        for(j in samp_Vec){
          for (m in c('D', 'r', 'Q')){
            legende1 <- legendeD
            if(m == 'r'){
              legende1 <- legendeR
            }
            if (m == 'Q') {
               legende1 <- legendeQ
            }
           
            df11 <- df %>%
                  filter(sample == j, shape == shape_typ[k], type %in% c(m, paste0('adh',m, '1'), paste0('adh',m, '2')), vers == 'estimate') %>%
                  droplevels()

            df12 <- df %>%
                  filter(sample == j, shape == shape_typ[k], type == m, vers == 'true') %>%
                  droplevels()

            df1 <- rbind(df11, df12)

            df1$psiVec <- psi(lbda_Vec)
            
            p <- ggplot(data = df1, aes(x=psiVec))
            p <- p + geom_line(aes(y = df1[,'ld'], color = type, linetype = vers), size=1.)
            p <- p + expand_limits(y = 0)
            p <- p + scale_y_continuous(limits = c(0, 1))
            p <- beautify(p, legende1, legende2, pos, cbPalette, lty2, NULL, NULL)
            if(name == 'Kenya'){
                p <- p + labs(x=expression(frac(lambda, 1 - e^-lambda)), y="", title=paste0("N = ", j, ", year = ", estim_Years[k]))
                outfile <- paste0(path,"plots/LD_plots_", dir, "/ldvsADH_SSize_",m,"_", j, "_year_", estim_Years[k], "_", name, ".pdf")
            }else{
                p <- p + labs(x=expression(frac(lambda, 1 - e^-lambda)), y="", title=paste0("N = ", j, ", m = ", gen[l,1], ", n = ", gen[l,2]))
                outfile <- paste0(path,"plots/LD_plots_", dir, "/ldvsADH_", shape_typ[k], "_SSize_",m,"_", j, "_nloci_", gen[l,1],"_",gen[l,2], "_", name, ".pdf")
            }

            pdf(outfile, height=5, width=8)
            print(p)
            dev.off()
          }
        }
      }
    }
  }

  legende1  <- samp_Vec
  
  if(1==1){ # Plotting bias for haplotype frequencies
    # Importing the data to plot
    freqbias <- readRDS(paste0(path, "dataset/freqbias", name, ".rds"))

      for(l in 1:n_Sim_Loci){ # 2 or 5 loci
        # Building the frequencies bias dataframe
        df <- dataframe_builder_Freqperf(freqbias, l)
        tru_freq <- sim_Param[[1]][[l]]

        for(k in 1:n_Freq_Distr){  # sym or asym
          trufreq_vec <- tru_freq[k,]
          for(i in 1:n_Hapl [l]){
            if(trufreq_vec[i] != 0){
              for(j in samp_Vec){
                df1 <- df %>%
                    filter(freq == i, shape == shape_typ[k]) %>%
                    droplevels()
                df1$psiVec <- psi(lbda_Vec)

                p <- ggplot(data = df1, aes(x=psiVec))
                p <- p + geom_line(aes(y = df1[,'bias'], color = sample), size=1.)
                p <- beautify(p, legende1, legende2, NULL, cbPalette, lty, 'N', 'italic')
                p <- p + expand_limits(y=0)
                p <- p + labs(x=expression(frac(lambda, 1 - e^-lambda)), y=paste0('bias frequencies (in %)'), title = paste0('p = ', round(trufreq_vec[i], 3)))
                if(name == 'Kenya'){
                  outfile <- paste0(path,"plots/Bias_cv_plots_", dir, "/frequency/", "bias_freq_", i, "_year_", estim_Years[k], "_", name, ".pdf")
                }else{
                  outfile <- paste0(path,"plots/Bias_cv_plots_", dir, "/frequency/", "bias_", shape_typ[k], "_freq_", i, "_nloci_", gen[l,1],"_",gen[l,2], "_",name, ".pdf")
                }
                pdf(outfile, height=5, width=8)
                print(p)
                dev.off()
              }
            }
          }
        }
    }
  }

  if(1==1){ # Plotting bias and coefficient of variation for MOI
    # Importing the data to plot
    moibias  <- readRDS(paste0(path, "dataset/moibias", name, ".rds"))
    moicv    <- readRDS(paste0(path, "dataset/moicv", name, ".rds"))

      for(l in 1:n_Sim_Loci){ # 2 or 5 loci
        # Building the prevalence dataframe
        df_cv   <- dataframe_builder_Moiperf(moicv, l)
        df_bias <- dataframe_builder_Moiperf(moibias, l)

        tru_freq <- sim_Param[[1]][[l]]

        for(k in 1:n_Freq_Distr){  # sym or asym
          trufreq_vec <- tru_freq[k,]
          df1 <- df_bias %>%
                filter(shape == shape_typ[k]) %>%
                droplevels()
          df1$psiVec <- psi(lbda_Vec)

          p <- ggplot(data = df1, aes(x=psiVec))
          p <- p + geom_line(aes(y = df1[,'bias'], color = sample), size=1.)
          p <- beautify(p, legende1, legende2, NULL, cbPalette, lty, 'N', 'italic')
          p <- p + expand_limits(y=0)
          p <- p + labs(x=expression(frac(lambda, 1 - e^-lambda)), y=paste0('bias MOI (in %)'))
          if(name == 'Kenya'){
              outfile <- paste0(path,"plots/Bias_cv_plots_", dir, "/moi/", "bias_moi_year_", estim_Years[k], "_", name, ".pdf")
          }else{
              outfile <- paste0(path,"plots/Bias_cv_plots_", dir, "/moi/", "bias_", shape_typ[k], "_moi_nloci_", gen[l,1],"_",gen[l,2], "_",name, ".pdf")
          }
          pdf(outfile, height=5, width=8)
          print(p)
          dev.off()

          df2 <- df_cv %>%
                filter(shape == shape_typ[k]) %>%
                droplevels()
          df2$lbd <- lbda_Vec

          p <- ggplot(data = df2, aes(x=lbd))
          p <- p + geom_line(aes(y = df2[,'bias'], color = sample), size=1.)
          p <- beautify(p, legende1, legende2, NULL, cbPalette, lty, 'N', 'italic')
          p <- p + expand_limits(y=0)
          p <- p + labs(x=expression(frac(lambda, 1 - e^-lambda)), y=paste0('Coef. variation MOI'))
           if(name == 'Kenya'){
              outfile <- paste0(path,"plots/Bias_cv_plots_", dir, "/moi/", "coefvar_moi_year_", estim_Years[k], "_", name, ".pdf")
          }else{
              outfile <- paste0(path,"plots/Bias_cv_plots_", dir, "/moi/", "coefvar_", shape_typ[k], "_moi_nloci_", gen[l,1],"_",gen[l,2], "_",name, ".pdf")
          }
          pdf(outfile, height=5, width=8)
          print(p)
          dev.off()
        }
      }
  }

  if(1==1){ # Plotting bias and standard deviation LD
    # Importing the data to plot
    ldtype <- c('D', 'r', 'Q')
    ldlabel <- c("D'", paste(expression(r^2)), "Q*")
    idx <- c('w', 's')

    # Plots parameters
    legendeD <- c("D'", expression(D[w]^"'"), expression(D[s]^"'"))
    legendeR <- c( expression(paste(r^2)), expression(r[w]^2), expression(r[s]^2))
    legendeQ <- c( expression(Q^"*"), expression(Q[w]^"*"), expression(Q[s]^"*"))
    
    legende <- list(legendeD, legendeR, legendeQ)
    for(i in 1:3){
      LDcv    <- readRDS(paste0(path, "dataset/LDcv_", paste0('',ldtype[i]), "_", 3, name, ".rds"))
      LDcv1   <- readRDS(paste0(path, "dataset/LDcv_", paste0('adhoc',ldtype[i]), "_", 1, name, ".rds"))
      LDcv2   <- readRDS(paste0(path, "dataset/LDcv_", paste0('adhoc',ldtype[i]), "_", 2, name, ".rds"))

      if(i == 1){
        legende2 <- legendeD
      }else if (i == 2) {
          legende2 <- legendeR
      }else{
        legende2 <- legendeQ
      }

      for(l in 1:n_Sim_Loci){ # 2 or 5 loci
        # Building the prevalence dataframe
        df_cv   <- dataframe_builder_Moiperf(LDcv, l)
        df_cv$vers <- ldtype[i]
        df_cv1   <- dataframe_builder_Moiperf(LDcv1, l)
        df_cv1$vers <- paste(ldtype[i],"w")
        df_cv2   <- dataframe_builder_Moiperf(LDcv2, l)
        df_cv2$vers <- paste(ldtype[i],"s")

        df <- rbind(df_cv, df_cv1, df_cv2)

        for(k in 1:n_Freq_Distr){  # sym or asym
          df1 <- df %>%
            filter(shape == shape_typ[k]) %>%
            droplevels()
          df1$lbd <- lbda_Vec

          p <- ggplot(data = df1, aes(x=lbd))
          p <- p + geom_line(aes(y = df1[,'bias'], color = sample, linetype = vers), size=1.)
          p <- beautify(p, legende1, legende2, NULL, cbPalette, lty2, 'N', 'italic')
          p <- p + expand_limits(y=0)
          p <- p + labs(x=expression(frac(lambda, 1 - e^-lambda)), y=paste('Variance ', ldlabel[i]), title=paste0("m = ", gen[l,1], ", n = ", gen[l,2]))
          if(name == 'Kenya'){
              outfile <- paste0(path,"plots/Bias_cv_plots_", dir, "/ld/", "coefvar_moi_year_", estim_Years[k], "_", name, ".pdf")
          }else{
              outfile <- paste0(path,"plots/Bias_cv_plots_", dir, "/ld/", "coefvar_ldvsadhoc_", ldtype[i], "_", shape_typ[k], "_moi_nloci_", gen[l,1],"_",gen[l,2], "_",name, ".pdf")
          }
          pdf(outfile, height=5, width=8)
          print(p)
          dev.off()
        }
      }
    }    
  }
}

# Relative path
path <- "/Volumes/GoogleDrive-117934057836063832284/My Drive/Maths against Malaria/Christian/Models/STRModel/"

# Define data origin ('' <- simualted data, kenya <- kenyan data)
namelist <- c('') #c('', 'Kenya')

for (name in namelist){
  
  # Loading true haplotype parameters for the simulation
  sim_Param <- readRDS(paste0(path, "dataset/true_Parameters", name, ".rds"))

  # Loading extra parameters
  extra_Sim_Param  <- readRDS(paste0(path, "dataset/extra_Parameters", name, ".rds"))

  # Simulation parameters
  n_Lbda        <- extra_Sim_Param[[1]]
  n_Sim_Loci    <- extra_Sim_Param[[2]]
  n_Hapl        <- extra_Sim_Param[[3]]
  n_Freq_Distr  <- extra_Sim_Param[[6]]
  gen           <- extra_Sim_Param[[7]]

  lbda_Vec      <- sim_Param[[2]]
  samp_Vec      <- sim_Param[[3]]

  estim_Years   <- c(2005, 2010)

  # Plotting
  main(sim_Param, name, gen)
}
