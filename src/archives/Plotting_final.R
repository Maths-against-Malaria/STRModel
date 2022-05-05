# Importing libraries
library(ggplot2)

# Importing the simulated dataset
Fdat <-read.csv("Dataset/Model2_simulated_data.csv")

# Defining the color palette (color-blind friendly) for the plots
cbPalette <-c(rgb(0 ,0 ,0 ,maxColorValue = 255), rgb(230 ,159 ,0 ,maxColorValue = 255),rgb(204 ,121 ,167 ,maxColorValue = 255),rgb(0 ,114 ,178 ,maxColorValue = 255),rgb(0 ,158 ,115 ,maxColorValue = 255), rgb(213 ,94 ,0 ,maxColorValue = 255))

# Plotting for each set of frequencies
for (k in 1:2){
  FDatsubs <- Fdat[Fdat$frequency == k,]
  
  if (k == 1){
    titre = "Symmetric distribution"
  }
  else{
    titre = "Skewed distribution"
  }
  
  # Plotting and saving the relative bias for frequencies
  plt.p1 <- ggplot(FDatsubs) +
    geom_line(aes(x = FDatsubs$lambda, y = FDatsubs$Bias_freq, colour = factor(FDatsubs$Sample)))+
    labs(x = expression(lambda), y = "Relative bias of Frequency in %", title = titre)+
    scale_colour_manual(values=cbPalette ,guide = guide_legend(title = "Sample size"))+
    theme(plot.title = element_text(hjust = 0.5,size=rel(1.5)),
          legend.text = element_text(size = rel(1.5)),
          legend.title = element_text(size = rel(1.5),face="plain"),
          legend.position = "right",
          aspect.ratio = 9/16,
          panel.background = element_rect(colour= NA, fill='white'),
          panel.border = element_rect(fill=NA),
          axis.text = element_text(size = rel(1.5)),axis.title = element_text(size = rel(1.5)),
          legend.key=element_blank(), legend.text.align = 0,legend.background = element_blank())
  
  ggsave(plt.p1, file =paste0("Plots/Relative_Bias_Freq_", k, "_.pdf"))
  
  # Plotting and saving the relative bias for moi
  plt.p1 <- ggplot(FDatsubs) +
    geom_line(aes(x = FDatsubs$lambda, y = FDatsubs$Bias_moi, colour = factor(FDatsubs$Sample)))+
    labs(x = expression(lambda), y = "Relative bias of MOI in %", title = titre)+
    scale_colour_manual(values=cbPalette ,guide = guide_legend(title = "N"))+
    theme(plot.title = element_text(hjust = 0.5,size=rel(1.5)),title= element_text(face ="italic"), 
          legend.text = element_text(size = rel(1.5)),
          legend.title = element_text(size = rel(1.5),face="plain"),
          legend.position = "right",
          aspect.ratio = 9/16,
          panel.background = element_rect(colour='black',fill='white'),
          panel.border = element_rect(fill=NA),
          axis.text = element_text(size = rel(1.5)),axis.title = element_text(size = rel(1.5)),
          legend.key=element_blank(), legend.text.align = 0,legend.background = element_blank())
  
  ggsave(plt.p1, file =paste0("Plots/Relative_Bias_MOI_", k, "_.pdf"))
  
  # Plotting and saving the standard deviation
  plt.s <- ggplot(FDatsubs) +
    geom_line(aes(x = FDatsubs$lambda, y = FDatsubs$Variation_coef_moi, colour = factor(FDatsubs$Sample)))+
    labs(x = expression(lambda), y = "Coefficient of variation MOI in %", title = titre)+
    scale_colour_manual(values=cbPalette ,guide = guide_legend(title = "N"))+
    theme(plot.title = element_text(hjust = 0.5,size=rel(1.5)),title= element_text(face ="italic"), 
          legend.text = element_text(size = rel(1.5)),
          legend.title = element_text(size = rel(1.5),face="plain"),
          legend.position = "right",
          aspect.ratio = 9/16,
          panel.background = element_rect(colour='black',fill='white'),
          panel.border = element_rect(fill=NA),
          axis.text = element_text(size = rel(1.5)),axis.title = element_text(size = rel(1.5)),
          legend.key=element_blank(), legend.text.align = 0,legend.background = element_blank())
  
  
  ggsave(plt.s, file =paste0("Plots/Coefficient_Variation_Moi_", k, "_.pdf"))
}

