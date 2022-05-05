# Importing libraries
library(ggplot2)

# Importing the simulated dataset
Fdat <- Model2_simulated_data

# Defining the color palette (color-blind friendly) for the plots
cbPalette <-c(rgb(0,.45,.7),rgb(.9,.6,0),rgb(0,0,0),rgb(0,.6,.5), rgb(.8,.4,0))

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
    geom_line(aes(x = FDatsubs$moi, y = FDatsubs$Bias_freq, colour = factor(FDatsubs$Sample)))+
    labs(x = "MOI", y = "Relative bias Frequency in %", title = titre)+
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
  
  ggsave(plt.p1, file =paste0("Relative_Bias_Freq_", k, "_.pdf"))
  
  # Plotting and saving the relative bias for moi
  plt.p1 <- ggplot(FDatsubs) +
    geom_line(aes(x = FDatsubs$moi, y = FDatsubs$Bias_moi, colour = factor(FDatsubs$Sample)))+
    labs(x = "MOI", y = "Relative bias MOI in %", title = titre)+
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
  
  ggsave(plt.p1, file =paste0("Relative_Bias_MOI_", k, "_.pdf"))
  
  # Plotting and saving the standard deviation
  plt.s <- ggplot(FDatsubs) +
    geom_line(aes(x = FDatsubs$moi, y = FDatsubs$Variation_coef_moi, colour = factor(FDatsubs$Sample)))+
    labs(x = "MOI", y = "Coefficient of variation MOI in %", title = titre)+
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
  
  
  ggsave(plt.s, file =paste0("Coefficient_Variation_Moi_", k, "_.pdf"))
}

