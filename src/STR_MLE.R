# Title        : Template script for the MLE of haplotype frequencies, MOI, and prevalence
#                from example dataset.
# Objective    : Estimate haplotype frequencies, MOI, and prevalence
#                from genomic/molecular data (microsatellite data)
# Created by   : Christian Tsoungui Obama
# Created on   : 05.05.22
# Last modified: 05.05.22

# Install the necessary packages if necessary
install.packages('xlsx')   # Comment this line if xlsx installed

# Loading libraries
library(xlsx)

# Import the dataset
DATA <- read.xlsx('/home/janedoe/Documents/example.xlsx', 1, header = TRUE)

# Load external resources
source("/home/janedoe/Documents/STRModel.R")

# Find the MLEs
est <- mle(DATA, id=TRUE)
