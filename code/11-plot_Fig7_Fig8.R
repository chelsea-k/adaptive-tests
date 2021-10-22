# Script for plotting Figures 7 and 8 from Krantsevich et al., 
# "BAYESIAN DECISION THEORY FOR TREE-BASED ADAPTIVE SCREENING TESTS 
#     WITH AN APPLICATION TO YOUTH DELINQUENCY"

library(rpart.plot)

# Set parameters
maxIPP <- 3

# Load fitted CART tree for subpopulation of Age >= 15
folder <- "output/in_sample/subpopulation/model_fits_CART"
savefile <- file.path(folder, paste0("fit.CART.XB.regression_maxIPP.", maxIPP, ".opt"))
load(savefile) # fitted tree was called "opt_tree"
opt_tree_sub <- opt_tree
# Plot the tree. Option 'tweak' determines the text size as a proportion
# above/below default; e.g. tweak=1.2 is 120% of default text size.
# Change it to suit your plotting window / tree size
rpart.plot(opt_tree_sub, tweak=1.3, roundint = FALSE)

# Load fitted CART tree and plot for population of all youth
folder <- "output/in_sample/all/model_fits_CART"
savefile <- file.path(folder, paste0("fit.CART.XB.regression_maxIPP.", maxIPP, ".opt"))
load(savefile)  
rpart.plot(opt_tree, tweak=1.1, roundint = FALSE)



