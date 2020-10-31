# Script for plotting Figures 7 and 8 from Krantsevich et al., 
# "BAYESIAN DECISION THEORY FOR TREE-BASED ADAPTIVE SCREENING TESTS 
#     WITH AN APPLICATION TO YOUTH DELINQUENCY"

library(rpart.plot)

# Set parameters
maxIPP = 3

# Load fitted CART tree for subpopulation of Age >= 15
folder = "data/in_sample/subpopulation/model_fits/"
savefile = paste0(folder, "fit.CART.XB.regression.maxIPP.", maxIPP, ".opt")
load(savefile) # fitted tree was called "opt.tree"

# Plot the tree. Option 'tweak' determines the text size as a proportion
# above/below default; e.g. tweak=1.2 is 120% of default text size.
# Change it to suit your plotting window / tree size
rpart.plot(opt.tree, box.palette = "Grays", tweak=1.1)

# Load fitted CART tree and plot for population of all youth
folder = "data/in_sample/all/model_fits/"
savefile = paste0(folder, "fit.CART.XB.regression.maxIPP.", maxIPP, ".opt")
load(savefile)  
rpart.plot(opt.tree, box.palette = "Grays", tweak=1.1)

