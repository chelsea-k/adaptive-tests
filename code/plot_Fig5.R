# Script for plotting Figure 5 from Krantsevich et al., 
# "BAYESIAN DECISION THEORY FOR TREE-BASED ADAPTIVE SCREENING TESTS 
#     WITH AN APPLICATION TO YOUTH DELINQUENCY"

source("code/util_functions.R")  # get_metrics()
library(ggplot2)
library(caret)

# Set parameters
maxIPP = 3
w.list = c(0.5, 0.65, 0.8)
post.idx.list = c(10,20, 30,40)

# Set up data folders and load item response data, cutoff dataframe
data.dir = "data/in_sample/all/synthetic_data/"
results.dir = "data/in_sample/all/results/"
synth.df = read.csv(paste0(data.dir, "synth.XB.csv"))
p.CART.df = read.csv(paste0(results.dir, "p.maxIPP.XB.synth.csv")) 
cutoffs = read.csv(paste0(results.dir, "cutoffs.csv"))

# Compute sensitivity/specificity for XBART, CART ROC curves
#   and (spec, sens) pairs for cutoff-specific points 

metrics.XBART = data.frame()
metrics.CART = data.frame()
ROC.points = data.frame()

for (t in seq_along(post.idx.list)) {
  post.idx = post.idx.list[[t]]
  ROC.samp = which(synth.df$post.idx == post.idx)
  col = which(colnames(p.CART.df)==paste0("m.",maxIPP))
  p.CART = p.CART.df[ROC.samp,col]
  y.class = synth.df[ROC.samp, "y"]
  p.XBART = synth.df[ROC.samp, "phat"]
  
  metrics.XBART.temp = get_metrics(p.XBART, y.class)
  metrics.XBART = rbind(metrics.XBART, cbind(group=0, post.idx = post.idx, 
                              Model="XBART", metrics.XBART.temp))
  metrics.CART.temp = get_metrics(p.CART, y.class)
  metrics.CART = rbind(metrics.CART, cbind(group=0, post.idx = post.idx, 
                              Model="maxIPP=6", metrics.CART.temp))
  
  ROC.points.temp = data.frame(matrix(0, nrow=2*length(w.list), ncol=5))
  colnames(ROC.points.temp) = c("Specificity", "Sensitivity", "w", "post.idx", "Model")
  for (k in 1:length(w.list)){
    w = w.list[[k]]
    ct.row = which(cutoffs$w==w)
    cm.XBART = confusionMatrix(as.factor(1*(p.XBART>=cutoffs[ct.row, "XBART"])), 
                               as.factor(y.class), positive="1")
    ROC.points.temp[2*k-1,] = list(cm.XBART$byClass[['Specificity']], 
                          cm.XBART$byClass[['Sensitivity']], w, post.idx, "XBART")
    CART.col = which(colnames(cutoffs)==paste0("CART.m.", maxIPP))
    cm.CART = confusionMatrix(as.factor(1*(p.CART>=cutoffs[ct.row,CART.col])), 
                               as.factor(y.class), positive="1")
    ROC.points.temp[2*k,] = list(cm.CART$byClass[['Specificity']], 
                          cm.CART$byClass[['Sensitivity']], w, post.idx, "CART")
  }
  ROC.points = rbind(ROC.points, ROC.points.temp)
}

# Post-process results and create plot for Figure 5
ROC.points$w = as.factor(ROC.points$w)
ROC.points$post.idx = as.integer(ROC.points$post.idx)
cbPalette = c("#009E73", "#0072B2", "#D55E00")

ggplot(ROC.points, aes(x=Specificity, y=Sensitivity, group=w)) +
  xlim(c(1,0)) +
  geom_line(data=metrics.XBART, color="black", 
            aes(x=Specificity, y=Sensitivity, linetype="XBART")) +
  geom_line(data=metrics.CART, color="black", 
            aes(x=Specificity, y=Sensitivity, linetype=paste0("maxIPP=", maxIPP,""))) +
  geom_point(size=3.5,  aes(color=w, shape=w)) +
  scale_linetype_manual(name="Action", values=c("dashed","solid")) +
  scale_colour_manual(values = cbPalette) + 
  scale_shape_manual(values=c(16,15,17)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_grid(cols=vars(post.idx), labeller = label_both) +
  theme_bw()




