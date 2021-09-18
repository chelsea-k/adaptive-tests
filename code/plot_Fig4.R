# Script for plotting Figure 4 from Krantsevich et al., 
# "BAYESIAN DECISION THEORY FOR TREE-BASED ADAPTIVE SCREENING TESTS 
#     WITH AN APPLICATION TO YOUTH DELINQUENCY"

source("code/util_functions.R")  # get_cutoff(), get_utility()
library(ggplot2)

# Set parameters
maxIPP.vals = c(2:15)
n.maxIPP = length(maxIPP.vals)
w.vals = c(0.5, 0.65, 0.8)

# Set up data folders and load data
data.dir = "data/in_sample/all/synthetic_data/"
model.dir = "data/in_sample/all/model_fits/"
results.dir = "data/in_sample/all/results/"

synth.df = read.csv(paste0(data.dir, "synth.XB.csv"))
p.CART.df = read.csv(paste0(results.dir, "p.maxIPP.XB.synth.csv")) 
n.mcmc = length(unique(synth.df$post.idx))
n.rows.total = n.maxIPP*n.mcmc

# Dataframe for storing results
delta.draws.all = data.frame(matrix(NA, nrow=3*n.rows.total, ncol=4))
colnames(delta.draws.all) = c("maxIPP", "post.idx", "Delta", "w")
cutoffs = data.frame(matrix(NA, nrow=length(w.vals), ncol=2+n.maxIPP))
colnames(cutoffs)=c("w", "XBART", paste0("CART.m.",maxIPP.vals))

# Compute delta results for all values of w and maxIPP
for (k in seq_along(w.vals)){
  cat(paste0("-------------- Predicting for cutoff #", k, " -------------------\n"))
  w = w.vals[[k]]  
  
  # compute optimal cutoffs
  XB.cutoff = get_cutoff(as.factor(synth.df$y), synth.df$phat, w)
  maxIPP.cutoff = rep(NA, n.maxIPP)
  for (i in seq_along(maxIPP.vals)){
    maxIPP = maxIPP.vals[[i]]
    col = which(colnames(p.CART.df)==paste0("m.", maxIPP))
    maxIPP.cutoff[[i]] = get_cutoff(as.factor(synth.df$y), p.CART.df[,col], w)
  }
  cutoffs[k,] = c(w, XB.cutoff, maxIPP.cutoff)
  
  # compute delta draws
  delta.draws = data.frame(matrix(NA, nrow=n.rows.total, ncol=3))
  colnames(delta.draws) = c("maxIPP", "post.idx", "Delta")
  
  for (j in 1:n.mcmc){
    if (j%%100==0) {cat(paste0("On iteration j = ", j, " out of ", n.mcmc, "\n"))}
    post.rows = which(synth.df$post.idx==j)
    temp.df = synth.df[post.rows,]
    XB.metrics = get_utility(temp.df$y,prob=temp.df$phat,cut=XB.cutoff,type="prob",w=w) 
    for (i in seq_along(maxIPP.vals)){
      temp.df2 = p.CART.df[post.rows,]
      maxIPP = maxIPP.vals[[i]]
      col = which(colnames(p.CART.df)==paste0("m.", maxIPP))
      CT.metrics = get_utility(temp.df$y, prob=temp.df2[,col], 
                               cut=maxIPP.cutoff[[i]], type="prob", w=w) 
      delta.draws[n.maxIPP*(j-1)+i,] = data.frame(maxIPP=maxIPP, 
                                                   post.idx=j, 
                                                   Delta=CT.metrics$util-XB.metrics$util)
    }
  }
  delta.draws.all[(k-1)*n.rows.total+(1:n.rows.total),] = cbind(delta.draws, w=w)
}
write.csv(delta.draws.all, paste0(results.dir, "delta.draws.all.csv"), row.names = FALSE)
write.csv(cutoffs, paste0(results.dir, "cutoffs.csv"), row.names = FALSE)

# Post-process results and create plot for Figure 4
delta.draws.all$maxIPP = as.factor(delta.draws.all$maxIPP)
delta.draws.all$w = as.factor(delta.draws.all$w)
cbPalette = c("#009E73", "#0072B2", "#D55E00")

ggplot(delta.draws.all, aes(x=maxIPP, y=Delta, fill=w)) +
 geom_boxplot(outlier.shape=NA, position=position_dodge2(padding=0.2, width=1.2)) + 
 theme(plot.title = element_text(hjust = 0.5))  + 
 ggtitle(paste0("Utility Differences")) + ylab(expression(Delta)) +
 scale_fill_manual(values=cbPalette) +
 theme_bw()


