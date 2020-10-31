# Script for plotting Figure 10 from Krantsevich et al., 
# "BAYESIAN DECISION THEORY FOR TREE-BASED ADAPTIVE SCREENING TESTS 
#     WITH AN APPLICATION TO YOUTH DELINQUENCY"

source('code/util_functions.R')
library(ggplot2)

# Set parameters
maxIPP.vals = c(2:15)
n.maxIPP = length(maxIPP.vals)
w = 0.5

# Set up folders and load data
data.dir = "data/out_of_sample/all/synthetic_data/"
results.dir = "data/out_of_sample/all/results/"

data.test = read.csv("data/item_response_data/item.response.data.test.csv")
  
synth.df = read.csv(paste0(data.dir, "synth.XB.csv"))
n.mcmc = length(unique(synth.df$post.idx))
n.rows.total = n.maxIPP*n.mcmc

p.maxIPP.XB.synth = read.csv(paste0(results.dir, "p.maxIPP.XB.synth.csv"))
p.maxIPP.XB.test = read.csv(paste0(results.dir, "p.maxIPP.XB.test.csv"))

load(paste0(data.dir, "XB.predict"))
p.XBART.test = XB.predict$test


# Compute optimal cutoffs and predicted classes 
XBART.cutoff = get_cutoff(as.factor(synth.df$y), synth.df$phat, w)
maxIPP.cutoff = rep(NA, n.maxIPP)
for (i in seq_along(maxIPP.vals)){
  maxIPP = maxIPP.vals[[i]]
  col = which(colnames(p.maxIPP.XB.synth)==paste0("m.", maxIPP))
  maxIPP.cutoff[[i]] = get_cutoff(as.factor(synth.df$y),p.maxIPP.XB.synth[,col],w)
}

# Compute delta results for trainig data, for all values of maxIPP
delta.draws.train = data.frame(matrix(NA, nrow=n.rows.total, ncol=3))
colnames(delta.draws.train) = c("maxIPP", "post.idx", "Delta")
for (j in 1:n.mcmc){
  if (j%%100==0) {cat(paste0("On iteration j = ", j, " out of ", n.mcmc, "\n"))}
  post.idx.rows = which(synth.df$post.idx==j)
  temp.df = synth.df[post.idx.rows,]
  y.temp = temp.df$y
  XB.metrics = get_utility(y.temp,prob=temp.df$phat,cut=XBART.cutoff,type="prob",w=w)
  for (i in seq_along(maxIPP.vals)){
    temp.df2 = p.maxIPP.XB.synth[post.idx.rows,]
    maxIPP = maxIPP.vals[i]
    col = which(colnames(p.maxIPP.XB.synth)==paste0("p.CART.m.", maxIPP))
    CT.metrics = get_utility(y.temp,prob=temp.df2[,col],cut=maxIPP.cutoff[[i]],type="prob",w=w)
    delta.draws.train[n.maxIPP*(j-1)+i,] = data.frame(maxIPP=maxIPP, 
                                                      post.idx=j, 
                                                      Delta=CT.metrics$util-XB.metrics$util)
  }
}
write.csv(delta.draws.train,paste0(results.dir,"delta.draws.train.",w,".csv"),row.names=F)

# Compute delta results for testing data, for all values of maxIPP
delta.draws.test = data.frame(matrix(NA, nrow=n.maxIPP, ncol=3))
colnames(delta.draws.test) = c("maxIPP", "post.idx", "Delta")
XB.met.test = get_utility(data.test$y,prob=p.XBART.test,cut=XBART.cutoff,type="prob",w=w)
for (i in seq_along(maxIPP.vals)){
  maxIPP = maxIPP.vals[[i]]
  col = which(colnames(p.maxIPP.XB.test)==paste0("p.CART.m.", maxIPP))
  CT.met.test = get_utility(data.test$y,prob=p.maxIPP.XB.test[,col],
                           cut=maxIPP.cutoff[[i]],type="prob",w=w)
  delta.draws.test[i,] = data.frame(maxIPP=maxIPP, 
                                           post.idx=0, 
                                           Delta=CT.met.test$util-XB.met.test$util)
}
write.csv(delta.draws.test,paste0(results.dir,"delta.draws.test.",w,".csv"),row.names=F)

# Post-process results and create plot for Figure 10
delta.draws.all = rbind(cbind(delta.draws.train, Type = "Predicted"), 
                        cbind(delta.draws.test, Type = "Actual"))

delta.draws.all$maxIPP = as.factor(delta.draws.all$maxIPP)
delta.draws.all$Type = factor(delta.draws.all$Type, levels=c("Predicted", "Actual"))
cbPalette = c("#000000","#D55E00", "#CC79A7")
ggplot(delta.draws.all, aes(x=maxIPP, y=Delta, color=Type)) +
  geom_boxplot(outlier.shape=NA, position="identity") + 
  theme(plot.title = element_text(hjust = 0.5))  + 
  ggtitle(paste0("Utility Differences")) + ylab(expression(Delta)) +
  scale_color_manual(values=cbPalette)


