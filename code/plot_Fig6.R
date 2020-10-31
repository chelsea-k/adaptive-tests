# Script for plotting Figure 6 from Krantsevich et al., 
# "BAYESIAN DECISION THEORY FOR TREE-BASED ADAPTIVE SCREENING TESTS 
#     WITH AN APPLICATION TO YOUTH DELINQUENCY"

source("code/util_functions.R")  # get_cutoff(), get_utility()
library(ggplot2)

# Set parameters
maxIPP.vals = c(2:15)
n.maxIPP = length(maxIPP.vals)
w = 0.5

# Set up data folders and load data
data.dir = "data/in_sample/subpopulation/synthetic_data/"
results.dir = "data/in_sample/subpopulation/results/"

synth.sub = read.csv(paste0(data.dir, "synth.XB.csv"))
p.CART.sub = read.csv(paste0(results.dir, "p.maxIPP.XB.synth.csv")) 

synth.all = read.csv("data/in_sample/all/synthetic_data/synth.XB.csv")
p.CART.all = read.csv("data/in_sample/all/results/p.maxIPP.XB.synth.csv")
n.mcmc = length(unique(synth.all$post.idx))
n.rows.total = n.maxIPP*n.mcmc

# Compute optimal cutoffs and predicted classes 
XB.cut.all = get_cutoff(as.factor(synth.all$y), synth.all$phat, w)
XB.cut.sub = get_cutoff(as.factor(synth.sub$y), synth.sub$phat, w)
maxIPP.cut.all = rep(NA, n.maxIPP)
maxIPP.cut.sub = rep(NA, n.maxIPP)
m.cols = which(colnames(p.CART.all) %in% paste0("m.", maxIPP.vals))
for (i in seq_along(maxIPP.vals)){
  maxIPP.cut.all[[i]] = get_cutoff(as.factor(synth.all$y), p.CART.all[,m.cols[[i]]], w)
  maxIPP.cut.sub[[i]] = get_cutoff(as.factor(synth.sub$y), p.CART.sub[,m.cols[[i]]], w)
}
y.XB.all = 1*(synth.all$phat >= XB.cut.all)
y.XB.sub = 1*(synth.sub$phat >= XB.cut.sub)
y.CART.all = get_class(p.CART.all, maxIPP.cut.all, maxIPP.vals, m.cols)
y.CART.sub = get_class(p.CART.sub, maxIPP.cut.sub, maxIPP.vals, m.cols)

# Compute delta results for all values of maxIPP for both populations
delta.draws.all = data.frame(matrix(NA, nrow=n.rows.total, ncol=4))
colnames(delta.draws.all) = c("maxIPP", "post.idx", "Delta", "Population")

for (j in 1:n.mcmc){
  if (j%%100==0) {cat(paste0("On iteration j = ", j, " out of ", n.mcmc, "\n"))}
  post.rows = which(synth.all$post.idx==j)
  temp.df.all = synth.all[post.rows,]
  temp.df.sub = synth.sub[post.rows,]
  XB.metrics.all = get_utility(temp.df.all$y,y.pred=y.XB.all[post.rows],type="class",w=w) 
  XB.metrics.sub = get_utility(temp.df.sub$y,y.pred=y.XB.sub[post.rows],type="class",w=w) 
  for (i in seq_along(maxIPP.vals)){
    maxIPP = maxIPP.vals[[i]]
    m.col.all = which(colnames(p.CART.all)==paste0("m.", maxIPP))
    m.col.sub = which(colnames(p.CART.sub)==paste0("m.", maxIPP))
    CART.metrics.all = get_utility(temp.df.all$y,y.pred=y.CART.all[post.rows,m.col.all],
                                   type="class",w=w) 
    CART.metrics.sub = get_utility(temp.df.sub$y,y.pred=y.CART.sub[post.rows,m.col.sub],
                                   type="class",w=w) 
    temp.results = data.frame(maxIPP = maxIPP, post.idx = j, 
                              Delta = c(CART.metrics.all$util - XB.metrics.all$util, 
                                        CART.metrics.sub$util - XB.metrics.sub$util),
                              Population=c("All", "Age >= 15"))
    delta.draws.all[c(1:2) + (2*(i-1)) + (2*n.maxIPP*(j-1)),] = temp.results
    }
}
write.csv(delta.draws.all, paste0(results.dir, "delta.draws.with.population.", w, ".csv"), row.names = FALSE)

# Post-process results and create plot for Figure 6
delta.draws.all$maxIPP = as.factor(delta.draws.all$maxIPP)
delta.draws.all$Population = factor(delta.draws.all$Population, levels=c("All", "Age >= 15"))
greyPalette = c("#A9A9A9","#000000")

ggplot(delta.draws.all, aes(x=maxIPP, y=Delta, fill=Population)) +
 geom_boxplot(outlier.shape=NA, position=position_dodge2(padding=0.2, width=10)) + 
 theme(plot.title = element_text(hjust = 0.5))  + 
 ggtitle(paste0("Utility Differences")) + 
 ylab(expression(Delta)) + xlab("Number of Items") +
 scale_fill_manual(values=greyPalette)


