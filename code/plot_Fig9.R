# Script for plotting Figure 9 from Krantsevich et al., 
# "BAYESIAN DECISION THEORY FOR TREE-BASED ADAPTIVE SCREENING TESTS 
#     WITH AN APPLICATION TO YOUTH DELINQUENCY"

source('code/util_functions.R')
library(ggplot2)

# Set parameters
maxIPP.vals = c(4,10)
n.maxIPP = length(maxIPP.vals)
w = 0.5

# Set up data folders and load data
data.dir = "data/in_sample/all/synthetic_data/"
results.dir = "data/in_sample/all/results/"

synth.df = read.csv(paste0(data.dir, "synth.XB.csv"))
n.mcmc = length(unique(synth.df$post.idx))
n.methods = 7
n.rows.total = n.methods*n.maxIPP*n.mcmc

p.maxIPP.XB = read.csv(paste0(results.dir, "p.maxIPP.XB.synth.csv"))
p.maxIPP.RF = read.csv(paste0(results.dir, "p.maxIPP.RF.synth.csv"))
p.maxdepth.XB = read.csv(paste0(results.dir, "p.maxDepth.XB.synth.csv"))
p.maxdepth.RF = read.csv(paste0(results.dir, "p.maxDepth.RF.synth.csv"))
y.class.real = read.csv(paste0(results.dir, "y.real.synth.csv"))
y.class.RF = read.csv(paste0(results.dir, "y.RF.synth.csv"))
y.class.XB = read.csv(paste0(results.dir, "y.XB.synth.csv"))

# Compute optimal cutoffs and predicted classes 
XB.cutoff = get_cutoff(as.factor(synth.df$y), synth.df$phat, w)
maxIPP.XB.cut = maxIPP.RF.cut = maxdepth.XB.cut = maxdepth.RF.cut = rep(NA, n.maxIPP)
m.colnames = paste0("m.", maxIPP.vals)
m.cols = which(colnames(p.maxIPP.XB) %in% m.colnames)
for (i in seq_along(maxIPP.vals)){
  maxIPP.XB.cut[[i]] = get_cutoff(as.factor(synth.df$y), p.maxIPP.XB[,m.cols[i]], w)
  maxIPP.RF.cut[[i]] = get_cutoff(as.factor(synth.df$y), p.maxIPP.RF[,m.cols[i]], w)
  maxdepth.XB.cut[[i]] = get_cutoff(as.factor(synth.df$y), p.maxdepth.XB[,m.cols[i]], w)
  maxdepth.RF.cut[[i]] = get_cutoff(as.factor(synth.df$y), p.maxdepth.RF[,m.cols[i]], w)
}
y.XBART = 1*(synth.df$phat >= XB.cutoff)
y.maxIPP.XB = get_class(p.maxIPP.XB, maxIPP.XB.cut, maxIPP.vals, m.cols)
y.maxIPP.RF = get_class(p.maxIPP.RF, maxIPP.RF.cut, maxIPP.vals, m.cols)
y.maxDepth.XB = get_class(p.maxdepth.XB, maxdepth.XB.cut, maxIPP.vals, m.cols)
y.maxDepth.RF = get_class(p.maxdepth.RF, maxdepth.RF.cut, maxIPP.vals, m.cols)
# Extract columns of maxIPP vals being used
y.class.real = y.class.real[,which(colnames(y.class.real) %in% m.colnames)]
y.class.RF = y.class.RF[,which(colnames(y.class.RF) %in% m.colnames)]
y.class.XB = y.class.XB[,which(colnames(y.class.XB) %in% m.colnames)]
  
# Compute delta results for all values of maxIPP for all methods
delta.draws.all = data.frame(matrix(NA, nrow=n.rows.total, ncol=4))
colnames(delta.draws.all) = c("num.items", "post.idx", "Delta", "Method")

for (j in 1:n.mcmc){
  if (j%%100==0) {cat(paste0("On iteration j = ", j, " out of ", n.mcmc, "\n"))}
  post.rows = which(synth.df$post.idx==j)
  temp.df = synth.df[post.rows,]
  y.temp = temp.df$y
  XB.metrics = get_utility(y.temp,y.pred=y.XBART[post.rows],type="class",w=w)
  for (i in seq_along(maxIPP.vals)){
    maxIPP = maxIPP.vals[[i]]
    maxIPP.XB.metrics = get_utility(y.temp,y.pred=y.maxIPP.XB[post.rows,i],type="class",w=w)
    maxdepth.XB.metrics = get_utility(y.temp,y.pred=y.maxDepth.XB[post.rows,i],type="class",w=w)
    maxIPP.RF.metrics = get_utility(y.temp,y.pred=y.maxIPP.RF[post.rows,i],type="class",w=w)
    maxdepth.RF.metrics = get_utility(y.temp,y.pred=y.maxDepth.RF[post.rows,i],type="class",w=w)
    class.XB.metrics  = get_utility(y.temp,y.pred=y.class.XB[post.rows,i],type="class",w=w)
    class.RF.metrics  = get_utility(y.temp,y.pred=y.class.RF[post.rows,i],type="class",w=w)
    class.real.metrics  = get_utility(y.temp,y.pred=y.class.real[post.rows,i],type="class",w=w)
    
    temp.results = data.frame(num.items = maxIPP.vals[[i]], post.idx = j, 
                              Delta = c(maxIPP.XB.metrics$util, maxdepth.XB.metrics$util, 
                                        maxIPP.RF.metrics$util, maxdepth.RF.metrics$util,
                                        class.XB.metrics$util, class.RF.metrics$util, 
                                        class.real.metrics$util) - XB.metrics$util,
                              Method = c("Regression (maxIPP, GCFM + XBART)", 
                                         "Regression (max depth, GCFM + XBART)", 
                                         "Regression (maxIPP, Perturb + RF)", 
                                         "Regression (max depth, Perturb + RF)",
                                         "Classification (GCFM + XBART)", 
                                         "Classification (Perturb + RF)", 
                                         "Classification (Real Data)"))
    
    delta.draws.all[c(1:n.methods)+(n.methods*(i-1))+(n.methods*n.maxIPP*(j-1)),]=temp.results
    }
}
write.csv(delta.draws.all,paste0(results.dir,"delta.draws.all.methods.",w,".csv"),row.names=F)

# Post-process results and create plot for Figure 9
cbPalette = c("#000000","#009E73", "#0072B2", "#D55E00","#F0E442",  "#56B4E9","#E69F00")
delta.draws.all$num.items = as.factor(delta.draws.all$num.items)
delta.draws.all$Method = factor(delta.draws.all$Method, 
                                 levels=c("Regression (maxIPP, GCFM + XBART)", 
                                          "Regression (max depth, GCFM + XBART)", 
                                          "Regression (maxIPP, Perturb + RF)", 
                                          "Regression (max depth, Perturb + RF)",
                                          "Classification (GCFM + XBART)", 
                                          "Classification (Perturb + RF)", 
                                          "Classification (Real Data)"))

ggplot(delta.draws.all, aes(y=Delta, fill=Method)) +
  geom_boxplot(outlier.shape=NA, position=position_dodge2(padding=0.1, width=10)) + 
  theme(plot.title = element_text(hjust = 0.5))  + 
  ggtitle(paste0("Utility Differences (Populating the Action Space)")) + 
  ylab(expression(Delta)) + xlab("Number of Items") +
  scale_fill_manual(values=cbPalette) +
  facet_grid(cols = vars(num.items), switch="both") +
  theme(axis.text.x=element_blank(), axis.ticks.x = element_blank())
