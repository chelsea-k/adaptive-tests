# Script for generating synthetic item response data using the  Gaussian 
# copula factor and Accelerated Bayesian Additive Regression Trees models

library(bfa.mod)
library(XBART)
library(plyr)

out.of.sample = TRUE
subpopulation = FALSE

########################## Hyperparamters ##############################

# Sampling parameters
n.mcmc = 100   # number of "sample populations" to draw, D
n.samp = 100   # number of data points per sample, N
n.prune.samp = 10   # n.prune.samp * n.mcmc is number of o.o.s data 
                    #    for pruning tree under maxIPP strategy

# BFA parameters for fitting f(X)
num.factor = 3
nburn = 500
if (subpopulation){
  cond.vars = list(Age=15)   # conditioning variables for subpopulation
  cond.type = list(">=")
} else {
  cond.vars=NA
  cond.type=NA
}

# XBART.multinomial parameters for fitting f(Y|X)
num_trees = 30
burnin = 100
Nmin = 4
max_depth = 250
num_cutpoints = 7
weight = 1

########################## Data Preparation ##############################

# read in data
if(out.of.sample){
  data.train = read.csv("data/item_response_data/item.response.data.train.csv")
  data.test = read.csv("data/item_response_data/item.response.data.test.csv")
} else {
  data.train = read.csv("data/item_response_data/item.response.data.all.csv")
  data.test = data.train
}
n.train = nrow(data.train)
n.test = nrow(data.test)
n.total = n.train + n.test
n.cond = ifelse(is.na(cond.vars), 0, length(cond.vars))
p = ncol(data.train) - 1 - n.cond # num columns in synthetic X matrix
item.demo.cols = which(colnames(data.train) != "y") 
item.cols = which(!(colnames(data.train) %in% c("y", "Age"))) 

# directories for model and data storage
folder = ifelse(out.of.sample, "data/out_of_sample/", "data/in_sample/")
folder = paste0(folder, ifelse(subpopulation, "subpopulation/", "all/"))
data.dir = paste0(folder, "synthetic_data/")
model.dir = paste0(folder, "model_fits/")
dir.create(data.dir, recursive = TRUE)
dir.create(model.dir, recursive = TRUE)

########################## Model Fitting ##############################

# fit Gaussian copula factor model
fit.BFA = bfa_copula(~., data=data.train[,item.demo.cols], 
                     num.factor = num.factor,
                     factor.scales = FALSE, 
                     keep.scores = FALSE, 
                     nburn = nburn, 
                     nsim = 2*n.mcmc, 
                     loading.prior = "gdp", 
                     imh = FALSE)
save(fit.BFA, file=paste0(model.dir, "fit.BFA"), ascii=TRUE)


# fit XBART model; XBART was giving strange results when all data were categorical, 
# added a dummy rnorm column as a temporary fix

# first compute extra params
p_categorical = length(item.cols)
mtry = p_categorical+1
fit.XBART = XBART.multinomial(y = as.matrix(data.train$y), num_class=2, 
                              X = as.matrix(cbind(rnorm(n.train), data.train[,item.cols])), 
                              Xtest = as.matrix(cbind(rnorm(n.train), data.test[,item.cols])),
                              num_trees = num_trees, 
                              num_sweeps = 2*n.mcmc + burnin, 
                              max_depth = max_depth, 
                              Nmin = Nmin, 
                              num_cutpoints = num_cutpoints, 
                              alpha = 0.95, 
                              beta = 1.25, 
                              tau = 1, 
                              no_split_penality = 1, 
                              weight = weight, 
                              burnin = burnin, 
                              mtry = mtry, 
                              p_categorical = p_categorical, 
                              kap = 1, 
                              s = 1, 
                              verbose = FALSE, 
                              parallel = FALSE, 
                              random_seed = NULL, 
                              sample_weights_flag = TRUE)
save(fit.XBART, file=paste0(model.dir, "fit.XBART"), ascii=TRUE)

# Predict on train and test data and save output
XB.predict = list(test = colMeans(fit.XBART$yhats_test[,,2]))
save(XB.predict, file=paste0(data.dir, "XB.predict"), ascii = TRUE)

########################## Sampling ##############################

# Draw samples from every other posterior draw of the fitted Gaussian copula factor model;
# match up posterior indices to draw probabilities and risk class from the fitted XBART model;
# Also compute posterior mean probability for fitting regression tree to
synth.data = array(NA, dim=c(n.mcmc, n.samp, p+3))
prune.data = array(NA, dim=c(n.mcmc, n.prune.samp, p+3))
for (j in 1:n.mcmc) {
  cat(sprintf("----------------------Predicting for iteration %d out of %d \n", j, n.mcmc))		
  # Draw samples
  Xtilde = predict_idx(fit.BFA, post.idx = 2L*j, n.samp = n.samp + n.prune.samp, 
                       cond.vars=cond.vars, cond.type =cond.type)
  X.item.cols = which(colnames(Xtilde) != "Age")
  p.XBART.draw = predict.XBARTmultinomial(fit.XBART, 
                                          X=as.matrix(cbind(rnorm(nrow(Xtilde)), 
                                                            Xtilde[,X.item.cols])), 
                                          iteration = 2L*j)
  p.XBART.draw = p.XBART.draw$yhats[,,2]
  Ytilde = rbinom(n=length(p.XBART.draw), size=1, prob=p.XBART.draw)
  p.XBART.mean = 0
  for (k in 1:n.mcmc){
    p.XBART = predict.XBARTmultinomial(fit.XBART,
                                       X=as.matrix(cbind(rnorm(nrow(Xtilde)), 
                                                         Xtilde[,X.item.cols])), 
                                       iteration=2L*k)
    p.XBART.mean =  p.XBART.mean + p.XBART$yhats[,,2]
  }
  p.XBART.mean = p.XBART.mean/n.mcmc
  
  # Store data
  temp = as.matrix(cbind(Xtilde, p.XBART.mean, p.XBART.draw, Ytilde))
  synth.data[j,,] = temp[1:n.samp,]
  prune.data[j,,] = temp[(n.samp+1):(n.samp + n.prune.samp),]
}

# Data postprocessing and saving 
synth.XB = adply(synth.data, .margins = 1, .id="post.idx")
prune.XB = adply(synth.data, .margins = 1, .id="post.idx")
Xtilde.cols = which(!(colnames(data.train) %in% c("y", names(cond.vars))))
Xtilde.colnames = colnames(data.train)[Xtilde.cols]
colnames(synth.XB) = colnames(prune.XB) = c('post.idx', Xtilde.colnames, 
                       'phat', 'phat.draw', 'y')
write.csv(synth.XB, paste0(data.dir, "synth.XB.csv"), row.names = FALSE)
write.csv(prune.XB, paste0(data.dir, "prune.XB.csv"), row.names = FALSE)
