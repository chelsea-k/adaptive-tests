# Script for generating synthetic item response data using the  Gaussian 
# copula factor and Accelerated Bayesian Additive Regression Trees models

library(bfa.mod)
library(XBART)
library(plyr)

out_of_sample <- FALSE
subpopulation <- FALSE

# paths to item response data for fitting models
IR_data_dir <- "simulated_data"
output_dir <- "output"
if (out_of_sample) {
  IR_data_train_path <- file.path(IR_data_dir, "item_response_data_train.csv")
  IR_data_test_path <- file.path(IR_data_dir, "item_response_data_test.csv")
} else {
  IR_data_path <- file.path(IR_data_dir, "item_response_data_all.csv")
}

# directories for model and data storage
folder <- ifelse(out_of_sample, 
                 file.path(output_dir, "out_of_sample"), 
                 file.path(output_dir, "in_sample"))
folder <- file.path(folder, ifelse(subpopulation, "subpopulation", "all"))
data_dir <- file.path(folder, "synthetic_data")
model_dir <- file.path(folder, "model_fits")
dir.create(data_dir, recursive = TRUE)
dir.create(model_dir, recursive = TRUE)

########################## Hyperparamters ##############################

# Sampling parameters
n_mcmc <- 1000   # number of "sample populations" to draw, D 
n_samp <- 1000  # number of data points in each sample population, N
n_prune_samp <- 100 

# BFA parameters for fitting f(X)
num_factor <- 3
nburn <- 5000
if (subpopulation){
  cond_vars <- list(Age = 15)   # conditioning variables for subpopulation
  cond_type <- list(">=")
} else {
  cond_vars <- NA
  cond_type <- NA
}

# XBART.multinomial parameters for fitting f(Y|X)
num_trees <- 30
burnin <- 100 
Nmin <- 4
max_depth <- 250
num_cutpoints <- 7
weight <- 1

########################## Data Preparation ##############################

# read in data
if(out_of_sample){
  data_train <- read.csv(IR_data_train_path)
  data_test <- read.csv(IR_data_test_path)
} else {
  data_train <- read.csv(IR_data_path)
  data_test <- data_train
}
n_train <- nrow(data_train)
n_test <- nrow(data_test)
n_cond <- ifelse(is.na(cond_vars), 0, length(cond_vars))
p <- ncol(data_train) - 1 - n_cond # num columns in synthetic X matrix
item_demo_cols <- which(colnames(data_train) != "y") 
item_cols <- which(!(colnames(data_train) %in% c("y", names(cond_vars)))) 

########################## Model Fitting ##############################

# check for Gaussian copula factor model & fit if not present
BFA_model_file <- file.path(model_dir, "fit_BFA")
if(file.exists(BFA_model_file)){
  cat("Loading BFA model...\n")
  load(BFA_model_file)
} else{
  cat("Fitting BFA model...\n")		
  fit_BFA <- bfa_copula(~., data=data_train[,item_demo_cols], 
                       num.factor = num_factor,
                       factor.scales = FALSE, 
                       keep.scores = FALSE, 
                       nburn = nburn, 
                       nsim = 2*n_mcmc, 
                       loading.prior = "gdp", 
                       imh = FALSE)
  save(fit_BFA, file = BFA_model_file, ascii=TRUE)
}

##### fit XBART model ##### 
# first compute needed params
p_categorical <- length(item_cols)
mtry <- p_categorical + 1
XB_num_sweeps <- 2*n_mcmc + burnin
XB_postburn_idx <- (burnin + 1):XB_num_sweeps

# fit model; XBART crashed with all categorical inputs -> added dummy rnorm column
cat("Fitting XBART model...\n")
fit_XBART <- XBART.multinomial(y = as.matrix(data_train$y), 
                               num_class = 2, 
                               X = as.matrix(cbind(rnorm(n_train), data_train[,item_cols])), 
                               Xtest = as.matrix(cbind(rnorm(n_test), data_test[,item_cols])),
                               num_trees = num_trees, 
                               num_sweeps = XB_num_sweeps, 
                               max_depth = max_depth, 
                               Nmin = Nmin, 
                               num_cutpoints = num_cutpoints, 
                               alpha = 0.95, 
                               beta = 1.25, 
                               tau_a = 1, 
                               tau_b = 1,
                               no_split_penality = 1, 
                               burnin = burnin, 
                               mtry = mtry, 
                               p_categorical = p_categorical, 
                               verbose = FALSE, 
                               parallel = FALSE, 
                               random_seed = NULL, 
                               sample_weights_flag = TRUE, 
                               separate_tree = FALSE,
                               weight = weight)
save(fit_XBART, file = file.path(model_dir, "fit_XBART"), ascii=TRUE)

# Predict on train and test data and save output
XB_pred_train <- colMeans(fit_XBART$yhats_train[XB_postburn_idx,,2])
XB_pred_test <- colMeans(fit_XBART$yhats_test[XB_postburn_idx,,2])

XB_predict <- list(train = XB_pred_train, test = XB_pred_test)
save(XB_predict, file = file.path(data_dir, "XB_predict"), ascii = TRUE)

########################## Sampling ##############################

# Draw posterior samples from the fitted Gaussian copula factor model;
# match up posterior indices to draw probabilities and risk class from the fitted XBART model;
# Also compute posterior mean probability \bar{E}(Y|x) for fitting regression tree 
synth_data_treefitting <- array(NA, dim=c(n_mcmc, n_samp, p+4))
synth_data_uncertainty <- array(NA, dim=c(n_mcmc, n_samp, p+4))
prune_data_treefitting <- array(NA, dim=c(n_mcmc, n_prune_samp, p+4))
prune_data_uncertainty <- array(NA, dim=c(n_mcmc, n_prune_samp, p+4))

for (j in 1:(2*n_mcmc)) {
  cat(sprintf("============ Predicting for iteration %d out of %d ============\n", j, 2*n_mcmc))		
  # Draw samples
  Xtilde <- predict_idx(fit_BFA, post.idx = j, n.samp = n_samp + n_prune_samp, 
                       cond.vars = cond_vars, cond.type = cond_type)
  X_item_cols <- which(colnames(Xtilde) != "Age")
  p_XBART_draw <- predict.XBARTmultinomial(fit_XBART, 
                                           X=as.matrix(cbind(rnorm(nrow(Xtilde)), 
                                                            Xtilde[,X_item_cols])), 
                                           iteration = as.integer(j+burnin))
  p_XBART_draw <- p_XBART_draw$yhats[,,2]
  Ytilde_draw <- rbinom(n = length(p_XBART_draw), size=1, prob=p_XBART_draw)
  
  p_XBART_mean_pred <- predict.XBARTmultinomial(fit_XBART, 
                                           X=as.matrix(cbind(rnorm(nrow(Xtilde)), 
                                                             Xtilde[,X_item_cols])))
  p_XBART_mean <- colMeans(p_XBART_mean_pred$yhats[XB_postburn_idx,,2])
  Ytilde_mean <- rbinom(n = length(p_XBART_mean), size=1, prob=p_XBART_mean)
  
  # Store data
  temp <- as.matrix(cbind(Xtilde, p_XBART_mean, Ytilde_mean, p_XBART_draw, Ytilde_draw))
  if(j%%2 == 0) {
    synth_data_treefitting[j/2,,] <- temp[1:n_samp,]
    prune_data_treefitting[j/2,,] <- temp[(n_samp+1):(n_samp + n_prune_samp),]
  } else {
    synth_data_uncertainty[(j+1)/2,,] <- temp[1:n_samp,]
    prune_data_uncertainty[(j+1)/2,,] <- temp[(n_samp+1):(n_samp + n_prune_samp),]
  }
}

# Data postprocessing and saving 
synth_treefitting_XB <- adply(synth_data_treefitting, .margins = 1, .id="post.idx")
synth_uncertainty_XB <- adply(synth_data_uncertainty, .margins = 1, .id="post.idx")
prune_treefitting_XB <- adply(prune_data_treefitting, .margins = 1, .id="post.idx")
prune_uncertainty_XB <- adply(prune_data_uncertainty, .margins = 1, .id="post.idx")
Xtilde_cols <- which(!(colnames(data_train) %in% c("y", names(cond_vars))))
Xtilde_colnames <- colnames(data_train)[Xtilde_cols]
colnames(synth_treefitting_XB) <- c('post.idx', Xtilde_colnames, 'phat.mean', 'y.mean', 'phat.draw', 'y.draw')
colnames(synth_uncertainty_XB) <- c('post.idx', Xtilde_colnames, 'phat.mean', 'y.mean', 'phat.draw', 'y.draw')
colnames(prune_treefitting_XB) <- c('post.idx', Xtilde_colnames, 'phat.mean', 'y.mean', 'phat.draw', 'y.draw')
colnames(prune_uncertainty_XB) <- c('post.idx', Xtilde_colnames, 'phat.mean', 'y.mean', 'phat.draw', 'y.draw')


write.csv(synth_treefitting_XB, file.path(data_dir, "synth_treefitting_XB.csv"), row.names = FALSE)
write.csv(synth_uncertainty_XB, file.path(data_dir, "synth_uncertainty_XB.csv"), row.names = FALSE)
write.csv(prune_treefitting_XB, file.path(data_dir, "prune_treefitting_XB.csv"), row.names = FALSE)
write.csv(prune_uncertainty_XB, file.path(data_dir, "prune_uncertainty_XB.csv"), row.names = FALSE)
