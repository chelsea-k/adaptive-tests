# Script for plotting Figure 9 from Krantsevich et al., 
# "BAYESIAN DECISION THEORY FOR TREE-BASED ADAPTIVE SCREENING TESTS 
#     WITH AN APPLICATION TO YOUTH DELINQUENCY"

source('adaptive_tests/code/util_functions.R')
library(ggplot2)
theme_set(theme_bw(base_size=14))

# Set parameters
maxIPP_vals <- c(2:15)
n_maxIPP <- length(maxIPP_vals)
w <- 0.5

# Set up folders and load data
data_dir <- "output/out_of_sample/all/synthetic_data"
results_dir <- "output/out_of_sample/all/results"
plots_dir <- "output/plots"
dir.create(plots_dir, recursive = TRUE)

data_train <- read.csv("preprocessed_original_data/IMC_data_train_preprocessed.csv")
data_test <- read.csv("preprocessed_original_data/IMC_data_test_preprocessed.csv")
  
synth_df <- read.csv(file.path(data_dir, "synth_uncertainty_XB.csv"))
n_mcmc <- length(unique(synth_df$post.idx))
n_rows_total <- n_maxIPP*n_mcmc

p_maxIPP_XB_synth <- read.csv(file.path(results_dir, "p_maxIPP_XB.synth_uncertainty_XB.csv"))
p_maxIPP_XB_test <- read.csv(file.path(results_dir, "p_maxIPP_XB.test.csv"))

load(file.path(data_dir, "XB_predict"))
p_XBART_test <- XB_predict$test


# Compute optimal cutoffs and predicted classes 
XBART_cutoff <- get_cutoff(as.factor(synth_df$y.draw), synth_df$phat.mean, w)
maxIPP_cutoff <- rep(NA, n_maxIPP)
for (i in seq_along(maxIPP_vals)){
  maxIPP <- maxIPP_vals[[i]]
  col <- which(colnames(p_maxIPP_XB_synth)==paste0("m.", maxIPP))
  maxIPP_cutoff[[i]] <- get_cutoff(as.factor(synth_df$y.draw),p_maxIPP_XB_synth[,col],w)
}

# Compute delta results for training data, for all values of maxIPP
delta_draws_train <- data.frame(matrix(NA, nrow=n_rows_total, ncol=3))
colnames(delta_draws_train) <- c("maxIPP", "post.idx", "Delta")
for (j in 1:n_mcmc){
  if (j%%100==0) {cat(paste0("On iteration j = ", j, " out of ", n_mcmc, "\n"))}
  post_idx_rows <- which(synth_df$post.idx==j)
  temp_df <- synth_df[post_idx_rows,]
  y_temp <- temp_df$y.draw
  XB_metrics <- get_utility(y_temp,prob=temp_df$phat.mean,cut=XBART_cutoff,type="prob",w=w)
  for (i in seq_along(maxIPP_vals)){
    temp_df2 <- p_maxIPP_XB_synth[post_idx_rows,]
    maxIPP <- maxIPP_vals[i]
    col <- which(colnames(p_maxIPP_XB_synth)==paste0("m.", maxIPP))
    CT_metrics <- get_utility(y_temp,prob=temp_df2[,col],cut=maxIPP_cutoff[[i]],type="prob",w=w)
    delta_draws_train[n_maxIPP*(j-1)+i,] <- data.frame(maxIPP=maxIPP, 
                                                      post.idx=j, 
                                                      Delta=CT_metrics$util-XB_metrics$util)
  }
}

# Compute delta results for testing data, for all values of maxIPP
delta_draws_test <- data.frame(matrix(NA, nrow=n_maxIPP, ncol=3))
colnames(delta_draws_test) <- c("maxIPP", "post.idx", "Delta")
XB_met_test <- get_utility(1*data_test$y,prob=p_XBART_test,cut=XBART_cutoff,type="prob",w=w)
for (i in seq_along(maxIPP_vals)){
  maxIPP <- maxIPP_vals[[i]]
  col <- which(colnames(p_maxIPP_XB_test)==paste0("m.", maxIPP))
  CT_met_test <- get_utility(data_test$y,prob=p_maxIPP_XB_test[,col],
                           cut=maxIPP_cutoff[[i]],type="prob",w=w)
  delta_draws_test[i,] <- data.frame(maxIPP=maxIPP, 
                                           post.idx=0, 
                                           Delta=CT_met_test$util-XB_met_test$util)
}

# Post-process results and create plot for Figure 9
delta_draws_all <- rbind(cbind(delta_draws_train, Type = "Predicted"), 
                        cbind(delta_draws_test, Type = "Actual"))

delta_draws_all$maxIPP <- as.factor(delta_draws_all$maxIPP)
delta_draws_all$Type <- factor(delta_draws_all$Type, levels=c("Predicted", "Actual"))
cbPalette <- c("#000000","#D55E00", "#CC79A7")
ggplot(delta_draws_all, aes(x=maxIPP, y=Delta, color=Type)) +
  geom_boxplot(outlier.shape=NA, position="identity") + 
  theme(plot.title = element_text(hjust = 0.5))  + 
  ggtitle("Out of Sample Utility Differences") + ylab(expression(Delta)) +
  scale_color_manual(values=cbPalette) +
  xlab("Number of Items")  

ggsave(file.path(plots_dir,"Fig9.png"), height = 3.25, width = 9, units = "in", dpi = 200)



