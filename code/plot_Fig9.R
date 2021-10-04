# Script for plotting Figure 9 from Krantsevich et al., 
# "BAYESIAN DECISION THEORY FOR TREE-BASED ADAPTIVE SCREENING TESTS 
#     WITH AN APPLICATION TO YOUTH DELINQUENCY"

source('adaptive_tests/code/util_functions.R')
library(ggplot2)
library(dplyr)
theme_set(theme_bw(base_size=14))

# Is this for Figure 9 or the appendix version with utility-based cutoffs?
appendix_version <- FALSE

# Should we use the fitting data to get the optimal cutoff?
use_fitdata_for_cutoff <- TRUE

# Set parameters
w <- 0.5
if(appendix_version == TRUE){
  maxIPP_vals <- c(2:15)
  n_methods <- 2
} else {
  maxIPP_vals <- c(3, 9, 15)
  n_methods <- 7
}
n_maxIPP <- length(maxIPP_vals)


# Set up data folders and load data
original_data_dir <- "preprocessed_original_data"
data_dir <- "output/in_sample/all/synthetic_data"
results_dir <- "output/in_sample/all/results"
plots_dir <- "output/plots"

dir.create(results_dir, recursive = TRUE)
dir.create(plots_dir, recursive = TRUE)

synth_uncertainty_XB <- read.csv(file.path(data_dir, "synth_uncertainty_XB.csv"))

n_mcmc <- length(unique(synth_uncertainty_XB$post.idx))
n_rows_total <- n_methods*n_maxIPP*n_mcmc

p_maxIPP_XB_uncertainty <- read.csv(file.path(results_dir, "p_maxIPP_XB.synth_uncertainty_XB.csv"))
p_maxIPP_RF_uncertainty <- read.csv(file.path(results_dir, "p_maxIPP_RF.synth_uncertainty_XB.csv"))
p_maxDepth_XB_uncertainty <- read.csv(file.path(results_dir, "p_maxDepth_XB.synth_uncertainty_XB.csv"))
p_maxDepth_RF_uncertainty <- read.csv(file.path(results_dir, "p_maxDepth_RF.synth_uncertainty_XB.csv"))
y_class_real_uncertainty <- read.csv(file.path(results_dir, "y_real.synth_uncertainty_XB.csv"))
y_class_RF_uncertainty <- read.csv(file.path(results_dir, "y_RF.synth_uncertainty_XB.csv"))
y_class_XB_uncertainty <- read.csv(file.path(results_dir, "y_XB.synth_uncertainty_XB.csv"))
y_class_util_uncertainty <- read.csv(file.path(results_dir, "y_util.synth_uncertainty_XB.csv"))

if(use_fitdata_for_cutoff == FALSE) {
  # assign fitting data for computing cutoffs
  synth_XB_cutoff <- synth_uncertainty_XB
  synth_RF_cutoff <- synth_uncertainty_XB
  
  # assign predicted probabilities for computing cutoffs
  p_maxIPP_XB_cutoff <- p_maxIPP_XB_uncertainty
  p_maxIPP_RF_cutoff <- p_maxIPP_RF_uncertainty
  p_maxDepth_XB_cutoff <- p_maxDepth_XB_uncertainty
  p_maxDepth_RF_cutoff <- p_maxDepth_RF_uncertainty
 } else {
   # read in fitting data for computing cutoffs
   synth_XB_cutoff <- read.csv(file.path(data_dir, "synth_treefitting_XB.csv"))
   synth_RF_cutoff <- read.csv(file.path(data_dir, "synth_treefitting_RF.csv"))
   synth_RF_cutoff <- rename(synth_RF_cutoff, y.draw = y, phat.mean = phat)
   
   # read in predicted probabilities for computing cutoffs
   p_maxIPP_XB_cutoff <- read.csv(file.path(results_dir, "p_maxIPP_XB.synth_treefitting_XB.csv"))
   p_maxIPP_RF_cutoff <- read.csv(file.path(results_dir, "p_maxIPP_RF.synth_treefitting_RF.csv"))
   p_maxDepth_XB_cutoff <- read.csv(file.path(results_dir, "p_maxDepth_XB.synth_treefitting_XB.csv"))
   p_maxDepth_RF_cutoff <- read.csv(file.path(results_dir, "p_maxDepth_RF.synth_treefitting_RF.csv"))
}

# Compute optimal cutoffs and predicted classes 
XB_cutoff <- get_cutoff(as.factor(synth_XB_cutoff$y.draw), synth_XB_cutoff$phat.mean, w)
maxIPP_XB_cut <- rep(NA, n_maxIPP)
maxIPP_RF_cut <- rep(NA, n_maxIPP)
maxdepth_XB_cut <- rep(NA, n_maxIPP)
maxdepth_RF_cut <- rep(NA, n_maxIPP)
m_colnames <- paste0("m.", maxIPP_vals)
m_cols <- which(colnames(p_maxIPP_XB_uncertainty) %in% m_colnames)
for (i in seq_along(maxIPP_vals)){
  maxIPP_XB_cut[[i]] <- get_cutoff(as.factor(synth_XB_cutoff$y.draw), p_maxIPP_XB_cutoff[,m_cols[i]], w)
  maxIPP_RF_cut[[i]] <- get_cutoff(as.factor(synth_RF_cutoff$y.draw), p_maxIPP_RF_cutoff[,m_cols[i]], w)
  maxdepth_XB_cut[[i]] <- get_cutoff(as.factor(synth_XB_cutoff$y.draw), p_maxDepth_XB_cutoff[,m_cols[i]], w)
  maxdepth_RF_cut[[i]] <- get_cutoff(as.factor(synth_RF_cutoff$y.draw), p_maxDepth_RF_cutoff[,m_cols[i]], w)
}
y_XBART <- 1*(synth_uncertainty_XB$phat.mean >= XB_cutoff)
y_maxIPP_XB_uncertainty <- get_class(p_maxIPP_XB_uncertainty, maxIPP_XB_cut, maxIPP_vals, m_cols)
y_maxIPP_RF_uncertainty <- get_class(p_maxIPP_RF_uncertainty, maxIPP_RF_cut, maxIPP_vals, m_cols)
y_maxDepth_XB_uncertainty <- get_class(p_maxDepth_XB_uncertainty, maxdepth_XB_cut, maxIPP_vals, m_cols)
y_maxDepth_RF_uncertainty <- get_class(p_maxDepth_RF_uncertainty, maxdepth_RF_cut, maxIPP_vals, m_cols)

# Extract columns of maxIPP vals being used
y_class_real_uncertainty <- y_class_real_uncertainty[,which(colnames(y_class_real_uncertainty) %in% m_colnames)]
y_class_RF_uncertainty <- y_class_RF_uncertainty[,which(colnames(y_class_RF_uncertainty) %in% m_colnames)]
y_class_XB_uncertainty <- y_class_XB_uncertainty[,which(colnames(y_class_XB_uncertainty) %in% m_colnames)]
y_class_util_uncertainty <- y_class_util_uncertainty[,which(colnames(y_class_util_uncertainty) %in% m_colnames)]

# Compute delta results for all values of maxIPP for all methods
delta_draws_all <- data.frame(matrix(NA, nrow=n_rows_total, ncol=4))
colnames(delta_draws_all) <- c("num.items", "post.idx", "Delta", "Method")

for (j in 1:n_mcmc){
  if (j%%100==0) {cat(paste0("On iteration j = ", j, " out of ", n_mcmc, "\n"))}
  post_rows <- which(synth_uncertainty_XB$post.idx==j)
  temp_df <- synth_uncertainty_XB[post_rows,]
  y_temp <- temp_df$y.draw
  XB_metrics <- get_utility(y_temp,y.pred=y_XBART[post_rows],type="class",w=w)
  for (i in seq_along(maxIPP_vals)){
    maxIPP <- maxIPP_vals[[i]]
    maxIPP_XB_metrics <- get_utility(y_temp,y.pred=y_maxIPP_XB_uncertainty[post_rows,i],type="class",w=w)
    maxdepth_XB_metrics <- get_utility(y_temp,y.pred=y_maxDepth_XB_uncertainty[post_rows,i],type="class",w=w)
    maxIPP_RF_metrics <- get_utility(y_temp,y.pred=y_maxIPP_RF_uncertainty[post_rows,i],type="class",w=w)
    maxdepth_RF_metrics <- get_utility(y_temp,y.pred=y_maxDepth_RF_uncertainty[post_rows,i],type="class",w=w)
    class_XB_metrics  <- get_utility(y_temp,y.pred=y_class_XB_uncertainty[post_rows,i],type="class",w=w)
    class_RF_metrics  <- get_utility(y_temp,y.pred=y_class_RF_uncertainty[post_rows,i],type="class",w=w)
    class_real_metrics  <- get_utility(y_temp,y.pred=y_class_real_uncertainty[post_rows,i],type="class",w=w)
    class_util_metrics  <- get_utility(y_temp,y.pred=y_class_util_uncertainty[post_rows,i],type="class",w=w)
    
    store_rows <- c(1:n_methods)+(n_methods*(i-1))+(n_methods*n_maxIPP*(j-1))
    if(appendix_version == TRUE){
      temp_results <- data.frame(num.items = maxIPP_vals[[i]], post.idx = j, 
                                 Delta = c(maxIPP_XB_metrics$util - XB_metrics$util,  
                                           class_util_metrics$util- XB_metrics$util), 
                                 Method = c("Regression (maxIPP, GCFM + XBART)",  
                                            "Classification (Utility Outcomes)"))
    } else {
        temp_results <- data.frame(num.items = maxIPP_vals[[i]], post.idx = j, 
                                  Delta = c(maxIPP_XB_metrics$util, maxdepth_XB_metrics$util, 
                                            maxIPP_RF_metrics$util, maxdepth_RF_metrics$util,
                                            class_XB_metrics$util, class_RF_metrics$util, 
                                            class_real_metrics$util) - XB_metrics$util,
                                  Method = c("Regression (maxIPP, GCFM + XBART)", 
                                             "Regression (max depth, GCFM + XBART)", 
                                             "Regression (maxIPP, Perturb + RF)", 
                                             "Regression (max depth, Perturb + RF)",
                                             "Classification (GCFM + XBART)", 
                                             "Classification (Perturb + RF)", 
                                             "Classification (Real Data)"))
    }
    delta_draws_all[store_rows,] <- temp_results
  }
}

# Post-process results and create plot for Figure 9
cbPalette <- c("#000000","#009E73", "#0072B2", "#D55E00","#F0E442",  "#56B4E9","#E69F00")
delta_draws_all$num.items <- as.factor(delta_draws_all$num.items)
if(appendix_version == TRUE) {
  levels = c("Regression (maxIPP, GCFM + XBART)", 
             "Classification (Utility Outcomes)")
} else {
  levels = c("Regression (maxIPP, GCFM + XBART)", 
             "Regression (max depth, GCFM + XBART)", 
             "Regression (maxIPP, Perturb + RF)", 
             "Regression (max depth, Perturb + RF)",
             "Classification (GCFM + XBART)", 
             "Classification (Perturb + RF)", 
             "Classification (Real Data)")
}
delta_draws_all$Method <- factor(delta_draws_all$Method, 
                                 levels= levels)

ggplot(delta_draws_all, aes(y=Delta, fill=Method)) +
  geom_boxplot(outlier.shape=NA, position=position_dodge2(padding=0.1, width=10)) + 
  theme(plot.title = element_text(hjust = 0.5))  + 
  ggtitle("Utility Differences (Changing Action Space)") + 
  ylab(expression(Delta)) + xlab("Number of Items") +
  scale_fill_manual(values=cbPalette) +
  facet_grid(cols = vars(num.items), switch="both") +
  theme(axis.text.x=element_blank(), axis.ticks.x = element_blank()) 

if(use_fitdata_for_cutoff == TRUE & appendix_version == TRUE) {
    plot_file <- file.path(plots_dir, "Appendix_Method_Compare_cutoff_from_fitdata.png")
} else if (use_fitdata_for_cutoff == TRUE & appendix_version == FALSE) {
    plot_file <- file.path(plots_dir, "Fig9_cutoff_from_fitdata.png")
} else if (use_fitdata_for_cutoff == FALSE & appendix_version == TRUE) {
    plot_file <- file.path(plots_dir, "Appendix_Method_Compare_cutoff_from_UQdata.png")
} else if (use_fitdata_for_cutoff == FALSE & appendix_version == FALSE) {
    plot_file <- file.path(plots_dir, "Fig9_cutoff_from_UQdata.png")
}

ggsave(file = plot_file, height = 3.25, width = 12, units = "in", dpi = 200)

