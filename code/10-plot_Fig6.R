# Script for plotting Figure 6 from Krantsevich et al., 
# "BAYESIAN DECISION THEORY FOR TREE-BASED ADAPTIVE SCREENING TESTS 
#     WITH AN APPLICATION TO YOUTH DELINQUENCY"

source("adaptive_tests/code/util_functions.R")  # get_cutoff(), get_utility()
library(ggplot2)
theme_set(theme_bw(base_size=14))

# Set parameters
maxIPP_vals <- c(2:15)
n_maxIPP <- length(maxIPP_vals)
w <- 0.5

# Set up data folders and load data
data_dir <- "output/in_sample/subpopulation/synthetic_data"
results_dir <- "output/in_sample/subpopulation/results"
plots_dir <- "output/plots"
dir.create(plots_dir, recursive = TRUE)

synth_sub_uncertainty <- read.csv(file.path(data_dir, "synth_uncertainty_XB.csv"))
p_CART_sub_uncertainty <- read.csv(file.path(results_dir, "p_maxIPP_XB.synth_uncertainty_XB.csv")) 
synth_all_uncertainty <- read.csv("output/in_sample/all/synthetic_data/synth_uncertainty_XB.csv")
p_CART_all_uncertainty <- read.csv("output/in_sample/all/results/p_maxIPP_XB.synth_uncertainty_XB.csv")

# hyperparameters
n_mcmc <- length(unique(synth_sub_uncertainty$post.idx))
n_rows_total <- n_maxIPP*n_mcmc

# Compute optimal cutoffs and predicted classes 
XB_cut_all <- get_cutoff(as.factor(synth_all_uncertainty$y.draw), synth_all_uncertainty$phat.mean, w)
XB_cut_sub <- get_cutoff(as.factor(synth_sub_uncertainty$y.draw), synth_sub_uncertainty$phat.mean, w)
maxIPP_cut_all <- rep(NA, n_maxIPP)
maxIPP_cut_sub <- rep(NA, n_maxIPP)
m_cols <- which(colnames(p_CART_sub_uncertainty) %in% paste0("m.", maxIPP_vals))
for (i in seq_along(maxIPP_vals)){
  maxIPP_cut_all[[i]] <- get_cutoff(as.factor(synth_all_uncertainty$y.draw), p_CART_all_uncertainty[,m_cols[[i]]], w)
  maxIPP_cut_sub[[i]] <- get_cutoff(as.factor(synth_sub_uncertainty$y.draw), p_CART_sub_uncertainty[,m_cols[[i]]], w)
}

# Get class labels
y_XB_all <- 1*(synth_all_uncertainty$phat.mean >= XB_cut_all)
y_XB_sub <- 1*(synth_sub_uncertainty$phat.mean >= XB_cut_sub)
y_CART_all <- get_class(p_CART_all_uncertainty, maxIPP_cut_all, maxIPP_vals, m_cols)
y_CART_sub <- get_class(p_CART_sub_uncertainty, maxIPP_cut_sub, maxIPP_vals, m_cols)

# Compute delta results for all values of maxIPP for both populations
delta_draws_all <- data.frame(matrix(NA, nrow=n_rows_total, ncol=4))
colnames(delta_draws_all) <- c("maxIPP", "post.idx", "Delta", "Population")

for (j in 1:n_mcmc){
  if (j%%100==0) {cat(paste0("On iteration j = ", j, " out of ", n_mcmc, "\n"))}
  post_rows <- which(synth_sub_uncertainty$post.idx==j)
  temp_df_all <- synth_all_uncertainty[post_rows,]
  temp_df_sub <- synth_sub_uncertainty[post_rows,]
  XB_metrics_all <- get_utility(temp_df_all$y.draw,y.pred=y_XB_all[post_rows],type="class",w=w) 
  XB_metrics_sub <- get_utility(temp_df_sub$y.draw,y.pred=y_XB_sub[post_rows],type="class",w=w) 
  for (i in seq_along(maxIPP_vals)){
    maxIPP <- maxIPP_vals[[i]]
    m_col_all <- which(colnames(p_CART_all_uncertainty)==paste0("m.", maxIPP))
    m_col_sub <- which(colnames(p_CART_sub_uncertainty)==paste0("m.", maxIPP))
    CART_metrics_all <- get_utility(temp_df_all$y.draw,y.pred=y_CART_all[post_rows,m_col_all],
                                   type="class",w=w) 
    CART_metrics_sub <- get_utility(temp_df_sub$y.draw,y.pred=y_CART_sub[post_rows,m_col_sub],
                                   type="class",w=w) 
    temp_results <- data.frame(maxIPP = maxIPP, post.idx = j, 
                              Delta = c(CART_metrics_all$util - XB_metrics_all$util, 
                                        CART_metrics_sub$util - XB_metrics_sub$util),
                              Population=c("All", "Ages 15+"))
    delta_draws_all[c(1:2) + (2*(i-1)) + (2*n_maxIPP*(j-1)),] = temp_results
  }
}

# Post-process results and create plot for Figure 6
delta_draws_all$maxIPP <- as.factor(delta_draws_all$maxIPP)
delta_draws_all$Population <- factor(delta_draws_all$Population, levels=c("All", "Ages 15+"))
greyPalette <- c("#808080","#D3D3D3")

ggplot(delta_draws_all, aes(x=maxIPP, y=Delta, fill=Population)) +
  geom_boxplot(outlier.shape=NA, position=position_dodge2(padding=0.2, width=10)) + 
  theme(plot.title = element_text(hjust = 0.5))  + 
  ggtitle("Utility Differences (Changing Target Population)") + 
  ylab(expression(Delta)) + xlab("Number of Items") +
  scale_fill_manual(values=greyPalette) +
  guides(fill=guide_legend(title="Population"))
ggsave(file.path(plots_dir, "Fig6.png"), height = 3.25, width = 8, units = "in", dpi = 200)


