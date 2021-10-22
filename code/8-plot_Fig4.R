# Script for plotting Figure 4 from Krantsevich et al., 
# "BAYESIAN DECISION THEORY FOR TREE-BASED ADAPTIVE SCREENING TESTS 
#     WITH AN APPLICATION TO YOUTH DELINQUENCY"

source("adaptive_tests/code/util_functions.R")  # get_cutoff(), get_utility()
library(ggplot2)
library(dplyr)
theme_set(theme_bw(base_size=14))

# Set parameters
maxIPP_vals <- c(2:15)
n_maxIPP <- length(maxIPP_vals)
w_vals <- c(0.25, 0.5, 0.75, 0.4, 0.6)
w_vals_plot <- c(0.25, 0.5, 0.75)

# Set up data folders and load data
data_dir <- "output/in_sample/all/synthetic_data"
results_dir <- "output/in_sample/all/results"
plots_dir <- "output/plots"
dir.create(plots_dir, recursive = TRUE)

synth_df_uncertainty <- read.csv(file.path(data_dir, "synth_uncertainty_XB.csv"))
p_CART_df_uncertainty <- read.csv(file.path(results_dir, "p_maxIPP_XB.synth_uncertainty_XB.csv")) 
n_mcmc <- length(unique(synth_df_uncertainty$post.idx))
n_rows_total <- n_maxIPP*n_mcmc

# Dataframe for storing results
delta_draws_all <- data.frame(matrix(NA, nrow=3*n_rows_total, ncol=4))
colnames(delta_draws_all) <- c("maxIPP", "post.idx", "Delta", "w")
cutoffs <- data.frame(matrix(NA, nrow=length(w_vals), ncol=2+n_maxIPP))
colnames(cutoffs)<-c("w", "XBART", paste0("CART.m.",maxIPP_vals))

# Compute delta results for all values of w and maxIPP
for (k in seq_along(w_vals)){
  cat(paste0("-------------- Predicting for cutoff #", k, " -------------------\n"))
  w <- w_vals[[k]]  
  
  # compute optimal cutoffs
  XB_cutoff <- get_cutoff(as.factor(synth_df_uncertainty$y.draw), synth_df_uncertainty$phat.mean, w)
  maxIPP_cutoff <- rep(NA, n_maxIPP)
  for (i in seq_along(maxIPP_vals)){
    maxIPP <- maxIPP_vals[[i]]
    col <- which(colnames(p_CART_df_uncertainty)==paste0("m.", maxIPP))
    maxIPP_cutoff[[i]] <- get_cutoff(as.factor(synth_df_uncertainty$y.draw), p_CART_df_uncertainty[,col], w)
  }
  cutoffs[k,] <- c(w, XB_cutoff, maxIPP_cutoff)
  
  # compute delta draws
  delta_draws <- data.frame(matrix(NA, nrow=n_rows_total, ncol=3))
  colnames(delta_draws) <- c("maxIPP", "post.idx", "Delta")
  
  for (j in 1:n_mcmc){
    if (j%%100==0) {cat(paste0("On iteration j = ", j, " out of ", n_mcmc, "\n"))}
    post_rows <- which(synth_df_uncertainty$post.idx==j)
    temp_df <- synth_df_uncertainty[post_rows,]
    XB_metrics <- get_utility(temp_df$y.draw,prob=temp_df$phat.mean,cut=XB_cutoff,type="prob",w=w) 
    for (i in seq_along(maxIPP_vals)){
      temp_df2 <- p_CART_df_uncertainty[post_rows,]
      maxIPP <- maxIPP_vals[[i]]
      col <- which(colnames(p_CART_df_uncertainty)==paste0("m.", maxIPP))
      CT_metrics <- get_utility(temp_df$y.draw, prob=temp_df2[,col], 
                               cut=maxIPP_cutoff[[i]], type="prob", w=w) 
      delta_draws[n_maxIPP*(j-1)+i,] <- data.frame(maxIPP=maxIPP, 
                                                   post.idx=j, 
                                                   Delta=CT_metrics$util-XB_metrics$util)
    }
  }
  delta_draws_all[(k-1)*n_rows_total+(1:n_rows_total),] <- cbind(delta_draws, w=w)
}
write.csv(delta_draws_all, file.path(results_dir, "delta_draws_all.csv"), row.names = FALSE)
write.csv(cutoffs, file.path(results_dir, "cutoffs.csv"), row.names = FALSE)


# Post-process results and create plot for Figure 4
delta_draws_plot <- delta_draws_all
delta_draws_plot <- delta_draws_plot %>%
  filter(w %in% w_vals_plot) %>%
  mutate(maxIPP = as.factor(maxIPP), w = as.factor(w))
cbPalette <- c("#009E73", "#0072B2", "#D55E00")

ggplot(delta_draws_plot, aes(x=maxIPP, y=Delta, fill=w)) +
 geom_boxplot(outlier.shape=NA, position=position_dodge2(padding=0.2, width=1.2)) + 
 theme(plot.title = element_text(hjust = 0.5))  + 
 ggtitle("Utility Differences (Changing Utility Function)") + ylab(expression(Delta)) +
 scale_fill_manual(values=cbPalette) +
 ylim(c(-0.15,0.05)) +
 xlab("Number of Items")  

ggsave(file.path(plots_dir, "Fig4.png"), height = 3.25, width = 10, units = "in", dpi = 200)

