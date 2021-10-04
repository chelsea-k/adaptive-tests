# Script for plotting Figure 5 from Krantsevich et al., 
# "BAYESIAN DECISION THEORY FOR TREE-BASED ADAPTIVE SCREENING TESTS 
#     WITH AN APPLICATION TO YOUTH DELINQUENCY"

source("adaptive_tests/code/util_functions.R")  # get_metrics()
library(ggplot2)
library(caret)
theme_set(theme_bw(base_size=14))

# Set parameters
maxIPP <- 3
w_list <- c(0.25, 0.5, 0.75)
post_idx_list <- c(100, 200, 300, 400)

# Set up data folders and load item response data, cutoff dataframe
data_dir <- "output/in_sample/all/synthetic_data"
results_dir <- "output/in_sample/all/results"
plots_dir <- "output/plots"
dir.create(plots_dir, recursive = TRUE)

synth_df_uncertainty <- read.csv(file.path(data_dir, "synth_uncertainty_XB.csv"))
p_CART_df_uncertainty <- read.csv(file.path(results_dir, "p_maxIPP_XB.synth_uncertainty_XB.csv")) 
n_mcmc <- length(unique(synth_df_uncertainty$post.idx))
cutoffs <- read.csv(file.path(results_dir, "cutoffs.csv"))

# Compute sensitivity/specificity for XBART, CART ROC curves
#   and (spec, sens) pairs for cutoff-specific points 

metrics_XBART <- data.frame()
metrics_CART <- data.frame()
ROC_points <- data.frame()

for (t in seq_along(post_idx_list)) {
  post_idx <- post_idx_list[[t]]
  ROC_samp <- which(synth_df_uncertainty$post.idx == post_idx)
  col <- which(colnames(p_CART_df_uncertainty)==paste0("m.",maxIPP))
  p_CART <- p_CART_df_uncertainty[ROC_samp,col]
  y_class <- synth_df_uncertainty[ROC_samp, "y.draw"]
  p_XBART <- synth_df_uncertainty[ROC_samp, "phat.mean"]
  
  metrics_XBART_temp <- get_metrics(p_XBART, y_class)
  metrics_XBART <- rbind(metrics_XBART, cbind(group=0, post_idx = post_idx, 
                              Model="XBART", metrics_XBART_temp))
  metrics_CART_temp <- get_metrics(p_CART, y_class)
  metrics_CART <- rbind(metrics_CART, cbind(group=0, post_idx = post_idx, 
                              Model="maxIPP=3", metrics_CART_temp))
  
  ROC_points_temp <- data.frame(matrix(0, nrow=2*length(w_list), ncol=5))
  colnames(ROC_points_temp) <- c("Specificity", "Sensitivity", "w", "post_idx", "Model")
  for (k in 1:length(w_list)){
    w <- w_list[[k]]
    ct_row <- which(cutoffs$w==w)
    p_XBART_factor <- factor(1*(p_XBART>=cutoffs[ct_row, "XBART"]), levels=c("0", "1"))
    y_class_factor <- factor(y_class, levels=c("0", "1"))
    cm_XBART <- confusionMatrix(p_XBART_factor,y_class_factor, positive="1")
    ROC_points_temp[2*k-1,] <- list(cm_XBART$byClass[['Specificity']], 
                          cm_XBART$byClass[['Sensitivity']], w, post_idx, "XBART")
    CART_col <- which(colnames(cutoffs)==paste0("CART.m.", maxIPP))
    p_CART_factor <- factor(1*(p_CART>=cutoffs[ct_row,CART_col]), levels=c("0", "1"))
    cm_CART <- confusionMatrix(p_CART_factor, y_class_factor, positive="1")
    ROC_points_temp[2*k,] <- list(cm_CART$byClass[['Specificity']], 
                          cm_CART$byClass[['Sensitivity']], w, post_idx, "CART")
  }
  ROC_points <- rbind(ROC_points, ROC_points_temp)
}

# Post-process results and create plot for Figure 5
ROC_points$w <- as.factor(ROC_points$w)
ROC_points$post_idx <- as.integer(ROC_points$post_idx)
cbPalette <- c("#009E73", "#0072B2", "#D55E00")

ggplot(ROC_points, aes(x=Specificity, y=Sensitivity, group=w)) +
  xlim(c(1,0)) +
  geom_line(data=metrics_XBART, color="black", 
            aes(x=Specificity, y=Sensitivity, linetype="XBART")) +
  geom_line(data=metrics_CART, color="black", 
            aes(x=Specificity, y=Sensitivity, linetype=paste0("maxIPP=", maxIPP,""))) +
  geom_point(size=3.5,  aes(color=w, shape=w)) +
  scale_linetype_manual(name="Action", values=c("dashed","solid")) +
  scale_colour_manual(values = cbPalette) + 
  scale_shape_manual(values=c(16,15,17)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_grid(cols=vars(post_idx), labeller = label_both) 

ggsave(file.path(plots_dir,"Fig5.png"), height = 3, width = 10, units = "in", dpi = 250)





