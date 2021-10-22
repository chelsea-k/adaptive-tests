# Script for recreating Tables in Appendix F from Krantsevich et al., 
# "BAYESIAN DECISION THEORY FOR TREE-BASED ADAPTIVE SCREENING TESTS 
#     WITH AN APPLICATION TO YOUTH DELINQUENCY"

source('adaptive_tests/code/util_functions.R')
library(xtable)

# Set parameters
maxIPP_vals <- c(2:15)
n_maxIPP <- length(maxIPP_vals)
w_vals <- c(0.4, 0.5, 0.6)

# Set up folders and load data
data_dir <- "output/out_of_sample/all/synthetic_data"
results_dir <- "output/out_of_sample/all/results"
tables_dir <- "output/tables"
dir.create(tables_dir, recursive = TRUE)

data_test <- read.csv("simulated_data/item_response_data_test.csv")
y_test <- as.integer(data_test$y)

synth_df_XB <- read.csv(file.path(data_dir, "synth_uncertainty_XB.csv"))
y_synth <- synth_df_XB$y.draw

p_maxIPP_XB_synth <- read.csv(file.path(results_dir, "p_maxIPP_XB.synth_uncertainty_XB.csv"))
p_maxIPP_RF_synth <- read.csv(file.path(results_dir, "p_maxIPP_RF.synth_uncertainty_XB.csv"))
p_maxdepth_XB_synth <- read.csv(file.path(results_dir, "p_maxDepth_XB.synth_uncertainty_XB.csv"))
p_maxdepth_RF_synth <- read.csv(file.path(results_dir, "p_maxDepth_RF.synth_uncertainty_XB.csv"))
#y_class_real_synth <- read.csv(file.path(results_dir, "y_real.synth_uncertainty_XB.csv"))
y_class_RF_synth <- read.csv(file.path(results_dir, "y_RF.synth_uncertainty_XB.csv"))
y_class_XB_synth <- read.csv(file.path(results_dir, "y_XB.synth_uncertainty_XB.csv"))
y_class_XB_util <- read.csv(file.path(results_dir, "y_util.synth_uncertainty_XB.csv"))

p_maxIPP_XB_test <- read.csv(file.path(results_dir, "p_maxIPP_XB.test.csv"))
p_maxIPP_RF_test <- read.csv(file.path(results_dir, "p_maxIPP_RF.test.csv"))
p_maxdepth_XB_test <- read.csv(file.path(results_dir, "p_maxDepth_XB.test.csv"))
p_maxdepth_RF_test <- read.csv(file.path(results_dir, "p_maxDepth_RF.test.csv"))
#y_class_real_test <- read.csv(file.path(results_dir, "y_real.test.csv"))
y_class_RF_test <- read.csv(file.path(results_dir, "y_RF.test.csv"))
y_class_XB_test <- read.csv(file.path(results_dir, "y_XB.test.csv"))
y_class_util_test_0.4 <- read.csv(file.path(results_dir, "y_util.test_w_0.4.csv"))
y_class_util_test_0.5 <- read.csv(file.path(results_dir, "y_util.test_w_0.5.csv"))
y_class_util_test_0.6 <- read.csv(file.path(results_dir, "y_util.test_w_0.6.csv"))



n_methods <- 7

# Compute sens/spec/util results for all values of w and maxIPP
sens_spec_results <- data.frame(matrix(NA, nrow=0, ncol=7))
colnames(sens_spec_results) <- c("Items", "Tree.Type", "Criterion", "w", 
                                 "Fitting.Data", "Sensitivity", "Specificity")

for (t in seq_along(w_vals)){
  w <- w_vals[[t]]
  
  # set up extra utility method
  if(w == 0.4) {  
    y_class_util_test <- y_class_util_test_0.4
  } else if (w == 0.5) {
    y_class_util_test <- y_class_util_test_0.5
  } else if (w == 0.6) {
    y_class_util_test <- y_class_util_test_0.6
  }

  # Compute optimal cutoffs and predicted classes 
  maxIPP_XB_cut <- rep(NA, n_maxIPP)
  maxIPP_RF_cut <- rep(NA, n_maxIPP)
  maxdepth_XB_cut <- rep(NA, n_maxIPP)
  maxdepth_RF_cut <- rep(NA, n_maxIPP)
  m_colnames <- paste0("m.", maxIPP_vals)
  m_cols <- which(colnames(p_maxIPP_XB_synth) %in% m_colnames)
  for (i in seq_along(maxIPP_vals)){
    maxIPP_XB_cut[[i]] <- get_cutoff(as.factor(y_synth), p_maxIPP_XB_synth[,m_cols[i]], w)
    maxIPP_RF_cut[[i]] <- get_cutoff(as.factor(y_synth), p_maxIPP_RF_synth[,m_cols[i]], w)
    maxdepth_XB_cut[[i]] <- get_cutoff(as.factor(y_synth), p_maxdepth_XB_synth[,m_cols[i]], w)
    maxdepth_RF_cut[[i]] <- get_cutoff(as.factor(y_synth), p_maxdepth_RF_synth[,m_cols[i]], w)
  }
  y_maxIPP_XB_test <- get_class(p_maxIPP_XB_test, maxIPP_XB_cut, maxIPP_vals, m_cols)
  y_maxIPP_RF_test <- get_class(p_maxIPP_RF_test, maxIPP_RF_cut, maxIPP_vals, m_cols)
  y_maxDepth_XB_test <- get_class(p_maxdepth_XB_test, maxdepth_XB_cut, maxIPP_vals, m_cols)
  y_maxDepth_RF_test <- get_class(p_maxdepth_RF_test, maxdepth_RF_cut, maxIPP_vals, m_cols)
  
  # Extract columns of maxIPP vals being used
  #y_class_real_test <- y_class_real_test[,which(colnames(y_class_real_test) %in% m_colnames)]
  y_class_RF_test <- y_class_RF_test[,which(colnames(y_class_RF_test) %in% m_colnames)]
  y_class_XB_test <- y_class_XB_test[,which(colnames(y_class_XB_test) %in% m_colnames)]
  y_class_util_test <- y_class_util_test[,which(colnames(y_class_util_test) %in% m_colnames)]
  
  # Compute sensitivity and specificity for all methods
  for (i in seq_along(maxIPP_vals)){
    maxIPP_XB_metrics <- get_utility(y_test,y.pred=y_maxIPP_XB_test[,i],type="class",w=w)
    maxdepth_XB_metrics <- get_utility(y_test,y.pred=y_maxDepth_XB_test[,i],type="class",w=w)
    maxIPP_RF_metrics <- get_utility(y_test,y.pred=y_maxIPP_RF_test[,i],type="class",w=w)
    maxdepth_RF_metrics <- get_utility(y_test,y.pred=y_maxDepth_RF_test[,i],type="class",w=w)
    class_XB_metrics  <- get_utility(y_test,y.pred=y_class_XB_test[,i],type="class",w=w)
    class_RF_metrics  <- get_utility(y_test,y.pred=y_class_RF_test[,i],type="class",w=w)
    #class_real_metrics  <- get_utility(y_test,y.pred=y_class_real_test[,i],type="class",w=w)
    class_util_metrics  <- get_utility(y_test,y.pred=y_class_util_test[,i],type="class",w=w)
    
    temp_results <- data.frame(Items = rep(maxIPP_vals[[i]], n_methods),
                               Tree.Type = c("Regression + Cutoff", 
                                             "Regression + Cutoff", 
                                             "Regression + Cutoff", 
                                             "Regression + Cutoff", 
                                             "Classification", 
                                             "Classification",
                                             "Classification"),
                               Criterion = c("maxIPP", 
                                             "maxDepth", 
                                             "maxIPP", 
                                             "maxDepth", 
                                             "maxDepth",
                                             "maxDepth", 
                                             "maxDepth"),
                               w = c(rep(w,5), rep("-", 2)),
                               Fitting.Data = c("GCFM + XBART", 
                                                "GCFM + XBART", 
                                                "Perturb + RF", 
                                                "Perturb + RF",
                                                "GCFM + Utility", 
                                                "GCFM + XBART", 
                                                "Perturb + RF"),
                               Sensitivity = c(maxIPP_XB_metrics$sens, 
                                        maxdepth_XB_metrics$sens, 
                                        maxIPP_RF_metrics$sens, 
                                        maxdepth_RF_metrics$sens,
                                        class_util_metrics$sens,
                                        class_XB_metrics$sens, 
                                        class_RF_metrics$sens),
                               Specificity = c(maxIPP_XB_metrics$spec, 
                                        maxdepth_XB_metrics$spec, 
                                        maxIPP_RF_metrics$spec, 
                                        maxdepth_RF_metrics$spec,
                                        class_util_metrics$spec,
                                        class_XB_metrics$spec, 
                                        class_RF_metrics$spec)
                              )
  
    sens_spec_results <- rbind(sens_spec_results, temp_results)
  }
}

# Post-processing and exporting Table from Appendix F to latex
sens_spec_results <- sens_spec_results[!duplicated(sens_spec_results),]
sens_spec_results <- sens_spec_results[order(sens_spec_results$w,
                                             sens_spec_results$Items, 
                                             desc(sens_spec_results$Tree.Type),
                                             desc(sens_spec_results$Fitting.Data)),]
write.csv(sens_spec_results, file.path(results_dir, "sens_spec_results.all.test.csv"), row.names = FALSE)

table_name <- paste0("table_AppendixF.tex")
print(xtable(sens_spec_results, type = "latex", digits=c(0,0,0,0,2,0,3,3)), 
      include.rownames=FALSE, file= file.path(tables_dir, table_name))

