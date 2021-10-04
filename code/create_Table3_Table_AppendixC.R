# Script for recreating Table 3 & Table in Appendix C from Krantsevich et al., 
# "BAYESIAN DECISION THEORY FOR TREE-BASED ADAPTIVE SCREENING TESTS 
#     WITH AN APPLICATION TO YOUTH DELINQUENCY"

source('adaptive_tests/code/util_functions.R')
library(xtable)

# Set parameters
which_table <- "table3" # "table3" or "app.C"
if(which_table == "table3"){
  maxIPP_vals <- c(3,9,15)
} else if (which_table == "app.C") {
  maxIPP_vals <- c(2:15)
} else if (which_table == "app.D") {
  maxIPP_vals <- c(2:15)
} 

n_maxIPP <- length(maxIPP_vals)
w_vals <- c(0.25,0.5,0.75)

# Set up folders and load data
data_dir <- "output/out_of_sample/all/synthetic_data"
results_dir <- "output/out_of_sample/all/results"
tables_dir <- "output/tables"
dir.create(tables_dir, recursive = TRUE)

data_test <- read.csv("preprocessed_original_data/IMC_data_test_preprocessed.csv")
y_test <- as.integer(data_test$y)

synth_df_XB <- read.csv(file.path(data_dir, "synth_uncertainty_XB.csv"))
y_synth <- synth_df_XB$y.draw

p_maxIPP_XB_synth <- read.csv(file.path(results_dir, "p_maxIPP_XB.synth_uncertainty_XB.csv"))
p_maxIPP_RF_synth <- read.csv(file.path(results_dir, "p_maxIPP_RF.synth_uncertainty_XB.csv"))
p_maxdepth_XB_synth <- read.csv(file.path(results_dir, "p_maxDepth_XB.synth_uncertainty_XB.csv"))
p_maxdepth_RF_synth <- read.csv(file.path(results_dir, "p_maxDepth_RF.synth_uncertainty_XB.csv"))
y_class_real_synth <- read.csv(file.path(results_dir, "y_real.synth_uncertainty_XB.csv"))
y_class_RF_synth <- read.csv(file.path(results_dir, "y_RF.synth_uncertainty_XB.csv"))
y_class_XB_synth <- read.csv(file.path(results_dir, "y_XB.synth_uncertainty_XB.csv"))
y_class_XB_util <- read.csv(file.path(results_dir, "y_util.synth_uncertainty_XB.csv"))

p_maxIPP_XB_test <- read.csv(file.path(results_dir, "p_maxIPP_XB.test.csv"))
p_maxIPP_RF_test <- read.csv(file.path(results_dir, "p_maxIPP_RF.test.csv"))
p_maxdepth_XB_test <- read.csv(file.path(results_dir, "p_maxDepth_XB.test.csv"))
p_maxdepth_RF_test <- read.csv(file.path(results_dir, "p_maxDepth_RF.test.csv"))
y_class_real_test <- read.csv(file.path(results_dir, "y_real.test.csv"))
y_class_RF_test <- read.csv(file.path(results_dir, "y_RF.test.csv"))
y_class_XB_test <- read.csv(file.path(results_dir, "y_XB.test.csv"))
y_class_util_test <- read.csv(file.path(results_dir, "y_util.test.csv"))


n_methods <- ifelse(which_table == "app.D", 2, 7)

# Compute sens/spec/util results for all values of w and maxIPP
sens_spec_results <- data.frame(matrix(NA, nrow=0, ncol=7))
colnames(sens_spec_results) <- c("Num.Items", "Tree.Type", "Stop.Criterion", "w", 
                               "Fitting.Data", "Sensitivity", "Specificity")

for (t in seq_along(w_vals)){
  w <- w_vals[[t]]
  
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
  y_class_real_test <- y_class_real_test[,which(colnames(y_class_real_test) %in% m_colnames)]
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
    class_real_metrics  <- get_utility(y_test,y.pred=y_class_real_test[,i],type="class",w=w)
    class_util_metrics  <- get_utility(y_test,y.pred=y_class_util_test[,i],type="class",w=w)
    
    if(which_table == "app.D") {
      temp_results <- data.frame(Num.Items = rep(maxIPP_vals[[i]], n_methods),
                                Tree.Type = c("Regression", "Classification"),
                                Stop.Criterion = c("maxIPP", "maxDepth"),
                                w = c(w, "NA"),
                                Fitting.Data = c("BFA + XBART", "Util Outcomes"),
                                Sensitivity = c(maxIPP_XB_metrics$sens, class_util_metrics$sens),
                                Specificity = c(maxIPP_XB_metrics$spec, class_util_metrics$spec)
      )
      
    } else {
        temp_results <- data.frame(Num.Items = rep(maxIPP_vals[[i]], n_methods),
                                  Tree.Type = c("Regression", 
                                                "Regression", 
                                                "Regression", 
                                                "Regression", 
                                                "Classification", 
                                                "Classification", 
                                                "Classification"),
                                  Stop.Criterion = c("maxIPP", 
                                                     "maxDepth", 
                                                     "maxIPP", 
                                                     "maxDepth", 
                                                     "maxDepth", 
                                                     "maxDepth", 
                                                     "maxDepth"),
                                  w = c(rep(w, 4), rep("NA", 3)),
                                  Fitting.Data = c("BFA + XBART", 
                                                   "BFA + XBART", 
                                                   "Perturb + RF", 
                                                   "Perturb + RF",
                                                   "BFA + XBART", 
                                                   "Perturb + RF", 
                                                   "Real Train Data"),
                                  Sensitivity = c(maxIPP_XB_metrics$sens, 
                                                  maxdepth_XB_metrics$sens, 
                                                  maxIPP_RF_metrics$sens, 
                                                  maxdepth_RF_metrics$sens,
                                                  class_XB_metrics$sens, 
                                                  class_RF_metrics$sens,
                                                  class_real_metrics$sens),
                                  Specificity = c(maxIPP_XB_metrics$spec, 
                                                  maxdepth_XB_metrics$spec, 
                                                  maxIPP_RF_metrics$spec, 
                                                  maxdepth_RF_metrics$spec,
                                                  class_XB_metrics$spec, 
                                                  class_RF_metrics$spec,
                                                  class_real_metrics$spec)
                                  )
    }
    sens_spec_results <- rbind(sens_spec_results, temp_results)
    }
}

# Post-processing and exporting Table 3 and Table from Appendix C to latex
sens_spec_results <- sens_spec_results[!duplicated(sens_spec_results),]
sens_spec_results <- sens_spec_results[order(sens_spec_results$Num.Items, sens_spec_results$w),]
write.csv(sens_spec_results, file.path(results_dir, "sens_spec_results.all.test.csv"), row.names = FALSE)

if (which_table == "table3") {
  sens_spec_results_summary <- sens_spec_results[which(!(sens_spec_results$Stop.Criterion == "maxDepth" & 
                                                          sens_spec_results$Tree.Type == "Regression")),]
  print(xtable(sens_spec_results_summary, type = "latex", digits=c(0,0,0,0,2,0,3,3)), 
        include.rownames=FALSE, file=file.path(tables_dir, "table3.tex"))
}else if(which_table == "app.C") {
    print(xtable(sens_spec_results, type = "latex", digits=c(0,0,0,0,2,0,3,3)), 
          include.rownames=FALSE, file= file.path(tables_dir, "table_AppendixC.tex"))
} else if(which_table == "app.D") {
  print(xtable(sens_spec_results, type = "latex", digits=c(0,0,0,0,2,0,3,3)), 
        include.rownames=FALSE, file= file.path(tables_dir, "table_AppendixD.tex"))
}
