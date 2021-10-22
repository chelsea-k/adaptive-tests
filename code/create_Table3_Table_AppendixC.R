# Script for recreating Table 3 & Table in Appendix C from Krantsevich et al., 
# "BAYESIAN DECISION THEORY FOR TREE-BASED ADAPTIVE SCREENING TESTS 
#     WITH AN APPLICATION TO YOUTH DELINQUENCY"

source('adaptive_tests/code/util_functions.R')
library(xtable)
library(tidyr)

# Set parameters
maxIPP_vals <- c(3,9,15)
n_maxIPP <- length(maxIPP_vals)
w_vals <- c(0.4, 0.5, 0.6)

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
y_class_RF_synth <- read.csv(file.path(results_dir, "y_RF.synth_uncertainty_XB.csv"))
y_class_XB_synth <- read.csv(file.path(results_dir, "y_XB.synth_uncertainty_XB.csv"))

p_maxIPP_XB_test <- read.csv(file.path(results_dir, "p_maxIPP_XB.test.csv"))
p_maxIPP_RF_test <- read.csv(file.path(results_dir, "p_maxIPP_RF.test.csv"))
y_class_RF_test <- read.csv(file.path(results_dir, "y_RF.test.csv"))
y_class_XB_test <- read.csv(file.path(results_dir, "y_XB.test.csv"))

n_methods <- 4 

# Compute sens/spec/util results for all values of w and maxIPP
sens_spec_results <- data.frame(matrix(NA, nrow=0, ncol=7))
colnames(sens_spec_results) <- c("Items", "Tree.Type",  "Criterion", "w",
                               "Fitting.Data", "Sens", "Spec")

for (t in seq_along(w_vals)){
  w <- w_vals[[t]]
  
  # Compute optimal cutoffs and predicted classes 
  maxIPP_XB_cut <- rep(NA, n_maxIPP)
  maxIPP_RF_cut <- rep(NA, n_maxIPP)
  m_colnames <- paste0("m.", maxIPP_vals)
  m_cols <- which(colnames(p_maxIPP_XB_synth) %in% m_colnames)
  for (i in seq_along(maxIPP_vals)){
    maxIPP_XB_cut[[i]] <- get_cutoff(as.factor(y_synth), p_maxIPP_XB_synth[,m_cols[i]], w)
    maxIPP_RF_cut[[i]] <- get_cutoff(as.factor(y_synth), p_maxIPP_RF_synth[,m_cols[i]], w)
  }
  y_maxIPP_XB_test <- get_class(p_maxIPP_XB_test, maxIPP_XB_cut, maxIPP_vals, m_cols)
  y_maxIPP_RF_test <- get_class(p_maxIPP_RF_test, maxIPP_RF_cut, maxIPP_vals, m_cols)

  # Extract columns of maxIPP vals being used
  y_class_RF_test <- y_class_RF_test[,which(colnames(y_class_RF_test) %in% m_colnames)]
  y_class_XB_test <- y_class_XB_test[,which(colnames(y_class_XB_test) %in% m_colnames)]

  # Compute sensitivity and specificity for all methods
  for (i in seq_along(maxIPP_vals)){
    maxIPP_XB_metrics <- get_utility(y_test,y.pred=y_maxIPP_XB_test[,i],type="class",w=w)
    maxIPP_RF_metrics <- get_utility(y_test,y.pred=y_maxIPP_RF_test[,i],type="class",w=w)
    class_XB_metrics  <- get_utility(y_test,y.pred=y_class_XB_test[,i],type="class",w=w)
    class_RF_metrics  <- get_utility(y_test,y.pred=y_class_RF_test[,i],type="class",w=w)

    # store results
    temp_results <- data.frame(Items = rep(maxIPP_vals[[i]], n_methods),
                              Tree.Type = c("Regression", 
                                            "Regression", 
                                            "Classification", 
                                            "Classification"),
                              Criterion = c("maxIPP", 
                                            "maxIPP", 
                                            "maxDepth", 
                                            "maxDepth"),
                              w = c(rep(w, 4)),
                              Fitting.Data = c("GCFM + XB", 
                                               "Perturb + RF",
                                               "GCFM + XB", 
                                               "Perturb + RF"),
                              Sens = c(maxIPP_XB_metrics$sens, 
                                              maxIPP_RF_metrics$sens, 
                                              class_XB_metrics$sens, 
                                              class_RF_metrics$sens),
                              Spec = c(maxIPP_XB_metrics$spec, 
                                              maxIPP_RF_metrics$spec, 
                                              class_XB_metrics$spec, 
                                              class_RF_metrics$spec)
                              )
    sens_spec_results <- rbind(sens_spec_results, temp_results)
    }
}

# Post-processing and exporting Table 2 to latex
sens_spec_results <- sens_spec_results[!duplicated(sens_spec_results),]
sens_spec_results <- sens_spec_results[order(sens_spec_results$Items),]
write.csv(sens_spec_results, file.path(results_dir, "sens_spec_results.all.test.csv"), row.names = FALSE)

sens_spec_results_wide <- pivot_wider(sens_spec_results,
                                      names_from = c("w", "Tree.Type"), 
                                      values_from = c("Sens", "Spec")) %>%
  mutate(Tree.Type = ifelse(Criterion == "maxIPP", "Regression", "Classification")) %>%
  select(-c(Sens_0.4_Classification, Spec_0.4_Classification,
            Spec_0.4_Classification, Spec_0.6_Classification)) %>%
  rename(Sens_Classification = Sens_0.5_Classification,
         Spec_Classification = Spec_0.5_Classification) %>%
  select(Items, Tree.Type, Fitting.Data, 
         Sens_Classification, Spec_Classification, 
         Sens_0.4_Regression, Spec_0.4_Regression,
         Sens_0.5_Regression, Spec_0.5_Regression,
         Sens_0.6_Regression, Spec_0.6_Regression,) %>%
  as.data.frame()

table_name <- paste0("table2.tex")
print(xtable(sens_spec_results_wide, type = "latex", digits=c(rep(0,4), rep(3, 8))), 
      include.rownames=FALSE, file=file.path(tables_dir, table_name))
