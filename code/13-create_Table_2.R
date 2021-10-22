# Script for recreating Table 2 from Krantsevich et al., 
# "BAYESIAN DECISION THEORY FOR TREE-BASED ADAPTIVE SCREENING TESTS 
#     WITH AN APPLICATION TO YOUTH DELINQUENCY"

source('adaptive_tests/code/util_functions.R')
library(xtable)
library(tidyr)

# Set parameters
maxIPP_vals <- c(3,6,9,12,15)
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
p_maxIPP_XB_test <- read.csv(file.path(results_dir, "p_maxIPP_XB.test.csv"))

# Compute sens/spec/util results for all values of w and maxIPP
sens_spec_results <- data.frame(matrix(NA, nrow=0, ncol=4))
colnames(sens_spec_results) <- c("Number.of.Items", "w","Sensitivity", "Specificity")

for (t in seq_along(w_vals)){
  w <- w_vals[[t]]
  
  # Compute optimal cutoffs and predicted classes 
  maxIPP_XB_cut <- rep(NA, n_maxIPP)
  m_colnames <- paste0("m.", maxIPP_vals)
  m_cols <- which(colnames(p_maxIPP_XB_synth) %in% m_colnames)
  for (i in seq_along(maxIPP_vals)){
    maxIPP_XB_cut[[i]] <- get_cutoff(as.factor(y_synth), p_maxIPP_XB_synth[,m_cols[i]], w)
  }
  y_maxIPP_XB_test <- get_class(p_maxIPP_XB_test, maxIPP_XB_cut, maxIPP_vals, m_cols)
  
  # Compute sensitivity and specificity for all methods
  for (i in seq_along(maxIPP_vals)){
    maxIPP_XB_metrics <- get_utility(y_test,y.pred=y_maxIPP_XB_test[,i],type="class",w=w)
   
    # store results
    temp_results <- data.frame(Number.of.Items = maxIPP_vals[[i]],
                               w = w,
                               Sensitivity = maxIPP_XB_metrics$sens,
                               Specificity = maxIPP_XB_metrics$spec)
  
    sens_spec_results <- rbind(sens_spec_results, temp_results)
  }
}

# Post-processing table, pivoting results to wide format 
sens_spec_results <- sens_spec_results[!duplicated(sens_spec_results),]
sens_spec_results <- sens_spec_results[order(sens_spec_results$Number.of.Items),]
write.csv(sens_spec_results, file.path(results_dir, "sens_spec_results.all.test.csv"), row.names = FALSE)

sens_spec_results_wide <- pivot_wider(sens_spec_results,
                                      names_from = c("w"), 
                                      values_from = c("Sensitivity", "Specificity")) %>%
  select(Number.of.Items, 
         Sensitivity_0.4, Specificity_0.4,
         Sensitivity_0.5, Specificity_0.5,
         Sensitivity_0.6, Specificity_0.6) %>%
  as.data.frame()

# Exporting Table 2 to latex
table_name <- paste0("table2.tex")
print(xtable(sens_spec_results_wide, type = "latex", digits=c(rep(0,2), rep(3, 6))), 
      include.rownames=FALSE, file=file.path(tables_dir, table_name))
