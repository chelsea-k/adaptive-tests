# Script for plotting Figure 10 and Table 4 from Krantsevich et al., 
# "BAYESIAN DECISION THEORY FOR TREE-BASED ADAPTIVE SCREENING TESTS 
#     WITH AN APPLICATION TO YOUTH DELINQUENCY"

source('adaptive_tests/code/util_functions.R')
library(caret)
library(tidyr)
library(dplyr)
library(cdata)
library(xtable)
library(ggplot2)
theme_set(theme_bw(base_size=14))

# parameters
maxIPP_list <- c(2:15)
w_list <- c(0.6)

# Set up folders and load data
data_test <- read.csv("simulated_data/item_response_data_test.csv")
data_test$y <- as.integer(data_test$y)
sub_rows <- which(data_test$Age>=15)

data_dir <- "output/out_of_sample/subpopulation/synthetic_data"
results_dir <- "output/out_of_sample/subpopulation/results"
plots_dir <- "output/plots"
tables_dir <- "output/tables"
dir.create(plots_dir, recursive = TRUE)
dir.create(tables_dir, recursive = TRUE)

synth_df_sub <- read.csv(file.path(data_dir, "synth_uncertainty_XB.csv"))
p_sub_synth <- read.csv(file.path(results_dir, "p_maxIPP_XB.synth_uncertainty_XB.csv"))
p_sub_test <- read.csv(file.path(results_dir, "p_maxIPP_XB.test.csv"))

synth_df_all <- read.csv("output/out_of_sample/all/synthetic_data/synth_uncertainty_XB.csv")
p_all_synth <- read.csv("output/out_of_sample/all/results/p_maxIPP_XB.synth_uncertainty_XB_other_pop.csv")
p_all_synth <- read.csv("output/out_of_sample/all/results/p_maxIPP_XB.synth_uncertainty_XB.csv")
p_all_test <- read.csv("output/out_of_sample/all/results/p_maxIPP_XB.test.csv")

results <- data.frame(matrix(NA, nrow=0, ncol=5))
colnames(results) <- c("maxIPP", "w", "Population", "Sensitivity", "Specificity")

for (i in seq_along(maxIPP_list)){
  maxIPP <- maxIPP_list[[i]]
  m_col <- which(colnames(p_sub_synth)==paste0("m.", maxIPP))
  p_sub_synth_temp <- p_sub_synth[,m_col]
  p_sub_test_temp <- p_sub_test[,m_col]
  p_all_synth_temp <- p_all_synth[,m_col]
  p_all_test_temp <- p_all_test[,m_col]
  
  cat(paste0("---------------- maxIPP = ", maxIPP, " ----------------\n"))
  for (j in seq_along(w_list)){
    w <- w_list[[j]]
    cutoff.sub <- get_cutoff(as.factor(synth_df_sub$y.draw), p_sub_synth_temp, w)
    y_sub <- as.integer(p_sub_test_temp >= cutoff.sub)
    cm_sub <- confusionMatrix(as.factor(y_sub[sub_rows]), 
                                as.factor(data_test$y[sub_rows]), positive="1")
    results <- rbind(results, data.frame(Num.Items=maxIPP, w=w, Population="Ages 15+",
                                        Sensitivity=cm_sub$byClass[['Sensitivity']], 
                                        Specificity=cm_sub$byClass[['Specificity']]))
    
    cutoff_all <- get_cutoff(as.factor(synth_df_sub$y.draw), p_all_synth_temp, w)
    y_all <- as.integer(p_all_test_temp >= cutoff_all)
    cm_all <- confusionMatrix(as.factor(y_all[sub_rows]), 
                                  as.factor(data_test$y[sub_rows]), positive="1")
    results <- rbind(results, data.frame(Num.Items=maxIPP, w=w, Population="All Youth",
                                        Sensitivity=cm_all$byClass[['Sensitivity']], 
                                        Specificity=cm_all$byClass[['Specificity']]))
    
  }
}
results$Utility <- (results$w)*(results$Sensitivity)+(1-results$w)*(results$Specificity)
write.csv(results, file.path(results_dir, "sens.spec.results.subpopulation.test.csv"), row.names=FALSE)

# Post-process results and create plot for Figure 10
results$Num.Items <- as.factor(results$Num.Items)
cbPalette <- c("#000000","#7a7a7a")

results_long <- gather(results, Quantity, Value, Sensitivity:Utility, factor_key=TRUE)
results_long_all <- results_long[which(results_long$Population=="All Youth"),] 
results_long_age <- results_long[which(results_long$Population=="Ages 15+"),]
results_long_diff <- results_long_all
results_long_diff$Value <- results_long_age$Value - results_long_all$Value
colnames(results_long_diff)[[5]] <- "Difference"
results_long_diff$Higher=ifelse(results_long_diff$Difference > 0,"Ages 15+","All Youth")

ggplot(results_long_diff, aes(x=Num.Items, y=Difference, fill=Higher)) + 
  geom_bar(stat = "identity") + 
  facet_grid(rows=vars(Quantity)) +
  scale_fill_manual(values=cbPalette) +
  xlab("Number of Items")  
ggsave(file.path(plots_dir,"Fig10.png"), height = 5, width = 7, units = "in", dpi = 200)

# Reformat results to be as in Table 4 
cT <- dplyr::tribble(
  ~Population, ~Sensitivity, ~Specificity, ~Utility,
  "Ages 15+",    "Sub_Sens", "Sub_Spec", "Sub_Util",
  "All Youth",    "All_Sens", "All_Spec", "All_Util"
)
results_wide <- blocks_to_rowrecs(results, cT, keyColumns = "Num.Items")
results_wide <- results_wide[,c(1,2,5,3,6,4,7)]

# Export table to latex format to create Table 4 (needs some reformatting of columns)
print(xtable(results_wide, type = "latex", digits=c(0,0,3,3,3,3,3,3)), 
      include.rownames=FALSE, file=file.path(tables_dir, "table4.tex"))
