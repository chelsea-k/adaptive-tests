# Script to grow regression and classification trees using maximum depth criterion

library(rpartMaxVPP) # can also use rpart for this one

out_of_sample = FALSE
subpopulation = FALSE

########################## Hyperparamters ##############################
CART_params = list(maxDepth_vals=c(2:15))
n_maxDepth = length(CART_params$maxDepth_vals)

########################## Data Preparation ##############################

# read in item response data
original_data_dir <- "preprocessed_original_data"
if(out_of_sample){
  data_train = read.csv(file.path(original_data_dir, "IMC_data_train_preprocessed.csv"))
  data_test = read.csv(file.path(original_data_dir, "IMC_data_test_preprocessed.csv"))
} else {
  data_train = read.csv(file.path(original_data_dir, "IMC_data_all_preprocessed.csv"))
  data_test = data_train
}

data_train$y <- as.integer(data_train$y)
data_test$y <- as.integer(data_test$y)

# read in generated data
cat("Loading synthetic data...\n")
synth_data_root <- "output"
folder = ifelse(out_of_sample, 
                file.path(synth_data_root, "out_of_sample"), 
                file.path(synth_data_root, "in_sample"))
folder =  file.path(folder, ifelse(subpopulation, "subpopulation", "all"))
synth_data_folder = file.path(folder, "synthetic_data")
synth_RF = read.csv(file.path(synth_data_folder, "synth_treefitting_RF.csv"))
synth_XB = read.csv(file.path(synth_data_folder, "synth_treefitting_XB.csv"))
synth_XB_util = read.csv(file.path(synth_data_folder, "synth_treefitting_XB_util_based_outcomes.csv"))
synth_XB_UQ = read.csv(file.path(synth_data_folder, "synth_uncertainty_XB.csv"))
n_synth = nrow(synth_XB)

# get column info for demo and item covariates
item_names = setdiff(colnames(data_train), c("y", "Age"))
RF_item_cols = which(colnames(synth_RF) %in% item_names)
XB_item_cols = which(colnames(synth_XB) %in% item_names)
real_item_cols = which(colnames(data_train) %in% item_names)

# set up storage folders
model_dir = file.path(folder, "model_fits_CART")
results_dir = file.path(folder, "results")
dir.create(model_dir, recursive = TRUE)
dir.create(results_dir, recursive = TRUE)


########################## Model Fitting ##############################

# Create dfs for saving RPART predictions; 
# "p." dfs are for regression trees, "y." dfs are for classification trees
cat("Creating dataframes for storing results...\n")
m_colnames = paste0("m.", CART_params$maxDepth_vals)

p_maxDepth_RF_synth = p_maxDepth_XB_synth = data.frame(matrix(NA, nrow=n_synth, ncol=n_maxDepth))
colnames(p_maxDepth_RF_synth) = colnames(p_maxDepth_XB_synth) = m_colnames

y_real_synth = y_RF_synth = y_XB_synth = y_util_synth = data.frame(matrix(NA, nrow=n_synth, ncol=n_maxDepth))
colnames(y_real_synth) = colnames(y_RF_synth) = colnames(y_XB_synth) = m_colnames
  
p_maxDepth_RF_test = p_maxDepth_XB_test = data.frame(matrix(NA, nrow=nrow(data_test), ncol=n_maxDepth))
colnames(p_maxDepth_RF_test) = colnames(p_maxDepth_XB_test) = m_colnames

y_real_test = y_RF_test = y_XB_test = y_util_test = data.frame(matrix(NA, nrow=nrow(data_test), ncol=n_maxDepth))
colnames(y_real_test) = colnames(y_RF_test) = colnames(y_XB_test) = m_colnames

# Define function for fitting CART via maxdepth (both regression and classification)
fit_CART_maxDepth = function(fitting_data, predict_data, test_data, 
                             item_cols, savefile, pred_df_synth, pred_df_test, 
                             tree_type, oos, params) {

maxDepth_vals = params$maxDepth_vals

for (i in seq_along(maxDepth_vals)){
  maxdepth=maxDepth_vals[[i]]
  m_col = which(colnames(pred_df_synth) == paste0("m.", maxdepth))
  cat(paste0("---------------- maxdepth = ", maxdepth, " ----------------\n"))
  if (tree_type == "classification") {
    tree_savefile <- paste0(savefile, ".classification_maxDepth.", maxdepth)
    if(file.exists(tree_savefile)) {
      cat("Loading classifcation tree \n")
      load(tree_savefile)
    } else {
      cat("Fitting classifcation tree \n")
      fit_rpart = rpartMaxVPP::rpart(y ~.,data = cbind(y=as.factor(fitting_data$y), fitting_data[,item_cols]), 
                        cp=0, method = "class", maxdepth=maxdepth, xval=0, maxcompete=0,
                        maxsurrogate=0, usesurrogate=0)
      save(fit_rpart, file=tree_savefile, ascii=TRUE)
    }
    cat("Predicting with classification tree \n")
    y_synth = predict(fit_rpart, predict_data, type="class")
    pred_df_synth[,m_col] = as.integer(as.character(y_synth))
    if(out_of_sample){
      y_test = predict(fit_rpart, test_data, type="class")
      pred_df_test[,m_col] = as.integer(as.character(y_test))
    }
  }
  
  if (tree_type == "regression") {
    tree_savefile <- paste0(savefile, ".regression_maxDepth.", maxdepth)
    if(file.exists(tree_savefile)){
      cat("Loading regression tree \n")
      load(tree_savefile)
    } else{
      cat("Fitting regression tree \n")
      fit_rpart = rpart(phat ~.,data = cbind(phat=fitting_data$phat, fitting_data[,item_cols]), 
                        cp=0, method = "anova", maxdepth=maxdepth, xval=0, maxcompete=0,
                        maxsurrogate=0, usesurrogate=0)
      save(fit_rpart, file=tree_savefile, ascii=TRUE)
    }
    cat("Predicting with regression tree \n")
    pred_df_synth[,m_col] = predict(fit_rpart, predict_data)
    if(out_of_sample){
      pred_df_test[,m_col] = predict(fit_rpart, test_data)
    }
  }
}
if (out_of_sample){
  results = list(pred_synth = pred_df_synth, pred_test = pred_df_test)
} else {
  results = list(pred_synth = pred_df_synth)
}
return(results)
}

# Run function to get maxDepth results and store results
cat("\n=== Fitting regression trees for RF synthetic data ===\n")
maxDepth_RF_results = fit_CART_maxDepth(fitting_data = synth_RF, predict_data = synth_XB_UQ,
                                        test_data = data_test, item_cols = RF_item_cols, 
                                        savefile = file.path(model_dir, "fit.CART.RF"),
                                        pred_df_synth = p_maxDepth_RF_synth, 
                                        pred_df_test = p_maxDepth_RF_test, 
                                        tree_type = "regression", oos = out_of_sample, 
                                        params = CART_params)
write.csv(maxDepth_RF_results$pred_synth, file.path(results_dir, "p_maxDepth_RF_synth.csv"), row.names = F)
if (out_of_sample) {
  write.csv(maxDepth_RF_results$pred_test, file.path(results_dir, "p_maxDepth_RF_test.csv"), row.names = F)
}

cat("\n=== Fitting regression trees for XBART synthetic data ===\n")
maxDepth_XB_results = fit_CART_maxDepth(fitting_data = synth_XB, predict_data = synth_XB_UQ, 
                                        test_data = data_test, item_cols = XB_item_cols, 
                                        savefile = file.path(model_dir, "fit.CART.XB"),
                                        pred_df_synth = p_maxDepth_XB_synth, 
                                        pred_df_test = p_maxDepth_XB_test, 
                                        tree_type = "regression", oos = out_of_sample, 
                                        params = CART_params)
write.csv(maxDepth_XB_results$pred_synth, file.path(results_dir, "p_maxDepth_XB_synth.csv"), row.names = F)
if (out_of_sample){
  write.csv(maxDepth_XB_results$pred_test, file.path(results_dir, "p_maxDepth_XB_test.csv"), row.names = F)
}

cat("\n=== Fitting classification trees for utility based outcome data ===\n")
class_util_results = fit_CART_maxDepth(fitting_data = synth_XB_util, predict_data = synth_XB_UQ,
                                       test_data = data_test, item_cols = XB_item_cols, 
                                       savefile = file.path(model_dir, "fit.CART.util"),
                                       pred_df_synth = y_real_synth, 
                                       pred_df_test = y_real_test, 
                                       tree_type = "classification", oos = out_of_sample, 
                                       params = CART_params)
write.csv(class_util_results$pred_synth, file.path(results_dir, "y_util_synth.csv"), row.names = F)
if (out_of_sample){
  write.csv(class_util_results$pred_test, file.path(results_dir, "y_util_test.csv"), row.names = F)
}

cat("\n=== Fitting classification trees for real IMC data ===\n")
class_real_results = fit_CART_maxDepth(fitting_data = data_train, predict_data = synth_XB_UQ,
                                       test_data = data_test, item_cols = real_item_cols, 
                                       savefile = file.path(model_dir, "fit.CART.real"),
                                       pred_df_synth = y_real_synth, 
                                       pred_df_test = y_real_test, 
                                       tree_type = "classification", oos = out_of_sample, 
                                       params = CART_params)
write.csv(class_real_results$pred_synth, file.path(results_dir, "y_real_synth.csv"), row.names = F)
if (out_of_sample){
  write.csv(class_real_results$pred_test, file.path(results_dir, "y_real_test.csv"), row.names = F)
}

cat("\n=== Fitting classification trees for RF synthetic data ===\n")
class_RF_results = fit_CART_maxDepth(fitting_data = synth_RF, predict_data = synth_XB_UQ,
                                     test_data = data_test, item_cols = RF_item_cols, 
                                     savefile = file.path(model_dir, "fit.CART.RF"),
                                     pred_df_synth = y_RF_synth, 
                                     pred_df_test = y_RF_test, 
                                     tree_type = "classification", oos = out_of_sample, 
                                     params = CART_params)
write.csv(class_RF_results$pred_synth, file.path(results_dir, "y_RF_synth.csv"), row.names = F)
if (out_of_sample){
  write.csv(class_RF_results$pred_test, file.path(results_dir, "y_RF_test.csv"), row.names = F)
}

cat("\n=== Fitting classification trees for XBART synthetic data ===\n\n")
class_XB_results = fit_CART_maxDepth(fitting_data = synth_XB, predict_data = synth_XB_UQ,
                                     test_data = data_test, item_cols = XB_item_cols, 
                                     savefile = file.path(model_dir, "fit.CART.XB"),
                                     pred_df_synth = y_XB_synth, 
                                     pred_df_test = y_XB_test, 
                                     tree_type = "classification", oos = out_of_sample, 
                                     params = CART_params)
write.csv(class_XB_results$pred_synth, file.path(results_dir, "y_XB_synth.csv"), row.names = F)
if (out_of_sample){
  write.csv(class_XB_results$pred_test, file.path(results_dir, "y_XB_test.csv"), row.names = F)
}

