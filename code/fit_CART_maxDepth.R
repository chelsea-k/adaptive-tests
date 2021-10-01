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
folder_root = ifelse(out_of_sample, 
                file.path(synth_data_root, "out_of_sample"), 
                file.path(synth_data_root, "in_sample"))
folder =  file.path(folder_root, ifelse(subpopulation, "subpopulation", "all"))
synth_data_folder = file.path(folder, "synthetic_data")
synth_treefitting_RF = read.csv(file.path(synth_data_folder, "synth_treefitting_RF.csv"))
synth_treefitting_XB = read.csv(file.path(synth_data_folder, "synth_treefitting_XB.csv"))
synth_treefitting_XB_util = read.csv(file.path(synth_data_folder, "synth_treefitting_XB_util_based_outcomes.csv"))
synth_uncertainty_XB = read.csv(file.path(synth_data_folder, "synth_uncertainty_XB.csv"))


folder_other_pop =  file.path(folder_root, ifelse(subpopulation, "all", "subpopulation"))
synth_data_folder_other_pop = file.path(folder_other_pop, "synthetic_data")
synth_uncertainty_XB_other_pop = read.csv(file.path(synth_data_folder_other_pop, "synth_uncertainty_XB.csv"))
n_synth = nrow(synth_treefitting_XB)
n_real = nrow(data_train)

# get item column info
item_names = setdiff(colnames(data_train), c("y", "Age"))

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

# set up results df for RF regression
p_maxDepth_RF.synth_treefitting_RF = data.frame(matrix(NA, nrow=n_synth, ncol=n_maxDepth))
p_maxDepth_RF.synth_uncertainty_XB = data.frame(matrix(NA, nrow=n_synth, ncol=n_maxDepth))
p_maxDepth_RF.synth_uncertainty_XB_other_pop = data.frame(matrix(NA, nrow=n_synth, ncol=n_maxDepth))

colnames(p_maxDepth_RF.synth_treefitting_RF) = m_colnames
colnames(p_maxDepth_RF.synth_uncertainty_XB) = m_colnames
colnames(p_maxDepth_RF.synth_uncertainty_XB_other_pop) = m_colnames

# set up results for XB regression
p_maxDepth_XB.synth_treefitting_XB = data.frame(matrix(NA, nrow=n_synth, ncol=n_maxDepth))
p_maxDepth_XB.synth_uncertainty_XB = data.frame(matrix(NA, nrow=n_synth, ncol=n_maxDepth))
p_maxDepth_XB.synth_uncertainty_XB_other_pop = data.frame(matrix(NA, nrow=n_synth, ncol=n_maxDepth))

colnames(p_maxDepth_XB.synth_treefitting_XB) = m_colnames
colnames(p_maxDepth_XB.synth_uncertainty_XB) = m_colnames
colnames(p_maxDepth_XB.synth_uncertainty_XB_other_pop) = m_colnames


# set up results for real data classification
y_real.real_data = data.frame(matrix(NA, nrow=n_real, ncol=n_maxDepth))
y_real.synth_uncertainty_XB = data.frame(matrix(NA, nrow=n_synth, ncol=n_maxDepth))
y_real.synth_uncertainty_XB_other_pop = data.frame(matrix(NA, nrow=n_synth, ncol=n_maxDepth))

colnames(y_real.real_data) = m_colnames
colnames(y_real.synth_uncertainty_XB) = m_colnames
colnames(y_real.synth_uncertainty_XB_other_pop) = m_colnames

# set up results for utility-based classification
y_util.synth_treefitting_XB_util = data.frame(matrix(NA, nrow=n_synth, ncol=n_maxDepth))
y_util.synth_uncertainty_XB = data.frame(matrix(NA, nrow=n_synth, ncol=n_maxDepth))
y_util.synth_uncertainty_XB_other_pop = data.frame(matrix(NA, nrow=n_synth, ncol=n_maxDepth))

colnames(y_util.synth_treefitting_XB_util) = m_colnames
colnames(y_util.synth_uncertainty_XB) = m_colnames
colnames(y_util.synth_uncertainty_XB_other_pop) = m_colnames

# set up results for RF classification
y_RF.synth_treefitting_RF = data.frame(matrix(NA, nrow=n_synth, ncol=n_maxDepth))
y_RF.synth_uncertainty_XB = data.frame(matrix(NA, nrow=n_synth, ncol=n_maxDepth))
y_RF.synth_uncertainty_XB_other_pop = data.frame(matrix(NA, nrow=n_synth, ncol=n_maxDepth))

colnames(y_RF.synth_treefitting_RF) = m_colnames
colnames(y_RF.synth_uncertainty_XB) = m_colnames
colnames(y_RF.synth_uncertainty_XB_other_pop) = m_colnames

# set up results for XB classification
y_XB.synth_treefitting_XB = data.frame(matrix(NA, nrow=n_synth, ncol=n_maxDepth))
y_XB.synth_uncertainty_XB = data.frame(matrix(NA, nrow=n_synth, ncol=n_maxDepth))
y_XB.synth_uncertainty_XB_other_pop = data.frame(matrix(NA, nrow=n_synth, ncol=n_maxDepth))

colnames(y_XB.synth_treefitting_XB) = m_colnames
colnames(y_XB.synth_uncertainty_XB) = m_colnames
colnames(y_XB.synth_uncertainty_XB_other_pop) = m_colnames


# set up results for test data
  
p_maxDepth.RF_test = p_maxDepth.XB_test = data.frame(matrix(NA, nrow=nrow(data_test), ncol=n_maxDepth))
colnames(p_maxDepth.RF_test) = colnames(p_maxDepth.XB_test) = m_colnames

y_real.test = y_RF.test = y_XB.test = y_util.test = data.frame(matrix(NA, nrow=nrow(data_test), ncol=n_maxDepth))
colnames(y_real.test) = colnames(y_RF.test) = colnames(y_XB.test) = colnames(y_util.test) = m_colnames

# Define function for fitting CART via maxdepth (both regression and classification)
fit_CART_maxDepth = function(fitting_data, 
                             predict_uncertainty_data,
                             predict_uncertainty_data_other_pop,
                             test_data, 
                             item_names,
                             savefile, 
                             pred_df_treefitting, 
                             pred_df_synth_uncertainty,
                             pred_df_synth_uncertainty_other_pop, 
                             pred_df_test, 
                             tree_type, 
                             oos, 
                             params) {

maxDepth_vals = params$maxDepth_vals

for (i in seq_along(maxDepth_vals)){
  maxdepth=maxDepth_vals[[i]]
  m_col = which(colnames(pred_df_treefitting) == paste0("m.", maxdepth))
  cat(paste0("---------------- maxdepth = ", maxdepth, " ----------------\n"))
  if (tree_type == "classification") {
    tree_savefile <- paste0(savefile, ".classification_maxDepth.", maxdepth)
    if(file.exists(tree_savefile)) {
      cat("Loading classifcation tree \n")
      load(tree_savefile)
    } else {
      cat("Fitting classifcation tree \n")
      item_cols = which(colnames(fitting_data) %in% item_names)
      fit_rpart = rpartMaxVPP::rpart(y ~.,data = cbind(y=as.factor(fitting_data$y), fitting_data[,item_cols]), 
                        cp=0, method = "class", maxdepth=maxdepth, xval=0, maxcompete=0,
                        maxsurrogate=0, usesurrogate=0)
	  #fit_rpart$where <- NULL
	  #fit_rpart$y <- NULL
      #save(fit_rpart, file=tree_savefile, ascii=TRUE)
    }
    cat("Predicting with classification tree \n")
    y_hat = predict(fit_rpart, fitting_data, type="class")
    pred_df_treefitting[,m_col] = as.integer(as.character(y_hat))
    
    y_hat = predict(fit_rpart, predict_uncertainty_data, type="class")
    pred_df_synth_uncertainty[,m_col] = as.integer(as.character(y_hat))
    
    y_hat = predict(fit_rpart, predict_uncertainty_data_other_pop, type="class")
    pred_df_synth_uncertainty_other_pop[,m_col] = as.integer(as.character(y_hat))
    
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
      item_cols = which(colnames(fitting_data) %in% item_names)
      fit_rpart = rpart(phat ~.,data = cbind(phat=fitting_data$phat, fitting_data[,item_cols]), 
                        cp=0, method = "anova", maxdepth=maxdepth, xval=0, maxcompete=0,
                        maxsurrogate=0, usesurrogate=0)
	  #fit_rpart$y <- NULL
	  #fit_rpart$where <- NULL
      #save(fit_rpart, file=tree_savefile, ascii=TRUE)
    }
    cat("Predicting with regression tree \n")
    pred_df_treefitting[,m_col] = predict(fit_rpart, fitting_data)
    pred_df_synth_uncertainty[,m_col] = predict(fit_rpart, predict_uncertainty_data)
    pred_df_synth_uncertainty_other_pop[,m_col] = predict(fit_rpart, predict_uncertainty_data_other_pop)
    if(out_of_sample){
      pred_df_test[,m_col] = predict(fit_rpart, test_data)
    }
  }
}
results = list(pred_synth_treefitting = pred_df_treefitting, 
               pred_synth_uncertainty = pred_df_synth_uncertainty, 
               pred_synth_uncertainty_other_pop = pred_df_synth_uncertainty_other_pop)

if (out_of_sample == TRUE) {
  results$pred_test = pred_df_test
}

return(results)
}


# Run function to get maxDepth results and store results
cat("\n=== Fitting regression trees for RF synthetic data ===\n")
maxDepth_RF_results = fit_CART_maxDepth(fitting_data = synth_treefitting_RF, 
                                        predict_uncertainty_data = synth_uncertainty_XB,
                                        predict_uncertainty_data_other_pop = synth_uncertainty_XB_other_pop, 
                                        test_data = data_test, 
                                        item_names = item_names, 
                                        savefile = file.path(model_dir, "fit.CART.RF"),
                                        pred_df_treefitting = p_maxDepth_RF.synth_treefitting_RF, 
                                        pred_df_synth_uncertainty = p_maxDepth_RF.synth_uncertainty_XB,
                                        pred_df_synth_uncertainty_other_pop = p_maxDepth_RF.synth_uncertainty_XB_other_pop,
                                        pred_df_test = p_maxDepth.RF_test, 
                                        tree_type = "regression", 
                                        oos = out_of_sample, 
                                        params = CART_params)
write.csv(maxDepth_RF_results$pred_synth_treefitting, file.path(results_dir, "p_maxDepth_RF.synth_treefitting_RF.csv"), row.names = F)
write.csv(maxDepth_RF_results$pred_synth_uncertainty, file.path(results_dir, "p_maxDepth_RF.synth_uncertainty_XB.csv"), row.names = F)
write.csv(maxDepth_RF_results$pred_synth_uncertainty_other_pop, file.path(results_dir, "p_maxDepth_RF.synth_uncertainty_XB_other_pop.csv"), row.names = F)
if (out_of_sample) {
  write.csv(maxDepth_RF_results$pred_test, file.path(results_dir, "p_maxDepth_RF.test.csv"), row.names = F)
}


cat("\n=== Fitting regression trees for XBART synthetic data ===\n")
maxDepth_XB_results = fit_CART_maxDepth(fitting_data = synth_treefitting_XB, 
                                        predict_uncertainty_data = synth_uncertainty_XB,
                                        predict_uncertainty_data_other_pop = synth_uncertainty_XB_other_pop, 
                                        test_data = data_test, 
                                        item_names = item_names,  
                                        savefile = file.path(model_dir, "fit.CART.XB"),
                                        pred_df_treefitting = p_maxDepth_XB.synth_treefitting_XB, 
                                        pred_df_synth_uncertainty = p_maxDepth_XB.synth_uncertainty_XB,
                                        pred_df_synth_uncertainty_other_pop = p_maxDepth_XB.synth_uncertainty_XB_other_pop,
                                        pred_df_test = p_maxDepth.XB_test, 
                                        tree_type = "regression", 
                                        oos = out_of_sample, 
                                        params = CART_params)
write.csv(maxDepth_XB_results$pred_synth_treefitting, file.path(results_dir, "p_maxDepth_XB.synth_treefitting_XB.csv"), row.names = F)
write.csv(maxDepth_XB_results$pred_synth_uncertainty, file.path(results_dir, "p_maxDepth_XB.synth_uncertainty_XB.csv"), row.names = F)
write.csv(maxDepth_XB_results$pred_synth_uncertainty_other_pop, file.path(results_dir, "p_maxDepth_XB.synth_uncertainty_XB_other_pop.csv"), row.names = F)
if (out_of_sample){
  write.csv(maxDepth_XB_results$pred_test, file.path(results_dir, "p_maxDepth_XB.test.csv"), row.names = F)
}


cat("\n=== Fitting classification trees for real IMC data ===\n")
class_real_results = fit_CART_maxDepth(fitting_data = data_train, 
                                       predict_uncertainty_data = synth_uncertainty_XB,
                                       predict_uncertainty_data_other_pop = synth_uncertainty_XB_other_pop, 
                                       test_data = data_test, 
                                       item_names = item_names, 
                                       savefile = file.path(model_dir, "fit.CART.real"),
                                       pred_df_treefitting = y_real.real_data,
                                       pred_df_synth_uncertainty = y_real.synth_uncertainty_XB,
                                       pred_df_synth_uncertainty_other_pop = y_real.synth_uncertainty_XB_other_pop,
                                       pred_df_test = y_real.test,
                                       tree_type = "classification", 
                                       oos = out_of_sample,
                                       params = CART_params)
write.csv(class_real_results$pred_synth_treefitting, file.path(results_dir, "y_real.real_data.csv"), row.names = F)
write.csv(class_real_results$pred_synth_uncertainty, file.path(results_dir, "y_real.synth_uncertainty_XB.csv"), row.names = F)
write.csv(class_real_results$pred_synth_uncertainty_other_pop, file.path(results_dir, "y_real.synth_uncertainty_XB_other_pop.csv"), row.names = F)

if (out_of_sample){
 write.csv(class_real_results$pred_test, file.path(results_dir, "y_real.test.csv"), row.names = F)
}


cat("\n=== Fitting classification trees for utility based outcome data ===\n")
class_util_results = fit_CART_maxDepth(fitting_data = synth_treefitting_XB_util, 
                                       predict_uncertainty_data = synth_uncertainty_XB,
                                       predict_uncertainty_data_other_pop = synth_uncertainty_XB_other_pop, 
                                       test_data = data_test, 
                                       item_names = item_names, 
                                       savefile = file.path(model_dir, "fit.CART.util"),
                                       pred_df_treefitting = y_util.synth_treefitting_XB_util,
                                       pred_df_synth_uncertainty = y_util.synth_uncertainty_XB,
                                       pred_df_synth_uncertainty_other_pop = y_util.synth_uncertainty_XB_other_pop,
                                       pred_df_test = y_util.test,
                                       tree_type = "classification", 
                                       oos = out_of_sample,
                                       params = CART_params)
write.csv(class_util_results$pred_synth_treefitting, file.path(results_dir, "y_util.synth_treefitting_XB_util.csv"), row.names = F)
write.csv(class_util_results$pred_synth_uncertainty, file.path(results_dir, "y_util.synth_uncertainty_XB.csv"), row.names = F)
write.csv(class_util_results$pred_synth_uncertainty_other_pop, file.path(results_dir, "y_util.synth_uncertainty_XB_other_pop.csv"), row.names = F)

if (out_of_sample){
  write.csv(class_util_results$pred_test, file.path(results_dir, "y_util.test.csv"), row.names = F)
}

cat("\n=== Fitting classification trees for RF synthetic data ===\n")
class_RF_results = fit_CART_maxDepth(fitting_data = synth_treefitting_RF, 
                                     predict_uncertainty_data = synth_uncertainty_XB,
                                     predict_uncertainty_data_other_pop = synth_uncertainty_XB_other_pop, 
                                     test_data = data_test, 
                                     item_names = item_names, 
                                     savefile = file.path(model_dir, "fit.CART.RF"),
                                     pred_df_treefitting = y_RF.synth_treefitting_RF,
                                     pred_df_synth_uncertainty = y_RF.synth_uncertainty_XB,
                                     pred_df_synth_uncertainty_other_pop = y_RF.synth_uncertainty_XB_other_pop,
                                     pred_df_test = y_RF.test,
                                     tree_type = "classification", oos = out_of_sample,
                                     params = CART_params)
write.csv(class_RF_results$pred_synth_treefitting, file.path(results_dir, "y_RF.synth_treefitting_RF.csv"), row.names = F)
write.csv(class_RF_results$pred_synth_uncertainty, file.path(results_dir, "y_RF.synth_uncertainty_XB.csv"), row.names = F)
write.csv(class_RF_results$pred_synth_uncertainty_other_pop, file.path(results_dir, "y_RF.synth_uncertainty_XB_other_pop.csv"), row.names = F)
if (out_of_sample){
 write.csv(class_RF_results$pred_test, file.path(results_dir, "y_RF.test.csv"), row.names = F)
}

cat("\n=== Fitting classification trees for XBART synthetic data ===\n\n")
class_XB_results = fit_CART_maxDepth(fitting_data = synth_treefitting_XB, 
                                     predict_uncertainty_data = synth_uncertainty_XB,
                                     predict_uncertainty_data_other_pop = synth_uncertainty_XB_other_pop, 
                                     test_data = data_test, 
                                     item_names = item_names, 
                                     savefile = file.path(model_dir, "fit.CART.XB"),
                                     pred_df_treefitting = y_XB.synth_treefitting_XB,
                                     pred_df_synth_uncertainty = y_XB.synth_uncertainty_XB,
                                     pred_df_synth_uncertainty_other_pop = y_XB.synth_uncertainty_XB_other_pop,
                                     pred_df_test = y_XB.test,
                                     tree_type = "classification", 
                                     oos = out_of_sample,
                                     params = CART_params)
write.csv(class_XB_results$pred_synth_treefitting, file.path(results_dir, "y_XB.synth_treefitting_XB.csv"), row.names = F)
write.csv(class_XB_results$pred_synth_uncertainty, file.path(results_dir, "y_XB.synth_uncertainty_XB.csv"), row.names = F)
write.csv(class_XB_results$pred_synth_uncertainty_other_pop, file.path(results_dir, "y_XB.synth_uncertainty_XB_other_pop.csv"), row.names = F)
if (out_of_sample){
 write.csv(class_XB_results$pred_test, file.path(results_dir, "y_XB.test.csv"), row.names = F)
}

