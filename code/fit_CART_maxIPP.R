# Script to fit regression trees using maxIPP growing method and 
# pruning based on out-of-sample RMSE

library(rpartMaxVPP)

out_of_sample = TRUE
subpopulation = FALSE

########################## Hyperparamters ##############################
CART_params = list(minbucket=100,
                   maxIPP_vals=c(2:15),
                   cp=1e-7,
                   patience = 10)
n_maxIPP = length(CART_params$maxIPP_vals)

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
#synth_treefitting_RF = read.csv(file.path(synth_data_folder, "synth_treefitting_RF.csv"))
synth_treefitting_XB = read.csv(file.path(synth_data_folder, "synth_treefitting_XB.csv"))
synth_uncertainty_XB = read.csv(file.path(synth_data_folder, "synth_uncertainty_XB.csv"))
#prune_RF = read.csv(file.path(synth_data_folder, "prune_RF.csv"))
prune_treefitting_XB = read.csv(file.path(synth_data_folder, "prune_treefitting_XB.csv"))
prune_uncertainty_XB = read.csv(file.path(synth_data_folder, "prune_uncertainty_XB.csv"))

# read in synthetic data from other population
folder_other_pop =  file.path(folder_root, ifelse(subpopulation, "all", "subpopulation"))
synth_data_folder_other_pop = file.path(folder_other_pop, "synthetic_data")
synth_uncertainty_XB_other_pop = read.csv(file.path(synth_data_folder_other_pop, "synth_uncertainty_XB.csv"))
n_synth = nrow(synth_treefitting_XB)
n_real = nrow(data_train)

# get column info for demo and item covariates
item_names = setdiff(colnames(data_train), c("y", "Age"))

# set up storage folders
model_dir = file.path(folder, "model_fits_CART/")
results_dir = file.path(folder, "results/")
dir.create(model_dir, recursive = TRUE)
dir.create(results_dir, recursive = TRUE)

########################## Model Fitting ##############################

# set up results df for RF regression
cat("Creating dataframes for storing results...\n")
m_colnames = paste0("m.", CART_params$maxIPP_vals)

p_maxIPP_RF.synth_treefitting_RF = data.frame(matrix(NA, nrow=n_synth, ncol=n_maxIPP))
p_maxIPP_RF.synth_uncertainty_XB = data.frame(matrix(NA, nrow=n_synth, ncol=n_maxIPP))
p_maxIPP_RF.synth_uncertainty_XB_other_pop = data.frame(matrix(NA, nrow=n_synth, ncol=n_maxIPP))

colnames(p_maxIPP_RF.synth_treefitting_RF) = m_colnames
colnames(p_maxIPP_RF.synth_uncertainty_XB) = m_colnames
colnames(p_maxIPP_RF.synth_uncertainty_XB_other_pop) = m_colnames

# set up results for XB regression
p_maxIPP_XB.synth_treefitting_XB = data.frame(matrix(NA, nrow=n_synth, ncol=n_maxIPP))
p_maxIPP_XB.synth_uncertainty_XB = data.frame(matrix(NA, nrow=n_synth, ncol=n_maxIPP))
p_maxIPP_XB.synth_uncertainty_XB_other_pop = data.frame(matrix(NA, nrow=n_synth, ncol=n_maxIPP))

colnames(p_maxIPP_XB.synth_treefitting_XB) = m_colnames
colnames(p_maxIPP_XB.synth_uncertainty_XB) = m_colnames
colnames(p_maxIPP_XB.synth_uncertainty_XB_other_pop) = m_colnames

p_maxIPP_RF.test = p_maxIPP_XB.test = data.frame(matrix(NA, nrow=nrow(data_test), ncol=n_maxIPP))
colnames(p_maxIPP_RF.test) = colnames(p_maxIPP_XB.test) = paste0("m.", CART_params$maxIPP_vals)

# Define function for fitting CART via maxIPP and pruning
fit_CART_maxIPP = function(fitting_data, 
                           prune_data, 
                           predict_uncertainty_data,
                           predict_uncertainty_data_other_pop, 
                           test_data, 
                           item_names, 
                           savefile, 
                           pred_df_treefitting, 
                           pred_df_synth_uncertainty,
                           pred_df_synth_uncertainty_other_pop, 
                           pred_df_test, 
                           oos, 
                           params, 
						               save_all_trees = FALSE,
						               save_model = FALSE) {
 
maxIPP_vals = params$maxIPP_vals   
cp = params$cp
minbucket = params$minbucket
patience = params$patience

for (i in seq_along(maxIPP_vals)){
  maxIPP=maxIPP_vals[[i]]
  cat(paste0("---------------- maxIPP = ", maxIPP, " ----------------\n"))
  
  # check if tree already exists
  tree_savefile <- paste0(savefile, ".regression_maxIPP.", maxIPP)
  opt_tree_savefile <- paste0(savefile, ".regression_maxIPP.", maxIPP, ".opt.fitboth")
  if(!file.exists(opt_tree_savefile)) {
    if(file.exists(tree_savefile)) {
      cat("Loading big regression tree... \n")
      load(tree_savefile)
    } else {
        cat("Fitting big regression tree... \n")
        # fit the large CART tree
        item_cols = which(colnames(fitting_data) %in% item_names)
        fit_rpart = rpartMaxVPP::rpart(p ~., data = cbind(p=fitting_data$phat.mean, fitting_data[,item_cols]),
                       method = "anova", cp=cp, maxdepth=25, minbucket = minbucket, maxvpp=maxIPP,
                       xval=0, maxcompete=0, maxsurrogate=0, usesurrogate=0, model=TRUE)
        if(save_all_trees == TRUE){
          if(save_model == FALSE) {fit_rpart$model <- NULL}
          environment(fit_rpart$terms) <- NULL
          fit_rpart$where <- NULL
          save(fit_rpart, file=opt_tree_savefile, ascii=TRUE)
        }
    }
    
    # find optimal pruning point based on threshold/patience
    cat("Pruning regression tree... \n")
    cplist = fit_rpart$cptable
    threshold =  ifelse(maxIPP<5, 1e-4, 1e-5)
    opt_cp = cplist[1,"CP"]
    root_tree = prune(fit_rpart, cp=opt_cp)
    p_CART = predict(root_tree, newdata = prune_data[,item_cols], type='vector')
    rmse_prev = sqrt(mean((p_CART - prune_data$phat.mean)^2))
    wait=0
    if(nrow(cplist) == 1) {
      opt_tree = fit_rpart
    } else {
        for (j in 2:nrow(cplist)) {
          if (j%%50==0) {cat("Predicting for cp: ", j, "out of", nrow(cplist),"\n")}
          this_cp = cplist[j,"CP"]
          pruned_tree = prune(fit_rpart, cp=this_cp)
          p_CART = predict(pruned_tree, newdata = prune_data[,item_cols], type='vector')
          rmse = sqrt(mean((p_CART - prune_data$phat.mean)^2))
          rmse_diff = rmse_prev - rmse
          if (rmse_diff > threshold){
              rmse_prev = rmse
              opt_cp = this_cp
              wait=0
          } else {
              wait = wait + 1
              rmse_prev = rmse
              if (wait >= patience){break}
          }  
        }
      
       # Prune to optimal cp and store optimally pruned tree
        opt_tree = prune(fit_rpart, cp=opt_cp)
    }
    if(save_all_trees == TRUE){
      if(save_model == FALSE) {opt_tree$model <- NULL}
      environment(opt_tree$terms) <- NULL
      opt_tree$where <- NULL
      save(opt_tree, file=opt_tree_savefile, ascii=TRUE)
    }
  } else {
    cat("Loading optimally pruned regression tree \n")
    load(opt_tree_savefile)
  }
  
  # Store CART predictions from pruned tree on synthetic XBART data
  opt_p_CART = predict(opt_tree, fitting_data)
  m_col = which(colnames(pred_df_treefitting) == paste0("m.", maxIPP))
  pred_df_treefitting[,m_col] = opt_p_CART
  
  opt_p_CART = predict(opt_tree, predict_uncertainty_data)
  m_col = which(colnames(pred_df_synth_uncertainty) == paste0("m.", maxIPP))
  pred_df_synth_uncertainty[,m_col] = opt_p_CART
  
  opt_p_CART = predict(opt_tree, predict_uncertainty_data_other_pop)
  m_col = which(colnames(pred_df_synth_uncertainty_other_pop) == paste0("m.", maxIPP))
  pred_df_synth_uncertainty_other_pop[,m_col] = opt_p_CART
  
  
  if(oos){
    # Store CART predictions from pruned tree on test data
    opt_p_test = predict(opt_tree, test_data)
    m_col = which(colnames(pred_df_test) == paste0("m.", maxIPP))
    pred_df_test[,m_col] = opt_p_test
  }
  
}

results = list(pred_synth_treefitting = pred_df_treefitting, 
               pred_synth_uncertainty = pred_df_synth_uncertainty, 
               pred_synth_uncertainty_other_pop = pred_df_synth_uncertainty_other_pop)

if (out_of_sample){
  results$pred_test <- pred_df_test
} 

return(results)
}

# Run function to get maxIPP results and store results
#maxIPP_RF_results = fit_CART_maxIPP(fitting_data = synth_treefitting_RF, 
#                                    prune_data = prune_RF,
#                                    predict_uncertainty_data = synth_uncertainty_XB,
#                                    predict_uncertainty_data_other_pop = synth_uncertainty_XB_other_pop, 
#                                    test_data = data_test,
#                                    item_names = item_names,
#                                    savefile = file.path(model_dir, "fit.CART.RF"),
#                                    pred_df_treefitting = p_maxIPP_RF.synth_treefitting_RF, 
#                                    pred_df_synth_uncertainty = p_maxIPP_RF.synth_uncertainty_XB,
#                                    pred_df_synth_uncertainty_other_pop = p_maxIPP_RF.synth_uncertainty_XB_other_pop, 
#                                    pred_df_test = p_maxIPP_RF.test,
#                                    oos = out_of_sample, 
#                                    params = CART_params)
#write.csv(maxIPP_RF_results$pred_synth_treefitting, file.path(results_dir, "p_maxIPP_RF.synth_treefitting_RF.csv"), row.names = F)
#write.csv(maxIPP_RF_results$pred_synth_uncertainty, file.path(results_dir, "p_maxIPP_RF.synth_uncertainty_XB.csv"), row.names = F)
#write.csv(maxIPP_RF_results$pred_synth_uncertainty_other_pop, file.path(results_dir, "p_maxIPP_RF.synth_uncertainty_XB_other_pop.csv"), row.names = F)
#if (out_of_sample) {
#  write.csv(maxIPP_RF_results$pred_test, file.path(results_dir, "p_maxIPP_RF.test.csv"), row.names = F)
#}

maxIPP_XB_results = fit_CART_maxIPP(fitting_data = synth_treefitting_XB,
                                    prune_data = prune_treefitting_XB, 
                                    predict_uncertainty_data = synth_uncertainty_XB,
                                    predict_uncertainty_data_other_pop = synth_uncertainty_XB_other_pop,
                                    test_data = data_test, 
                                    item_names = item_names, 
                                    savefile = file.path(model_dir, "fit.CART.XB"),
                                    pred_df_treefitting = p_maxIPP_XB.synth_treefitting_XB, 
                                    pred_df_synth_uncertainty = p_maxIPP_XB.synth_uncertainty_XB,
                                    pred_df_synth_uncertainty_other_pop = p_maxIPP_XB.synth_uncertainty_XB_other_pop, 
                                    pred_df_test = p_maxIPP_XB.test, 
                                    oos = out_of_sample, 
                                    params = CART_params,
									                  save_all_trees = TRUE,
									                  save_model = FALSE)

write.csv(maxIPP_XB_results$pred_synth_treefitting, file.path(results_dir, "p_maxIPP_XB.synth_treefitting_XB.csv"), row.names = F)
write.csv(maxIPP_XB_results$pred_synth_uncertainty, file.path(results_dir, "p_maxIPP_XB.synth_uncertainty_XB.csv"), row.names = F)
write.csv(maxIPP_XB_results$pred_synth_uncertainty_other_pop, file.path(results_dir, "p_maxIPP_XB.synth_uncertainty_XB_other_pop.csv"), row.names = F)
if (out_of_sample){
  write.csv(maxIPP_XB_results$pred_test, file.path(results_dir, "p_maxIPP_XB.test.csv"), row.names = F)
}
