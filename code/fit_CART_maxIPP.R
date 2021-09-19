# Script to fit regression trees using maxIPP growing method and 
# pruning based on out-of-sample RMSE

library(rpartMaxVPP)

out_of_sample = FALSE
subpopulation = FALSE

########################## Hyperparamters ##############################
CART_params = list(minbucket=100,
                   maxIPP_vals=c(2:15),
                   cp=1e-7,
                   patience = 10)
n_maxIPP = length(CART_params$maxIPP_vals)

########################## Data Preparation ##############################

# read in item response data
# if(out_of_sample){
#   data_train = read.csv("data/item_response_data/item.response.data_train.csv")
#   data_test = read.csv("data/item_response_data/item.response.data_test.csv")
# } else {
#   data_train = read.csv("data/item_response_data/item.response.data.all.csv")
#   data_test = data_train
# }
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
# folder = ifelse(out_of_sample, "data/out_of_sample/", "data/in_sample/")
# folder =  paste0(folder, ifelse(subpopulation, "subpopulation/", "all/"))
# synth_RF = read.csv(paste0(folder, "synthetic_data/synth_RF.csv"))
# prune_RF = read.csv(paste0(folder, "synthetic_data/prune_RF.csv"))
# synth_XB = read.csv(paste0(folder, "synthetic_data/synth_XB.csv"))
# prune_XB = read.csv(paste0(folder, "synthetic_data/prune_XB.csv"))
synth_data_root <- "output"
folder = ifelse(out_of_sample, 
                file.path(synth_data_root, "out_of_sample"), 
                file.path(synth_data_root, "in_sample"))
folder =  file.path(folder, ifelse(subpopulation, "subpopulation", "all"))
synth_data_folder = file.path(folder, "synthetic_data")
synth_RF = read.csv(file.path(synth_data_folder, "synth_treefitting_RF.csv"))
prune_RF = read.csv(file.path(synth_data_folder, "prune_RF.csv"))
synth_XB = read.csv(file.path(synth_data_folder, "synth_treefitting_XB.csv"))
synth_XB_UQ = read.csv(file.path(synth_data_folder, "synth_uncertainty_XB.csv"))
prune_XB = read.csv(file.path(synth_data_folder, "prune_XB.csv"))
n_synth = nrow(synth_XB)

# get column info for demo and item covariates
item_names = setdiff(colnames(data_train), c("y", "Age"))
RF_item_cols = which(colnames(synth_RF) %in% item_names)
XB_item_cols = which(colnames(synth_XB) %in% item_names)

# set up storage folders
model_dir = file.path(folder, "model_fits_CART/")
results_dir = file.path(folder, "results/")
dir.create(model_dir, recursive = TRUE)
dir.create(results_dir, recursive = TRUE)

########################## Model Fitting ##############################

# Create df for saving RPART predictions 
p_maxIPP_RF_synth = p_maxIPP_XB_synth = data.frame(matrix(NA, nrow=nrow(synth_RF), ncol=n_maxIPP))
colnames(p_maxIPP_RF_synth) = colnames(p_maxIPP_XB_synth) = paste0("m.", CART_params$maxIPP_vals)

p_maxIPP_RF_test = p_maxIPP_XB_test = data.frame(matrix(NA, nrow=nrow(data_test), ncol=n_maxIPP))
colnames(p_maxIPP_RF_test) = colnames(p_maxIPP_XB_test) = paste0("m.", CART_params$maxIPP_vals)

# Define function for fitting CART via maxIPP and pruning
fit_CART_maxIPP = function(fitting_data, prune_data, predict_data, test_data, 
                           item_cols, savefile, p_df_synth, p_df_test, oos, params) {
 
maxIPP_vals = params$maxIPP_vals   
cp = params$cp
minbucket = params$minbucket
patience = params$patience

for (i in seq_along(maxIPP_vals)){
  maxIPP=maxIPP_vals[[i]]
  cat(paste0("---------------- maxIPP = ", maxIPP, " ----------------\n"))
  
  # check if tree already exists
  tree_savefile <- paste0(savefile, ".regression_maxIPP.", maxIPP)
  opt_tree_savefile <- paste0(savefile, ".regression_maxIPP.", maxIPP, ".opt")
  if(!file.exists(opt_tree_savefile)) {
    if(file.exists(tree_savefile)) {
      cat("Loading big regression tree \n")
      load(tree_savefile)
    } else {
        # fit the large CART tree
        fit_rpart = rpartMaxVPP::rpart(p ~., data = cbind(p=fitting_data$phat, fitting_data[,item_cols]),
                       method = "anova", cp=cp, maxdepth=25, minbucket = minbucket, maxvpp=maxIPP,
                       xval=0, maxcompete=0, maxsurrogate=0, usesurrogate=0, model=TRUE)
        save(fit_rpart, file=tree_savefile, ascii=TRUE)
    }
    
    # check if optimum tree exists 
    # find optimal pruning point based on threshold/patience
    cplist = fit_rpart$cptable
    threshold =  ifelse(maxIPP<5, 1e-4, 1e-5)
    opt_cp = cplist[1,"CP"]
    root_tree = prune(fit_rpart, cp=opt_cp)
    p_CART = predict(root_tree, newdata = prune_data[,item_cols], type='vector')
    rmse_prev = sqrt(mean((p_CART - prune_data$phat)^2))
    wait=0
    if(nrow(cplist) == 1) {
      opt_tree = fit_rpart
    } else {
        for (j in 2:nrow(cplist)) {
          if (j%%50==0) {cat("Predicting for cp: ", j, "out of", nrow(cplist),"\n")}
          this_cp = cplist[j,"CP"]
          pruned_tree = prune(fit_rpart, cp=this_cp)
          p_CART = predict(pruned_tree, newdata = prune_data[,item_cols], type='vector')
          rmse = sqrt(mean((p_CART - prune_data$phat)^2))
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
    save(opt_tree, file=opt_tree_savefile, ascii=TRUE)
  } else {
    cat("Loading optimally pruned regression tree \n")
    load(opt_tree_savefile)
  }
  
  # Store CART predictions from pruned tree on synthetic XBART data
  opt_p_CART = predict(opt_tree, predict_data)
  m_col = which(colnames(p_df_synth) == paste0("m.", maxIPP))
  p_df_synth[,m_col] = opt_p_CART
  
  if(oos){
    # Store CART predictions from pruned tree on test data
    opt_p_test = predict(opt_tree, test_data)
    m_col = which(colnames(p_df_test) == paste0("m.", maxIPP))
    p_df_test[,m_col] = opt_p_test
  }
  
}
if (out_of_sample){
  results = list(p_synth = p_df_synth, p_test = p_df_test)
} else {
  results = list(p_synth = p_df_synth)
}
return(results)
}

# Run function to get maxIPP results and store results
maxIPP_RF_results = fit_CART_maxIPP(fitting_data = synth_RF, prune_data = prune_RF, 
                             predict_data = synth_XB_UQ, test_data = data_test, 
                             item_cols = RF_item_cols, 
                             savefile = file.path(model_dir, "fit.CART.RF"),
                             p_df_synth = p_maxIPP_RF_synth, p_df_test = p_maxIPP_RF_test, 
                             oos = out_of_sample, params = CART_params)
write.csv(maxIPP_RF_results$p_synth, file.path(results_dir, "p_maxIPP_RF_synth.csv"), row.names = F)
if(out_of_sample) {
  write.csv(maxIPP_RF_results$p_test, file.path(results_dir, "p_maxIPP_RF_test.csv"), row.names = F)
}

maxIPP_XB_results = fit_CART_maxIPP(fitting_data = synth_XB, prune_data = prune_XB, 
                             predict_data = synth_XB_UQ, test_data = data_test, 
                             item_cols = XB_item_cols, 
                             savefile = file.path(model_dir, "fit.CART.XB"),
                             p_df_synth = p_maxIPP_XB_synth, p_df_test = p_maxIPP_XB_test, 
                             oos = out_of_sample, params = CART_params)

write.csv(maxIPP_XB_results$p_synth, file.path(results_dir, "p_maxIPP_XB_synth.csv"), row.names = F)
if (out_of_sample){
  write.csv(maxIPP_XB_results$p_test, file.path(results_dir, "p_maxIPP_XB_test.csv"), row.names = F)
}