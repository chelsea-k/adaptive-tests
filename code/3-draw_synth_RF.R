# Script for generating synthetic data using local perturbations and 
# a Random Forest model

library(randomForest)

out_of_sample = TRUE
subpopulation = TRUE

########################## Hyperparamters ##############################

n_samp = 1000000
n_prune = 100000

########################## Data Preparation ##############################

# read in data
# if(out_of_sample){
#   data_train = read.csv("data/item_response_data/item.response.data_train.csv")
#   data_test = read.csv("data/item_response_data/item.response.data_test.csv")
# } else {
#   data_train = read.csv("data/item_response_data/item.response.data.all.csv")
#   data_test = data_train
# }
original_data_dir <- "~/Desktop/ASU/Research/CAT_project/preprocessed_original_data"
if(out_of_sample){
  data_train = read.csv(file.path(original_data_dir, "IMC_data_train_preprocessed.csv"))
  data_test = read.csv(file.path(original_data_dir, "IMC_data_test_preprocessed.csv"))
} else {
  data_train = read.csv(file.path(original_data_dir, "IMC_data_all_preprocessed.csv"))
  data_test = data_train
}

if (subpopulation){
  data_train = data_train[which(data_train$Age >= 15),]
  data_test = data_test[which(data_test$Age >= 15),]
}

item_cols = which(!(colnames(data_train) %in% c("y", "Age"))) 
item_names = colnames(data_train)[item_cols]

# directories for model and data storage
# folder = ifelse(out_of_sample, "data/out_of_sample/", "data/in_sample/")
# folder = paste0(folder, ifelse(subpopulation, "subpopulation/", "all/"))
# data_dir = paste0(folder, "synthetic_data/")
# model_dir = paste0(folder, "model_fits/")
# dir.create(data_dir, recursive = TRUE)
# dir.create(model_dir, recursive = TRUE)
output_dir <- "~/Desktop/ASU/Research/CAT_project/output"
folder <- ifelse(out_of_sample, 
                 file.path(output_dir, "out_of_sample"), 
                 file.path(output_dir, "in_sample"))
folder <- file.path(folder, ifelse(subpopulation, "subpopulation", "all"))
data_dir <- file.path(folder, "synthetic_data")
model_dir <- file.path(folder, "model_fits")
dir.create(data_dir, recursive = TRUE)
dir.create(model_dir, recursive = TRUE)

########################## Model Fitting, Data Generation ##############################

# create artificial item responses via local perturbations
cat("Generating synthetic item response data via local perturbations...\n") 
n_copies = floor((2*n_samp + n_prune)/nrow(data_train)) + 1
synth_data = as.matrix(data_train) %x% matrix(1,n_copies,1)  # duplicate original data 
synth_data = synth_data + sample(c(-1,0,1), length(synth_data), replace=TRUE)    # add 0,1,-1 with equal prob
cat("Truncating synthetic data to be in observed range...\n")
for(i in 1:ncol(synth_data)){    # truncate if outside observed range
  synth_data[synth_data[,i]<min(data_train[,i]),i] = min(data_train[,i])
  synth_data[synth_data[,i]>max(data_train[,i]),i] = max(data_train[,i])
}
synth_data = data.frame(synth_data[,item_cols])
colnames(synth_data) = item_names   # turn into a dataframe with proper column names

# Fit Random Forest model for predicting "at-risk" probability
cat("Fitting Random Forest object on IMC data...\n")
fit_RF = randomForest(y ~., data = cbind(y=as.factor(data_train$y), data_train[,item_cols]))
save(fit_RF, file=file.path(model_dir, "fit_randomForest"), ascii = TRUE)

# Generate synthetic "at-risk" probs by predicting from fitted model
cat("Predicting probabilities for synthetic item responses...\n")
synth_data$phat = predict(fit_RF, newdata=synth_data, type="prob")[,"TRUE"]
cat("Predicting classes for synthetic item responses...\n")
synth_data$y = as.integer(as.logical(predict(fit_RF, newdata=synth_data, type="response")))
synth_data = synth_data[sample(nrow(synth_data)),] 

# Save generated data
cat("Saving data...\n")
synth_treefitting <- synth_data[1:n_samp,]
synth_uncertainty <- synth_data[(n_samp+1):(2*n_samp),]
prune_data <- synth_data[(2*n_samp+1):(2*n_samp+n_prune),]

#write.csv(synth_treefitting, file=file.path(data_dir, "synth_treefitting_RF.csv"), row.names = F)
#write.csv(synth_uncertainty, file=file.path(data_dir, "synth_uncertainty_RF.csv"), row.names = F)
#write.csv(prune_data, file=file.path(data_dir, "prune_RF.csv"), row.names = F)
