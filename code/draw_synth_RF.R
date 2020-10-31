# Script for generating synthetic data using local perturbations and 
# a Random Forest model

library(randomForest)

out.of.sample = TRUE
subpopulation = TRUE

########################## Hyperparamters ##############################

n.samp = 10000
n.prune = 1000

########################## Data Preparation ##############################

# read in data
if(out.of.sample){
  data.train = read.csv("data/item_response_data/item.response.data.train.csv")
  data.test = read.csv("data/item_response_data/item.response.data.test.csv")
} else {
  data.train = read.csv("data/item_response_data/item.response.data.all.csv")
  data.test = data.train
}

if (subpopulation){
  data.train = data.train[which(data.train$Age >= 15),]
  data.test = data.test[which(data.test$Age >= 15),]
}

item.cols = which(!(colnames(data.train) %in% c("y", "Age"))) 
item.names = colnames(data.train)[item.cols]

# directories for model and data storage
folder = ifelse(out.of.sample, "data/out_of_sample/", "data/in_sample/")
folder = paste0(folder, ifelse(subpopulation, "subpopulation/", "all/"))
data.dir = paste0(folder, "synthetic_data/")
model.dir = paste0(folder, "model_fits/")
dir.create(data.dir, recursive = TRUE)
dir.create(model.dir, recursive = TRUE)

########################## Model Fitting, Data Generation ##############################

# create artificial item responses via local perturbations 
n.copies = floor((n.samp + n.prune)/nrow(data.train)) + 1
synth.data = as.matrix(data.train) %x% matrix(1,n.copies,1)  # duplicate original data 
synth.data = synth.data + sample(c(-1,0,1), length(synth.data), replace=TRUE)    # add 0,1,-1 with equal prob
for(i in 1:ncol(synth.data)){    # truncate if outside observed range
  synth.data[synth.data[,i]<min(data.train[,i]),i] = min(data.train[,i])
  synth.data[synth.data[,i]>max(data.train[,i]),i] = max(data.train[,i])
}
synth.data = data.frame(synth.data[,item.cols])
colnames(synth.data) = item.names   # turn into a dataframe with proper column names

# Fit Random Forest model for predicting "at-risk" probability
fit.RF = randomForest(y ~., data = cbind(y=as.factor(data.train$y), data.train[,item.cols]))
save(fit.RF, file=paste0(model.dir, "fit.randomForest"), ascii = TRUE)

# Generate synthetic "at-risk" probs by predicting from fitted model
synth.data$phat = predict(fit.RF, newdata=synth.data, type="prob")[,"1"]
synth.data$y = predict(fit.RF, newdata=synth.data, type="response")
synth.data = synth.data[sample(nrow(synth.data)),] 

# Save generated data
write.csv(synth.data[1:n.samp,], file=paste0(data.dir, "synth.RF.csv"), row.names = F)
write.csv(synth.data[(n.samp+1):(n.samp+n.prune),], file=paste0(data.dir, "prune.RF.csv"), row.names = F)
