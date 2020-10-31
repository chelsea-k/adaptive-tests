# Script to grow regression and classification trees using maximum depth criterion

library(rpartMaxVPP) # can also use rpart for this one

out.of.sample = FALSE
subpopulation = FALSE

########################## Hyperparamters ##############################
CART.params = list(maxDepth.vals=c(2:15))
n.maxDepth = length(CART.params$maxDepth.vals)

########################## Data Preparation ##############################

# read in item response data
if(out.of.sample){
  data.train = read.csv("data/item_response_data/item.response.data.train.csv")
  data.test = read.csv("data/item_response_data/item.response.data.test.csv")
} else {
  data.train = read.csv("data/item_response_data/item.response.data.all.csv")
  data.test = data.train
}

# read in generated data
folder = ifelse(out.of.sample, "data/out_of_sample/", "data/in_sample/")
folder =  paste0(folder, ifelse(subpopulation, "subpopulation/", "all/"))
synth.RF = read.csv(paste0(folder, "synthetic_data/synth.RF.csv"))
synth.RF$y[1] = 1
synth.XB = read.csv(paste0(folder, "synthetic_data/synth.XB.csv"))
n.synth = nrow(synth.RF)

# get column info for demo and item covariates
item.names = setdiff(colnames(data.train), c("y", "Age"))
RF.item.cols = which(colnames(synth.RF) %in% item.names)
XB.item.cols = which(colnames(synth.XB) %in% item.names)
real.item.cols = which(colnames(data.train) %in% item.names)

# set up storage folders
model.dir = paste0(folder, "model_fits/")
results.dir = paste0(folder, "results/")
dir.create(model.dir, recursive = TRUE)
dir.create(results.dir, recursive = TRUE)


########################## Model Fitting ##############################

# Create dfs for saving RPART predictions; 
# "p." dfs are for regression trees, "y." dfs are for classification trees
m.colnames = paste0("m.", CART.params$maxDepth.vals)

p.maxDepth.RF.synth = p.maxDepth.XB.synth = data.frame(matrix(NA, nrow=n.synth, ncol=n.maxDepth))
colnames(p.maxDepth.RF.synth) = colnames(p.maxDepth.XB.synth) = m.colnames

y.real.synth = y.RF.synth = y.XB.synth = data.frame(matrix(NA, nrow=n.synth, ncol=n.maxDepth))
colnames(y.real.synth) = colnames(y.RF.synth) = colnames(y.XB.synth) = m.colnames
  
p.maxDepth.RF.test = p.maxDepth.XB.test = data.frame(matrix(NA, nrow=nrow(data.test), ncol=n.maxDepth))
colnames(p.maxDepth.RF.test) = colnames(p.maxDepth.XB.test) = m.colnames

y.real.test = y.RF.test = y.XB.test = data.frame(matrix(NA, nrow=nrow(data.test), ncol=n.maxDepth))
colnames(y.real.test) = colnames(y.RF.test) = colnames(y.XB.test) = m.colnames

# Define function for fitting CART via maxdepth (both regression and classification)
fit_CART_maxDepth = function(fitting.data, predict.data, test.data, 
                             item.cols, savefile, pred.df.synth, pred.df.test, 
                             tree.type, oos, params) {

maxDepth.vals = params$maxDepth.vals

for (i in seq_along(maxDepth.vals)){
  maxdepth=maxDepth.vals[[i]]
  m.col = which(colnames(pred.df.synth) == paste0("m.", maxdepth))
  cat(paste0("---------------- maxdepth = ", maxdepth, " ----------------\n"))
  if (tree.type == "classification") {
    cat("Fitting classifcation tree \n")
    fit.rpart = rpart(y ~.,data = cbind(y=as.factor(fitting.data$y), fitting.data[,item.cols]), 
                      cp=0, method = "class", maxdepth=maxdepth, xval=0, maxcompete=0,
                      maxsurrogate=0, usesurrogate=0)
    save(fit.rpart, file=paste0(savefile, ".classification.maxDepth.", maxdepth), ascii=TRUE)
    
    cat("Predicting with classification tree \n")
    y.synth = predict(fit.rpart, predict.data, type="class")
    pred.df.synth[,m.col] = as.integer(as.character(y.synth))
    if(out.of.sample){
      y.test = predict(fit.rpart, test.data, type="class")
      pred.df.test[,m.col] = as.integer(as.character(y.test))
    }
  }
  
  if (tree.type == "regression") {
    cat("Fitting regression tree \n")
    fit.rpart = rpart(phat ~.,data = cbind(phat=fitting.data$phat, fitting.data[,item.cols]), 
                      cp=0, method = "anova", maxdepth=maxdepth, xval=0, maxcompete=0,
                      maxsurrogate=0, usesurrogate=0)
    save(fit.rpart, file=paste0(savefile, ".regression.maxDepth.", maxdepth), ascii=TRUE)
    
    cat("Predicting with regression tree \n")
    pred.df.synth[,m.col] = predict(fit.rpart, predict.data)
    if(out.of.sample){
      pred.df.test[,m.col] = predict(fit.rpart, test.data)
    }
  }
}
if (out.of.sample){
  results = list(pred.synth = pred.df.synth, pred.test = pred.df.test)
} else {
  results = list(pred.synth = pred.df.synth)
}
return(results)
}

# Run function to get maxDepth results and store results
maxDepth.RF.results = fit_CART_maxDepth(fitting.data = synth.RF, predict.data = synth.XB,
                                        test.data = data.test, item.cols = RF.item.cols, 
                                        savefile = paste0(model.dir, "fit.CART.RF"),
                                        pred.df.synth = p.maxDepth.RF.synth, 
                                        pred.df.test = p.maxDepth.RF.test, 
                                        tree.type = "regression", oos = out.of.sample, 
                                        params = CART.params)

maxDepth.XB.results = fit_CART_maxDepth(fitting.data = synth.XB, predict.data = synth.XB, 
                                        test.data = data.test, item.cols = XB.item.cols, 
                                        savefile = paste0(model.dir, "fit.CART.XB"),
                                        pred.df.synth = p.maxDepth.XB.synth, 
                                        pred.df.test = p.maxDepth.XB.test, 
                                        tree.type = "regression", oos = out.of.sample, 
                                        params = CART.params)

class.real.results = fit_CART_maxDepth(fitting.data = data.train, predict.data = synth.XB,
                                       test.data = data.test, item.cols = real.item.cols, 
                                       savefile = paste0(model.dir, "fit.CART.real"),
                                       pred.df.synth = y.real.synth, 
                                       pred.df.test = y.real.test, 
                                       tree.type = "classification", oos = out.of.sample, 
                                       params = CART.params)

class.RF.results = fit_CART_maxDepth(fitting.data = synth.RF, predict.data = synth.XB,
                                     test.data = data.test, item.cols = RF.item.cols, 
                                     savefile = paste0(model.dir, "fit.CART.RF"),
                                     pred.df.synth = y.RF.synth, 
                                     pred.df.test = y.RF.test, 
                                     tree.type = "classification", oos = out.of.sample, 
                                     params = CART.params)

class.XB.results = fit_CART_maxDepth(fitting.data = synth.XB, predict.data = synth.XB,
                                     test.data = data.test, item.cols = XB.item.cols, 
                                     savefile = paste0(model.dir, "fit.CART.XB"),
                                     pred.df.synth = y.XB.synth, 
                                     pred.df.test = y.XB.test, 
                                     tree.type = "classification", oos = out.of.sample, 
                                     params = CART.params)

write.csv(maxDepth.RF.results$pred.synth, paste0(results.dir, "p.maxDepth.RF.synth.csv"), row.names = F)
write.csv(maxDepth.XB.results$pred.synth, paste0(results.dir, "p.maxDepth.XB.synth.csv"), row.names = F)
write.csv(class.real.results$pred.synth, paste0(results.dir, "y.real.synth.csv"), row.names = F)
write.csv(class.RF.results$pred.synth, paste0(results.dir, "y.RF.synth.csv"), row.names = F)
write.csv(class.XB.results$pred.synth, paste0(results.dir, "y.XB.synth.csv"), row.names = F)

if (out.of.sample){
  write.csv(maxDepth.RF.results$pred.test, paste0(results.dir, "p.maxDepth.RF.test.csv"), row.names = F)
  write.csv(maxDepth.XB.results$pred.test, paste0(results.dir, "p.maxDepth.XB.test.csv"), row.names = F)
  write.csv(class.real.results$pred.test, paste0(results.dir, "y.real.test.csv"), row.names = F)
  write.csv(class.RF.results$pred.test, paste0(results.dir, "y.RF.test.csv"), row.names = F)
  write.csv(class.XB.results$pred.test, paste0(results.dir, "y.XB.test.csv"), row.names = F)
}

