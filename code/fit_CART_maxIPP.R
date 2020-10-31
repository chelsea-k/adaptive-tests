# Script to fit regression trees using maxIPP growing method and 
# pruning based on out-of-sample RMSE

library(rpartMaxVPP)

out.of.sample = TRUE
subpopulation = TRUE

########################## Hyperparamters ##############################
CART.params = list(minbucket=100,
                   maxIPP.vals=c(2:15),
                   cp=1e-7,
                   patience = 10)
n.maxIPP = length(CART.params$maxIPP.vals)

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
prune.RF = read.csv(paste0(folder, "synthetic_data/prune.RF.csv"))
synth.XB = read.csv(paste0(folder, "synthetic_data/synth.XB.csv"))
prune.XB = read.csv(paste0(folder, "synthetic_data/prune.XB.csv"))

# get column info for demo and item covariates
item.names = setdiff(colnames(data.train), c("y", "Age"))
RF.item.cols = which(colnames(synth.RF) %in% item.names)
XB.item.cols = which(colnames(synth.XB) %in% item.names)

# set up storage folders
model.dir = paste0(folder, "model_fits/")
results.dir = paste0(folder, "results/")
dir.create(model.dir, recursive = TRUE)
dir.create(results.dir, recursive = TRUE)

########################## Model Fitting ##############################

# Create df for saving RPART predictions 
p.maxIPP.RF.synth = p.maxIPP.XB.synth = data.frame(matrix(NA, nrow=nrow(synth.RF), ncol=n.maxIPP))
colnames(p.maxIPP.RF.synth) = colnames(p.maxIPP.XB.synth) = paste0("m.", CART.params$maxIPP.vals)

p.maxIPP.RF.test = p.maxIPP.XB.test = data.frame(matrix(NA, nrow=nrow(data.test), ncol=n.maxIPP))
colnames(p.maxIPP.RF.test) = colnames(p.maxIPP.XB.test) = paste0("m.", CART.params$maxIPP.vals)

# Define function for fitting CART via maxIPP and pruning
fit_CART_maxIPP = function(fitting.data, prune.data, predict.data, test.data, 
                           item.cols, savefile, p.df.synth, p.df.test, oos, params) {
 
maxIPP.vals = params$maxIPP.vals   
cp = params$cp
minbucket = params$minbucket
patience = params$patience

for (i in seq_along(maxIPP.vals)){
  maxIPP=maxIPP.vals[[i]]
  cat(paste0("---------------- maxIPP = ", maxIPP, " ----------------\n"))
  # fit the large CART tree
  fit.rpart = rpart(p ~., data = cbind(p=fitting.data$phat, fitting.data[,item.cols]),
                 method = "anova", cp=cp, maxdepth=25, minbucket = minbucket, maxvpp=maxIPP,
                 xval=0, maxcompete=0, maxsurrogate=0, usesurrogate=0, model=TRUE)
  save(fit.rpart, file=paste0(savefile, ".regression.maxIPP.", maxIPP), ascii=TRUE)
  
  # find optimal pruning point based on threshold/patience
  cplist = fit.rpart$cptable
  threshold =  ifelse(maxIPP<5, 1e-4, 1e-5)
  opt.cp = cplist[1,"CP"]
  root.tree = prune(fit.rpart, cp=opt.cp)
  p.CART = predict(root.tree, newdata = prune.data[,item.cols], type='vector')
  rmse.prev = sqrt(mean((p.CART - prune.data$phat)^2))
  wait=0
  for (j in 2:nrow(cplist)) {
    if (j%%50==0) {cat("Predicting for cp: ", j, "out of", nrow(cplist),"\n")}
    this.cp = cplist[j,"CP"]
    pruned.tree = prune(fit.rpart, cp=this.cp)
    p.CART = predict(pruned.tree, newdata = prune.data[,item.cols], type='vector')
    rmse = sqrt(mean((p.CART - prune.data$phat)^2))
    rmse.diff = rmse.prev - rmse
    if (rmse.diff > threshold){
        rmse.prev = rmse
        opt.cp = this.cp
        wait=0
    } else {
        wait = wait + 1
        rmse.prev = rmse
        if (wait >= patience){break}
    }  
  }
  
  # Prune to optimal cp and store optimally pruned tree
  opt.tree = prune(fit.rpart, cp=opt.cp)
  save(opt.tree, file=paste0(savefile, ".regression.maxIPP.", maxIPP, ".opt"), ascii=TRUE)
  
  # Store CART predictions from pruned tree on synthetic XBART data
  opt.p.CART = predict(opt.tree, predict.data)
  m.col = which(colnames(p.df.synth) == paste0("m.", maxIPP))
  p.df.synth[,m.col] = opt.p.CART
  
  if(oos){
    # Store CART predictions from pruned tree on test data
    opt.p.test = predict(opt.tree, test.data)
    m.col = which(colnames(p.df.test) == paste0("m.", maxIPP))
    p.df.test[,m.col] = opt.p.test
  }
  
}
if (out.of.sample){
  results = list(p.synth = p.df.synth, p.test = p.df.test)
} else {
  results = list(p.synth = p.df.synth)
}
return(results)
}

# Run function to get maxIPP results and store results
maxIPP.RF.results = fit_CART_maxIPP(fitting.data = synth.RF, prune.data = prune.RF, 
                             predict.data = synth.XB, test.data = data.test, 
                             item.cols = RF.item.cols, 
                             savefile = paste0(model.dir, "fit.CART.RF"),
                             p.df.synth = p.maxIPP.RF.synth, p.df.test = p.maxIPP.RF.test, 
                             oos = out.of.sample, params = CART.params)

maxIPP.XB.results = fit_CART_maxIPP(fitting.data = synth.XB, prune.data = prune.XB, 
                             predict.data = synth.XB, test.data = data.test, 
                             item.cols = XB.item.cols, 
                             savefile = paste0(model.dir, "fit.CART.XB"),
                             p.df.synth = p.maxIPP.XB.synth, p.df.test = p.maxIPP.XB.test, 
                             oos = out.of.sample, params = CART.params)

write.csv(maxIPP.RF.results$p.synth, paste0(results.dir, "p.maxIPP.RF.synth.csv"), row.names = F)
write.csv(maxIPP.XB.results$p.synth, paste0(results.dir, "p.maxIPP.XB.synth.csv"), row.names = F)

if (out.of.sample){
  write.csv(maxIPP.RF.results$p.test, paste0(results.dir, "p.maxIPP.RF.test.csv"), row.names = F)
  write.csv(maxIPP.XB.results$p.test, paste0(results.dir, "p.maxIPP.XB.test.csv"), row.names = F)
}