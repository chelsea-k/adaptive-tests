# Script for recreating Table 3/Table in Appendix C from Krantsevich et al., 
# "BAYESIAN DECISION THEORY FOR TREE-BASED ADAPTIVE SCREENING TESTS 
#     WITH AN APPLICATION TO YOUTH DELINQUENCY"

source('code/util_functions.R')
library(xtable)

# Set parameters
which.table = "table3" # "table3" or "app.C"
if(which.table == "table3"){
  maxIPP.vals = c(4,8,12)
} else if (which.table == "app.C") {
  maxIPP.vals = c(2:15)
} else {
  stop(" 'which.table' must be one of 'table3' or 'app.C' ")
}
n.maxIPP = length(maxIPP.vals)
w.list = c(0.5,0.65,0.8)

# Set up folders and load data
data.dir = "data/out_of_sample/all/synthetic_data/"
results.dir = "data/out_of_sample/all/results/"

data.test = read.csv("data/item_response_data/item.response.data.test.csv")
y.test = data.test$y

synth.df.XB = read.csv(paste0(data.dir, "synth.XB.csv"))
y.synth = synth.df.XB$y

p.maxIPP.XB.synth = read.csv(paste0(results.dir, "p.maxIPP.XB.synth.csv"))
p.maxIPP.RF.synth = read.csv(paste0(results.dir, "p.maxIPP.RF.synth.csv"))
p.maxdepth.XB.synth = read.csv(paste0(results.dir, "p.maxDepth.XB.synth.csv"))
p.maxdepth.RF.synth = read.csv(paste0(results.dir, "p.maxDepth.RF.synth.csv"))
y.class.real.synth = read.csv(paste0(results.dir, "y.real.synth.csv"))
y.class.RF.synth = read.csv(paste0(results.dir, "y.RF.synth.csv"))
y.class.XB.synth = read.csv(paste0(results.dir, "y.XB.synth.csv"))

p.maxIPP.XB.test = read.csv(paste0(results.dir, "p.maxIPP.XB.test.csv"))
p.maxIPP.RF.test = read.csv(paste0(results.dir, "p.maxIPP.RF.test.csv"))
p.maxdepth.XB.test = read.csv(paste0(results.dir, "p.maxDepth.XB.test.csv"))
p.maxdepth.RF.test = read.csv(paste0(results.dir, "p.maxDepth.RF.test.csv"))
y.class.real.test = read.csv(paste0(results.dir, "y.real.test.csv"))
y.class.RF.test = read.csv(paste0(results.dir, "y.RF.test.csv"))
y.class.XB.test = read.csv(paste0(results.dir, "y.XB.test.csv"))

n.methods = 7

# Compute sens/spec/util results for all values of w and maxIPP
sens.spec.results <- data.frame(matrix(NA, nrow=0, ncol=7))
colnames(sens.spec.results) <- c("Num.Items", "Tree.Type", "Stop.Criterion", "w", 
                               "Fitting.Data", "Sensitivity", "Specificity")

for (t in seq_along(w.list)){
  w = w.list[[t]]
  
  # Compute optimal cutoffs and predicted classes 
  maxIPP.XB.cut = maxIPP.RF.cut = maxdepth.XB.cut = maxdepth.RF.cut = rep(NA, n.maxIPP)
  m.colnames = paste0("m.", maxIPP.vals)
  m.cols = which(colnames(p.maxIPP.XB.synth) %in% m.colnames)
  for (i in seq_along(maxIPP.vals)){
    maxIPP.XB.cut[[i]] = get_cutoff(as.factor(y.synth), p.maxIPP.XB.synth[,m.cols[i]], w)
    maxIPP.RF.cut[[i]] = get_cutoff(as.factor(y.synth), p.maxIPP.RF.synth[,m.cols[i]], w)
    maxdepth.XB.cut[[i]] = get_cutoff(as.factor(y.synth), p.maxdepth.XB.synth[,m.cols[i]], w)
    maxdepth.RF.cut[[i]] = get_cutoff(as.factor(y.synth), p.maxdepth.RF.synth[,m.cols[i]], w)
  }
  y.maxIPP.XB.test = get_class(p.maxIPP.XB.test, maxIPP.XB.cut, maxIPP.vals, m.cols)
  y.maxIPP.RF.test = get_class(p.maxIPP.RF.test, maxIPP.RF.cut, maxIPP.vals, m.cols)
  y.maxDepth.XB.test = get_class(p.maxdepth.XB.test, maxdepth.XB.cut, maxIPP.vals, m.cols)
  y.maxDepth.RF.test = get_class(p.maxdepth.RF.test, maxdepth.RF.cut, maxIPP.vals, m.cols)
  # Extract columns of maxIPP vals being used
  y.class.real.test = y.class.real.test[,which(colnames(y.class.real.test) %in% m.colnames)]
  y.class.RF.test = y.class.RF.test[,which(colnames(y.class.RF.test) %in% m.colnames)]
  y.class.XB.test = y.class.XB.test[,which(colnames(y.class.XB.test) %in% m.colnames)]
  
  # Compute sensitivity and specificity for all methods
  for (i in seq_along(maxIPP.vals)){
    maxIPP.XB.metrics = get_utility(y.test,y.pred=y.maxIPP.XB.test[,i],type="class",w=w)
    maxdepth.XB.metrics = get_utility(y.test,y.pred=y.maxDepth.XB.test[,i],type="class",w=w)
    maxIPP.RF.metrics = get_utility(y.test,y.pred=y.maxIPP.RF.test[,i],type="class",w=w)
    maxdepth.RF.metrics = get_utility(y.test,y.pred=y.maxDepth.RF.test[,i],type="class",w=w)
    class.XB.metrics  = get_utility(y.test,y.pred=y.class.XB.test[,i],type="class",w=w)
    class.RF.metrics  = get_utility(y.test,y.pred=y.class.RF.test[,i],type="class",w=w)
    class.real.metrics  = get_utility(y.test,y.pred=y.class.real.test[,i],type="class",w=w)
    
    temp.results = data.frame(Num.Items = rep(maxIPP.vals[[i]], n.methods),
                              Tree.Type = c("Regression", "Regression", "Regression", "Regression", 
                                            "Classification", "Classification", "Classification"),
                              Stop.Criterion = c("maxIPP", "maxDepth", "maxIPP", "maxDepth", 
                                                 "maxDepth", "maxDepth", "maxDepth"),
                              w = c(rep(w, 4), rep("NA", 3)),
                              Fitting.Data = c("BFA + XBART", "BFA + XBART", 
                                               "Perturb + RF", "Perturb + RF",
                                               "BFA + XBART", "Perturb + RF", "Real Training Data"),
                              Sensitivity = c(maxIPP.XB.metrics$sens, maxdepth.XB.metrics$sens, 
                                              maxIPP.RF.metrics$sens, maxdepth.RF.metrics$sens,
                                              class.XB.metrics$sens, class.RF.metrics$sens,
                                              class.real.metrics$sens),
                              Specificity = c(maxIPP.XB.metrics$spec, maxdepth.XB.metrics$spec, 
                                              maxIPP.RF.metrics$spec, maxdepth.RF.metrics$spec,
                                              class.XB.metrics$spec, class.RF.metrics$spec,
                                              class.real.metrics$spec)
                              )
    sens.spec.results = rbind(sens.spec.results, temp.results)
    }
}

# Post-processing and exporting Table 3 and Table from Appendix C to latex
sens.spec.results = sens.spec.results[!duplicated(sens.spec.results),]
sens.spec.results = sens.spec.results[order(sens.spec.results$Num.Items, sens.spec.results$w),]
write.csv(sens.spec.results, paste0(results.dir, "sens.spec.results.all.test.csv"), row.names = FALSE)

dir.create("tables", recursive = TRUE) # in case doesn't exist
if (which.table == "table3") {
  sens.spec.results.summary = sens.spec.results[which(!(sens.spec.results$Stop.Criterion == "maxDepth" & 
                                                          sens.spec.results$Tree.Type == "Regression")),]
  print(xtable(sens.spec.results.summary, type = "latex", digits=c(0,0,0,0,2,0,3,3)), 
        include.rownames=FALSE, file="tables/table3.tex")
}else if(which.table == "app.C") {
    print(xtable(sens.spec.results, type = "latex", digits=c(0,0,0,0,2,0,3,3)), 
          include.rownames=FALSE, file="tables/table_AppendixC.tex")
} else {
  stop(" 'which.table' must be one of 'table3' or 'app.C' ")
}
