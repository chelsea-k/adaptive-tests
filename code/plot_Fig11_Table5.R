# Script for plotting Figure 11 and Table 5 from Krantsevich et al., 
# "BAYESIAN DECISION THEORY FOR TREE-BASED ADAPTIVE SCREENING TESTS 
#     WITH AN APPLICATION TO YOUTH DELINQUENCY"

source('code/util_functions.R')
library(caret)
library(tidyr)
library(dplyr)
library(cdata)
library(xtable)

# paraemeters
maxIPP.list=c(2:15)
w.list = c(0.65)

# Set up folders and load data
data.test = read.csv("data/item_response_data/item.response.data.test.csv")
sub.rows = which(data.test$Age>=15)

data.dir = "data/out_of_sample/subpopulation/synthetic_data/"
results.dir = "data/out_of_sample/subpopulation/results/"

synth.df.sub = read.csv(paste0(data.dir, "synth.XB.csv"))
p.sub.synth = read.csv(paste0(results.dir, "p.maxIPP.XB.synth.csv"))
p.sub.test = read.csv(paste0(results.dir, "p.maxIPP.XB.test.csv"))

synth.df.all = read.csv("data/out_of_sample/all/synthetic_data/synth.XB.csv")
p.all.synth = read.csv("data/out_of_sample/all/results/p.maxIPP.XB.synth.csv")
p.all.test = read.csv("data/out_of_sample/all/results/p.maxIPP.XB.test.csv")

results = data.frame(matrix(NA, nrow=0, ncol=5))
colnames(results) = c("maxIPP", "w", "Population", "Sensitivity", "Specificity")

for (i in seq_along(maxIPP.list)){
  maxIPP=maxIPP.list[[i]]
  m.col = which(colnames(p.sub.synth)==paste0("m.", maxIPP))
  p.sub.synth.temp = p.sub.synth[,m.col]
  p.sub.test.temp = p.sub.test[,m.col]
  p.all.synth.temp = p.all.synth[,m.col]
  p.all.test.temp = p.all.test[,m.col]
  
  cat(paste0("---------------- maxIPP = ", maxIPP, " ----------------\n"))
  for (j in seq_along(w.list)){
    w = w.list[[j]]
    cutoff.sub = get_cutoff(as.factor(synth.df.sub$y), p.sub.synth.temp, w)
    y.sub = as.integer(p.sub.test.temp >= cutoff.sub)
    cm.sub = confusionMatrix(as.factor(y.sub[sub.rows]), 
                                as.factor(data.test$y[sub.rows]), positive="1")
    results = rbind(results, data.frame(Num.Items=maxIPP, w=w, Population="Age >= 15",
                                        Sensitivity=cm.sub$byClass[['Sensitivity']], 
                                        Specificity=cm.sub$byClass[['Specificity']]))
    
    cutoff.all = get_cutoff(as.factor(synth.df.all$y), p.all.synth.temp, w)
    y.all = as.integer(p.all.test.temp >= cutoff.all)
    cm.all = confusionMatrix(as.factor(y.all[sub.rows]), 
                                  as.factor(data.test$y[sub.rows]), positive="1")
    results = rbind(results, data.frame(Num.Items=maxIPP, w=w, Population="All Youth",
                                        Sensitivity=cm.all$byClass[['Sensitivity']], 
                                        Specificity=cm.all$byClass[['Specificity']]))
    
  }
}
results$Utility = (results$w)*(results$Sensitivity)+(1-results$w)*(results$Specificity)
write.csv(results, paste0(results.dir, "sens.spec.results.subpopulation.test.csv"), row.names=FALSE)

# Post-process results and create plot for Figure 11
results$Num.Items = as.factor(results$Num.Items)
cbPalette <- c("#000000","#7a7a7a")

results_long = gather(results, Quantity, Value, Sensitivity:Utility, factor_key=TRUE)
results_long_all = results_long[which(results_long$Population=="All Youth"),] 
results_long_age = results_long[which(results_long$Population=="Age >= 15"),]
results_long_diff = results_long_all
results_long_diff$Value = results_long_age$Value - results_long_all$Value
colnames(results_long_diff)[[5]] = "Difference"
results_long_diff$Higher=ifelse(results_long_diff$Difference > 0,"Age >= 15","All Youth")

ggplot(results_long_diff, aes(x=Num.Items, y=Difference, fill=Higher)) + 
  geom_bar(stat = "identity") + 
  facet_grid(rows=vars(Quantity)) +
  xlab("maxIPP") + scale_fill_manual(values=cbPalette) +
  theme(legend.position = "none")

# Reformat results to be as in Table 5 
cT <- dplyr::tribble(
  ~Population, ~Sensitivity, ~Specificity, ~Utility,
  "Age >= 15",    "Sub_Sens", "Sub_Spec", "Sub_Util",
  "All Youth",    "All_Sens", "All_Spec", "All_Util"
)
results_wide = blocks_to_rowrecs(results, cT, keyColumns = "Num.Items")
results_wide = results_wide[,c(1,2,5,3,6,4,7)]

# Export table to latex format to create Table 5 (needs some reformatting of columns)
dir.create("tables", recursive = TRUE) # in case doesn't exist
print(xtable(results_wide, type = "latex", digits=c(0,0,3,3,3,3,3,3)), 
      include.rownames=FALSE, file="tables/table5.tex")
