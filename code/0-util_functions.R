# Utility functions for plotting figures from Krantsevich et al., 
# "BAYESIAN DECISION THEORY FOR TREE-BASED ADAPTIVE SCREENING TESTS 
#     WITH AN APPLICATION TO YOUTH DELINQUENCY"

library(pROC) 

get_metrics <- function(pred.prob, resp.y, return.cutoff=FALSE) {
  roc.fit = roc(response=as.factor(resp.y), predictor=pred.prob, 
                levels=c('0', '1'), direction="<")
  if(return.cutoff){
    metrics = data.frame(Specificity=roc.fit$specificities, 
                         Sensitivity=roc.fit$sensitivities,
                         Cutoff=roc.fit$thresholds)
  } else {
    metrics = data.frame(Specificity=roc.fit$specificities, 
                         Sensitivity=roc.fit$sensitivities)
  }
  return(metrics)
}

get_cutoff = function(response, predictor, w){
  roc.fit = roc(as.factor(response), predictor, levels=c("0", "1"), 
                direction="<", ret="all_coords")
  wtd.avg = w*roc.fit$sensitivities+(1-w)*roc.fit$specificities
  cutoff = roc.fit$thresholds[which.max(wtd.avg)]
  return(cutoff)
}

get_utility = function(y.true, prob=NA, cut=NA, y.pred=NA, type="prob", w) {
  if(type=="prob"){
    sens = sum(prob >= cut & y.true == 1)/sum(y.true == 1)
    spec = sum(prob < cut & y.true == 0)/sum(y.true == 0)
    util = w*sens + (1-w)*spec
  } else if (type=="class") {
    sens = sum(y.pred == 1 & y.true == 1)/sum(y.true == 1)
    spec = sum(y.pred == 0 & y.true == 0)/sum(y.true == 0)
    util = w*sens + (1-w)*spec
  } else {
    stop("'type' must be one of 'prob', 'class' ")
  }
  return(list(sens=sens, spec=spec, util=util))
} 

get_class = function(p.df, cutoff, num.items, cols){
  y.class = apply(matrix(seq_along(num.items), nrow=1), MARGIN=2, 
                  function(i) {1*(p.df[,cols[i]] >= cutoff[i])})
  return(y.class)
}


