# Script for plotting Figures 2 and 3 from Krantsevich et al., 
# "BAYESIAN DECISION THEORY FOR TREE-BASED ADAPTIVE SCREENING TESTS 
#     WITH AN APPLICATION TO YOUTH DELINQUENCY"

source("code/util_functions.R") # get_metrics()
library(ggplot2)

# Choose cutoffs C for which we'll plot corresponding (spec, sens) pairs (in Fig 2)
cutoffs = c(C1 = 0.10, C2 = 0.25, C3 = 0.50)

# Define function for getting sensitivity/specificity pairs


# Read in the data for ROC curve; choose a few post. idx for faster plotting
synth.df = read.csv(paste0("data/in_sample/all/synthetic_data/synth_treefitting_XB.csv"))
temp.df = synth.df[which(synth.df$post.idx %in% c(1:100)),] 

# Compute sensitivity/specificity/cutoffs for synth data
metrics = get_metrics(temp.df$phat, temp.df$y, return.cutoff = TRUE)

# Create dataframe of (spec, sens) pairs and corresponding cutoffs for Fig 2
cutoffs.idx = lapply(cutoffs, function(x) {which.min(abs(metrics$Cutoff-x))})
cutoffs.df = metrics[unlist(cutoffs.idx),]
cutoffs.df$Cutoff = factor(format(cutoffs, digits = 2))

# Create plot for Figure 2
ggplot(metrics, aes(x=Specificity, y=Sensitivity)) +
  geom_line(color="gray45") + xlim(c(1,0)) + 
  geom_point(data=cutoffs.df, size=4, aes(shape=Cutoff)) +
  scale_shape_manual(values=c(15,16,17))+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("ROC Curve Example") +
  theme_bw()

# Create dataframes for point, rectangle, and line segments in Figure 3
point.df = metrics[which.max(metrics$Specificity+metrics$Sensitivity),]
max.spec = point.df$Specificity
max.sens = point.df$Sensitivity
rect.df = data.frame(x1 = 0,  y1 = 0,
                     x2 = max.spec, y2 = max.sens)
seg.df = data.frame(x1 = c(max.spec, max.spec), y1 = c(0, max.sens), 
                    x2 = c(max.spec, 0), y2 = c(max.sens, max.sens), 
                    Quantity = c("Sensitivity", "Specificity"))

# Create plot for Figure 3
ggplot(metrics, aes(x=Specificity, y=Sensitivity)) +
  geom_line() + xlim(c(1,0)) + 
  geom_rect(data = rect.df, inherit.aes=FALSE, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2),
            color="black", fill="grey75") +
  geom_segment(data=seg.df, inherit.aes=FALSE,color="black", size=1.5,
               aes(x = x1, xend = x2, y = y1, yend = y2)) +
  geom_point(data=point.df, size=4) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Utility Function and ROC Curve") +
  theme_bw()




