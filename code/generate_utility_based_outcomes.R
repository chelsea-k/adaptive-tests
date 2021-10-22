# set parameters
out_of_sample <- FALSE
subpopulation <- FALSE 

# read in original data
IMC_data_all <- read.csv("preprocessed_original_data/IMC_data_all_preprocessed.csv")

# read in generated data
synth_data_root <- "output"
folder = ifelse(out_of_sample, 
                file.path(synth_data_root, "out_of_sample"), 
                file.path(synth_data_root, "in_sample"))
folder =  file.path(folder, ifelse(subpopulation, "subpopulation", "all"))
synth_data_folder = file.path(folder, "synthetic_data")
synth_treefitting_df = read.csv(file.path(synth_data_folder, "synth_treefitting_XB.csv"))
synth_uncertainty_df = read.csv(file.path(synth_data_folder, "synth_uncertainty_XB.csv"))

# compute threshold to determine gamma^*_k
pr_1 <- mean(synth_uncertainty_df$phat.mean)
pr_0 <- 1 - pr_1
w <- 0.6
U0 <- (1 - w) / pr_0
U1 <- w / pr_1
cutoff <- U0 / (U0 + U1)

# compute utility function based outcomes
gamma_star <- as.integer(synth_treefitting_df$phat.mean >= cutoff)
synth_treefitting_df$y.mean <- gamma_star

gamma_star_draw <- as.integer(synth_treefitting_df$phat.draw >= cutoff)
synth_treefitting_df$y.draw <- gamma_star_draw


# store data
synth_data_store_name <- paste0("synth_treefitting_XB_util_based_outcomes_w_", w, ".csv")
write.csv(synth_treefitting_df, 
          file.path(synth_data_folder, synth_data_store_name),
          row.names = FALSE)




