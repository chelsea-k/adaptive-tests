library(stringr)
root <- "~/Desktop/output/"
synth_files <- list.files(root, recursive = TRUE)

for (i in seq_along(synth_files)) {
  file_name <- synth_files[[i]]
  if(grepl("synth_treefitting|synth_uncertainty|prune", file_name)){
    cat(paste(file_name, "\n", sep=""))
    df <- read.csv(file.path(root, file_name))
    df$phat <- df$phat/2
    write.csv(df, file.path(root, file_name), row.names = FALSE)
  }
}

