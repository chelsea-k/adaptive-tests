library(stringr)
output_files <- list.files("output", recursive = TRUE)
for (i in seq_along(output_files)) {
  n_files <- length(output_files)
  if (i %% 10 ==0 ) {
    cat(str_glue("Processing file {i} out of {n_files}\n\n"))
  }
  my_file <- output_files[[i]]
  path_to_file <- file.path("output", my_file)
  if (str_detect(my_file, "fit.CART.")) {
	cat("trying_to_load..\n")	  
    load(path_to_file)
    if (str_detect(my_file, "opt")) {
      opt_tree$y <- NULL
      opt_tree$where <- NULL
      save(opt_tree, file = path_to_file)
    } else {
      fit_rpart$y <- NULL
      fit_rpart$where <- NULL
      save(fit_rpart, file = path_to_file)
    }
  }
}

