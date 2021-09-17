library(dplyr)

# Define where the item response (IR) data is stored
IR_data_dir <- "~/Desktop/ASU/Research/CAT_project/original_data"
IR_data_path <- file.path(IR_data_dir, "Honduras_complete_data_w_demo.csv")
store_data_dir <- "~/Desktop/ASU/Research/CAT_project/preprocessed_original_data/"

IR_data <- read.csv(IR_data_path)

# Define column names for different outcome types
outcome_col <- "vio_beh"
item_cols <- setdiff(colnames(IR_data), c("vio_beh","pro_cri","gan_inv",             
                                                   "dru_use", "dru_sel", "car_wea","tru", 
                                                   "Municipality", "Zone", "Gender",
                                                   "Age", "Birthplace", "Religion", 
                                                   "Importance.Religion", "Interview.date"))
demo_cols <- "Age"
y <- "y"


# Select only columns we need for the analysis
IR_data_all <- IR_data %>% 
  select(all_of(c(outcome_col, demo_cols, item_cols, "Interview.date"))) 
colnames(IR_data_all)[which(colnames(IR_data_all)==outcome_col)] <- "y"

# Create train / test splits based on interview date
IR_data_all <- IR_data_all %>%
  mutate(Interview.date = anydate(Interview.date)) %>%
  mutate(in_first_data_collection = Interview.date < anydate("2017-11-30"))

IR_data_train <- IR_data_all %>% 
  filter(in_first_data_collection) %>% select(-c(in_first_data_collection, Interview.date))

IR_data_test <- IR_data_all %>% 
  filter(!in_first_data_collection) %>% select(-c(in_first_data_collection, Interview.date))

IR_data_all <- IR_data_all %>%
  select(-c(in_first_data_collection, Interview.date))#
  

write.csv(IR_data_all, 
          file.path(store_data_dir, "IMC_data_all_preprocessed.csv"), 
          row.names = FALSE)

write.csv(IR_data_train, 
          file.path(store_data_dir, "IMC_data_train_preprocessed.csv"), 
          row.names = FALSE)

write.csv(IR_data_test, 
          file.path(store_data_dir, "IMC_data_test_preprocessed.csv"), 
          row.names = FALSE)




  



