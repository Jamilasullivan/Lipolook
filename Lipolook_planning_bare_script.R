################################################################################
############################ VARIABLES TO FILL #################################
################################################################################

working_directory <- "C:/Users/jamsu/OneDrive - Cardiff University/University/Masters/Big Data Biology/Modules/Dissertation/R Programme/Dissertation/test_data"

raw_data_file_name <- "raw_data_1.csv"

lipids_tested_file_name <- "lipids_tested_1.csv"

################################################################################
############################ PACKAGES TO INSTALL ###############################
################################################################################

library(tidyr)
library(dplyr)
library(ggplot2)

################################################################################
######################### SETTING WORK DIRECTOTY ###############################
################################################################################

setwd(working_directory)

################################################################################
####################### CREATING DIRECTORY STRUCTURE ###########################
################################################################################

if(!dir.exists("outputs")){
  cat("'Outputs' directory created")
  dir.create("outputs")
} else {
  cat("'Outputs' directory already exists")
}

if(!dir.exists("outputs/total_lipids")){
  cat("'total_lipids' directory created")
  dir.create("outputs/total_lipids", recursive = TRUE)
} else {
  cat("'total_lipids' directory already exists")
}

if(!dir.exists("outputs/error_files")){
  cat("'error_files' directory created")
  dir.create("outputs/error_files", recursive = TRUE)
} else {
  cat("'error_files' directory already exists")
}

################################################################################
######################### INITIAL DATA TIDYING #################################
################################################################################

raw_data <- read.csv(raw_data_file_name, header = T)
names(raw_data)[c(1,2,3,4,5)] <- c("Filename", "CMaLL Sample Name", "Group", "code", "Cell Count")
raw_data <- raw_data[-1,]
raw_data$code <- NULL
raw_data$`Cell Count` <- NULL
str(raw_data)
raw_data[, 4:ncol(raw_data)] <- lapply(raw_data[, 4:ncol(raw_data)], function(x)
  as.numeric(as.character(x)))
sum(is.na(raw_data[, 4:ncol(raw_data)]))
str(raw_data)
sapply(raw_data, is.numeric)

################################################################################
######################### RAW DATA MANIPULATION ################################
################################################################################

raw_data_lipids <- raw_data
rownames(raw_data_lipids) <- raw_data_lipids[[2]]
raw_data_lipids <- raw_data_lipids[, -(1:2)]
duplicated_columns <- duplicated(as.list(raw_data_lipids))
any(duplicated_columns)
duplicated_column_names <- names(raw_data_lipids)[duplicated_columns]
print(duplicated_column_names)
duplicated_column_names <- data.frame(duplicated_columns = duplicated_column_names)
write.csv(duplicated_column_names, "outputs/error_files/duplicated_columns.csv", row.names = FALSE)

################################################################################
######################## LIPIDS TESTED MANIPULATION ############################
################################################################################

lipids_tested <- read.csv(lipids_tested_file_name, header = T)
lipids_tested <- pivot_longer(lipids_tested,
    cols = everything(),
    names_to = "family",
    values_to = "lipid"
  )
lipids_tested <- lipids_tested[order(lipids_tested$family), ]
lipids_tested <- lipids_tested[!(is.na(lipids_tested$lipid) | lipids_tested$lipid == ""), ]
lipids_tested$lipid <- make.names(lipids_tested$lipid)

################################################################################
####################### CHECKING FOR MISMATCHED COLUMNS ########################
######################## I.E. COLUMNS TESTED THAT THERE ########################
########################## WAS NO RESULTANT DATA FOR ###########################
################################################################################

raw_data_columns <- colnames(raw_data_lipids[2:ncol(raw_data_lipids)])
lipids <- lipids_tested$lipid
unmatched_cols <- setdiff(raw_data_columns, lipids)
number_unmatched <- length(unmatched_cols)
unmatched_cols <- data.frame(Unmatched_columns = unmatched_cols)
write.csv(unmatched_cols, "outputs/error_files/unmatched_columns.csv", row.names = FALSE)

################################################################################
######################## SUBSET DATA BY LIPID FAMILY ###########################
################################################################################

lipid_families <- split(lipids_tested$lipid, lipids_tested$family)
raw_data_by_family <- list()

for(family in names(lipid_families)) {
  cols <- lipid_families[[family]]
  cols <- intersect(cols, colnames(raw_data_lipids))
  raw_data_by_family[[family]] <- raw_data[, cols, drop = FALSE]
}

for(family in names(raw_data_by_family)) {
  assign(paste0("raw_data_", family), raw_data_by_family[[family]])
}

raw_data_names <- ls(pattern = "^raw_data_")

for (name in raw_data_names) {
  df <- get(name)
  if (is.data.frame(df) && ncol(df) == 0) {
    message("Removing empty data frame: ", name)
    rm(list = name, envir = .GlobalEnv)
  }
}

################################################################################
################### CREATING FOLDERS FOR ALL LIPID FAMILIES ####################
################################################################################

raw_data_names <- ls(pattern = "^raw_data_")
raw_data_names <- raw_data_names[sapply(raw_data_names, function(x) {
  df <- get(x)
  is.data.frame(df) && all(sapply(df, is.numeric))
})]
lipid_families_folder <- file.path("outputs", "lipid_families")

if (!dir.exists(lipid_families_folder)) {
  dir.create(lipid_families_folder)
}

for (name in raw_data_names) {
  cat("Creating folder for:", name, "\n")
  lipid_family <- sub("^raw_data_", "", name)
  folder_path <- file.path("outputs/lipid_families", lipid_family)
  if (!dir.exists(folder_path)) {
    dir.create(folder_path)
  }
}

################################################################################
########################## END OF DATA ORGANISATION ############################
################################################################################
################################################################################
########################## END OF DATA ORGANISATION ############################
################################################################################
################################################################################
########################## END OF DATA ORGANISATION ############################
################################################################################

################################################################################
############################## TOTAL DATA ANOVA ################################
################################################################################

raw_data_lipids$Group <- as.factor(raw_data_lipids$Group)

################################################################################
######################### HISTOGRAMS FOR NORMALITY #############################
################################################################################

for (df_name in raw_data_names) {
  lipid_family <- sub("^raw_data_", "", df_name)
  folder_path <- file.path("outputs", "lipid_families", 
                           lipid_family, "histograms")
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = T)
  }
  cat("Processing histograms for:", df_name, "\n")
  df <- get(df_name)
  for (col in setdiff(names(df), "groups")) {
    png(file.path(folder_path, paste0(col, "_hist.png")))
    hist(df[[col]],
         main = paste("Histogram of", col),
         xlab = col)
    dev.off()
  }
}

























































































