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
############################## WAS NO DATA FOR #################################
################################################################################



















