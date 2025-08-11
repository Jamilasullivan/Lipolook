################################################################################
#################### PARAMETERS TO SET #########################################
################################################################################

working_directory <- "C:/Users/jamsu/OneDrive - Cardiff University/University/Masters/Big Data Biology/Modules/Dissertation/R Programme/Dissertation/test_data"

raw_data <- "raw_data_1.csv"






################################################################################
################ START OF THE SCRIPT ###########################################
################################################################################


## packages ####################################################################

library(tidyr)
library(dplyr)

## Data given in .xlsx format ##################################################

setwd(working_directory)

################################################################################
## data provided to the client #################################################
################################################################################

raw_data <- read.csv(raw_data, header = T) # reading in the raw data

names(raw_data)[c(1,2,3,4,5)] <- c("Filename", "CMaLL Sample Name", "Group", "code", "Cell Count") # renaming columns

raw_data <- raw_data[-1,] # deleting the first row 

raw_data$code <- NULL # deleting 'code' column

raw_data$`Cell Count` <- NULL # deleting 'Cell count' column

#View(raw_data) # viewing raw data 

summary(raw_data) # summary of results

str(raw_data) # looking at data structure. They were all characters here.

raw_data[, 4:ncol(raw_data)] <- lapply(raw_data[, 4:ncol(raw_data)], function(x)
  as.numeric(as.character(x))) # changing columns 4-n in raw data to numeric variables

sum(is.na(raw_data[, 4:ncol(raw_data)])) # checking that no NA values were introduced by coercion

str(raw_data) # checking the structure again

summary(raw_data) # summarising the data 

summary(raw_data[, 4:10]) # summaries for a subset of data

##### setting row names to sample names (not decided long term) and deleting all columns other than lipids tested. ##############################################

raw_data_lipids <- raw_data # copying the data frame to manipulate

rownames(raw_data_lipids) <- raw_data_lipids[[2]] # renaming the row names as the sample names

raw_data_lipids <- raw_data_lipids[, -(1:2)] # Deleting columns 1-2/3 to just keep the lipid data and remove unnecessary metadata. To 2 is keeping grouping, to 3 is without grouping.

duplicated_columns <- duplicated(as.list(raw_data_lipids)) # creates logical list of whether any columns are completely duplicated values

any(duplicated_columns) # tells you if any in the list are true