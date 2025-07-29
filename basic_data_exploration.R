## packages ####################################################################

library(tidyr)
library(dplyr)

## Data given in .xlsx format ##################################################

setwd("C:/Users/jamsu/OneDrive - Cardiff University/University/Masters/Big Data Biology/Modules/Dissertation/R Programme/Dissertation/test_data")

################################################################################
## data provided to the client #################################################
################################################################################

raw_data <- read.csv("raw_data_1.csv", header = T)
names(raw_data)[c(1,2,3,4,5)] <- c("Filename", "CMaLL Sample Name", "Group", "code", "Cell Count")
raw_data <- raw_data[-1,]
raw_data$code <- NULL
raw_data$`Cell Count` <- NULL
#View(raw_data)
summary(raw_data)
str(raw_data) # they were all characters here
# Convert columns 4 to end (since the first 3 are metadata)
raw_data[, 4:ncol(raw_data)] <- lapply(raw_data[, 4:ncol(raw_data)], function(x) as.numeric(as.character(x)))
sum(is.na(raw_data[, 4:ncol(raw_data)])) # no NA values introduced by coercion
str(raw_data)
summary(raw_data)

summary(raw_data[, 4:10]) # summaries for a subset of data

## setting row names to sample names (not decided long term) and deleting all columns other than lipids tested.

raw_data_lipids <- raw_data

rownames(raw_data_lipids) <- raw_data_lipids[[2]]

raw_data_lipids <- raw_data_lipids[, -(1:2)] # to 2 is keeping grouping, to 3 is without grouping

duplicated_columns <- duplicated(as.list(raw_data_lipids)) # creates logical list of whether any columns are completely duplicated values
any(duplicated_columns) # tells you if any in the list are true

################################################################################
## lipids tested for the client ################################################
################################################################################

lipids_tested <- read.csv("lipids_tested_1.csv", header = T)
#View(lipids_tested)

lipids_tested <- lipids_tested %>% 
  pivot_longer(
    cols = everything(),
    names_to = "family",
    values_to = "lipid"
  ) %>% 
  filter(!is.na(lipid))

lipids_tested <- lipids_tested[order(lipids_tested$family), ] # ordering the lipid families alphabetically

summary(lipids_tested) # 4318 rows with empty lipid cells

lipids_tested <- lipids_tested[!(is.na(lipids_tested$lipid) | lipids_tested$lipid == ""), ] # getting rid of the rows without relevant lipids attached (i.e. blank cells for lipid names)

lipids_tested$lipid <- make.names(lipids_tested$lipid) # changing the (): characters to . to match the syntax of the column names in the other data frame.

#View(lipids_tested)
summary(lipids_tested) # 1602 rows without empty lipid cells

## checking for mismatched columns ############################################

#all_cols <- unique(unlist(lapply(raw_data_by_family, colnames))) 
#number_of_lipids_saved <- length(all_cols)  # total number of unique column names. Should be the same as the number of variables in the raw_data_lipids. If it's not then there is an issue with mismatched columns somewhere and the below script will help with figuring out where those are.

raw_data_columns <- colnames(raw_data_lipids)
lipids <- lipids_tested$lipid

unmatched_cols <- setdiff(raw_data_columns, lipids)
number_unmatched <- length(unmatched_cols)
print(unmatched_cols[1:number_unmatched])

## subset data by lipid family #################################################

lipid_families <- split(lipids_tested$lipid, lipids_tested$family) # split the lipids by their family and save them as such

raw_data_by_family <- list() # creating an empty list to store the separated data frames

for(family in names(lipid_families)) {
  # Get lipid columns for this family
  cols <- lipid_families[[family]]
  
  # Only keep columns that actually exist in raw_data
  cols <- intersect(cols, colnames(raw_data_lipids))
  
  # Subset raw_data by these columns
  raw_data_by_family[[family]] <- raw_data[, cols, drop = FALSE]
} # the column numbers (489 added) do not match the number of columns in the raw data (500 columns) because of the mismatched columns previously detected. This means that this program will only separate columns in raw data that are identically found in the list of tested lipids. Anything named any differently will be ignored. Therefore, data could be lost this way and is something to be very careful of.

for(family in names(raw_data_by_family)) {
  assign(paste0("raw_data_", family), raw_data_by_family[[family]])
} # creating different data frames of the separated 

## separate groups out #########################################################

groups <- data.frame(raw_data[3])










