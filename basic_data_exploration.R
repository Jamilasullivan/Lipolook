## packages ####################################################################

library(tidyr)
library(dplyr)

## Data given in .xlsx format ##################################################

setwd("C:/Users/jamsu/OneDrive - Cardiff University/University/Masters/Big Data Biology/Modules/Dissertation/R Programme/Dissertation/test_data")

################################################################################
## data provided to the client #################################################
################################################################################

raw_data <- read.csv("raw_data_1.csv", header = T) # reading in the raw data

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

################################################################################
## lipids tested for the client ################################################
################################################################################

lipids_tested <- read.csv("lipids_tested_1.csv", header = T) # reading in the list of lipids tested
#View(lipids_tested)

lipids_tested <- lipids_tested %>% 
  pivot_longer(
    cols = everything(),
    names_to = "family",
    values_to = "lipid"
  ) %>% 
  filter(!is.na(lipid)) # putting all column values (lipid) into rows with their associated column name (family)

lipids_tested <- lipids_tested[order(lipids_tested$family), ] # ordering the lipid families alphabetically

summary(lipids_tested) # 4318 rows with empty lipid cells

lipids_tested <- lipids_tested[!(is.na(lipids_tested$lipid) | lipids_tested$lipid == ""), ] # getting rid of the rows without relevant lipids attached (i.e. blank cells for lipid names)

lipids_tested$lipid <- make.names(lipids_tested$lipid) # changing the (): characters to . to match the syntax of the column names in the other data frame.

#View(lipids_tested)
summary(lipids_tested) # 1602 rows without empty lipid cells

## checking for mismatched columns #############################################

#all_cols <- unique(unlist(lapply(raw_data_by_family, colnames))) 
#number_of_lipids_saved <- length(all_cols)  # total number of unique column names. Should be the same as the number of variables in the raw_data_lipids. If it's not then there is an issue with mismatched columns somewhere and the below script will help with figuring out where those are.

raw_data_columns <- colnames(raw_data_lipids) # saving all of the column names (lipids)
lipids <- lipids_tested$lipid # saving the names of all of the lipids tested.

unmatched_cols <- setdiff(raw_data_columns, lipids) # looking at how the two above lists match by outputting which ones didn't
number_unmatched <- length(unmatched_cols) # asking how many columns didn't match
print(unmatched_cols[1:number_unmatched]) # to tell me what the mismatched columns are

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

## remove data frames with 0 variables before averaging ########################
## it causes issues with loops to keep them there and they are unnecessary #####

raw_data_names <- ls(pattern = "^raw_data_") # making a list of all the data frames I have with names starting 'raw_data_'

for (name in raw_data_names) {
  df <- get(name)
  
  if (is.data.frame(df) && ncol(df) == 0) {
    message("Removing empty data frame: ", name)
    rm(list = name, envir = .GlobalEnv)
  }
} # checking for empty data frames and removing them. This loops says, check for every name in the list I've made, get the data frame that's related to it, if the data frame has no columns then give a message that it's empty and then remove it. 

################################################################################
############ Primary data exploration ##########################################
################################################################################

## my idea is to firstly do this for one single lipid family to see what I'd like as an output and then loop this through all for the programme.

## below is the test using just 

## FIRSTLY: COMPARING THE AVERAGE VALUES OF EACH COLUMN (LIPID) WITHIN THE FAMILY

avg_row <- colMeans(raw_data_Triacyl.glycerols, na.rm = TRUE) # make a data frame with the average values

raw_data_Triacyl.glycerols <- rbind(raw_data_Triacyl.glycerols, Average = avg_row) # adding the average row to the original data frame

avg_values <- as.numeric(raw_data_Triacyl.glycerols["Average", ]) # removing the average values to plot

names(avg_values) <- colnames(raw_data_Triacyl.glycerols)

pal <- colorRampPalette(c("mistyrose", "darkred")) # creating the colour palette for bars
bar_colours <- pal(length(avg_values))[rank(avg_values)] # 

bp <- barplot(
  avg_values,
  main = "Average Lipid Values",
  horiz = TRUE,
  las = 2,                # Rotate x-axis labels
  cex.names = 0.7,        # Shrink label size if needed
  col = bar_colours
)

text(
  x = avg_values + max(avg_values) * 0.02,  # position slightly to the right
  y = bp,
  labels = round(avg_values, 2),            # round to 2 decimals
  cex = 0.7,
  adj = 0
)

## SECONDLY: COMPARING AVERAGE VALUES BETWEEN EXPERIMENTAL GROUPS




