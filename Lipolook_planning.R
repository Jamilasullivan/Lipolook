## packages ####################################################################

library(tidyr)
library(dplyr)
library(ggplot2)

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

## FIRSTLY: COMPARING THE AVERAGE VALUES OF EACH COLUMN (LIPID) WITHIN THE FAMILY #########################################################################

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

# I now want to make a forest plot for all of these TGs too

raw_data_Triacyl.glycerols <- raw_data_Triacyl.glycerols[1:57,]

raw_data_Triacyl.glycerols$Sample <- rownames(raw_data_Triacyl.glycerols)

raw_data_Triacyl.glycerols_long <- raw_data_Triacyl.glycerols %>%
  pivot_longer(
    cols = -Sample, 
    names_to = "Lipid", 
    values_to = "Value"
  )

TG_summary_stats <- raw_data_Triacyl.glycerols_long %>%
  group_by(Lipid) %>%
  summarise(
    mean = mean(Value, na.rm = TRUE),
    sd = sd(Value, na.rm = TRUE),
    n = n(),
    se = sd / sqrt(n),
    lower = mean - 1.96 * se,
    upper = mean + 1.96 * se
  )

ggplot(TG_summary_stats, aes(x = mean, y = reorder(Lipid, mean))) +
  geom_point(color = "red", size = 3) +             # mean points
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.3, color = "red") +  # 95% CI
  labs(
    x = "Mean Lipid Value with 95% CI",
    y = "Lipid",
    title = "Forest Plot: Lipid Distribution Across Samples"
  ) +
  theme_minimal() # this plots all 63 lipids in the TG family

# I'd like to now just plot the top 10%

TG_summary_stats <- TG_summary_stats %>%
  mutate(ci_width = upper - lower)

cutoff <- quantile(TG_summary_stats$ci_width, 0.9)

top_10_percent_CI <- TG_summary_stats %>%
  filter(ci_width >= cutoff)

ggplot(top_10_percent_CI, aes(x = mean, y = reorder(Lipid, mean))) +
  geom_point(color = "darkred", size = 3) +                        # mean points
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.3,    # horizontal error bars
                 color = "darkred") +
  labs(
    x = "Mean Lipid Value with 95% CI",
    y = "Lipid",
    title = "Forest Plot: Lipids with Top 10% Widest 95% CI"
  ) +
  theme_minimal()

## SECONDLY: COMPARING AVERAGE VALUES BETWEEN EXPERIMENTAL GROUPS ##############

# this is dependent on having the groups back in the data frame

groupings <- raw_data_lipids[,1] # making a list of all the groupings associated with the data

raw_data_Triacyl.glycerols <- raw_data_Triacyl.glycerols[1:57,] # removing the average column so that the next step of binding has the same number of columns to work with 

TG_grouping <- cbind(groupings,raw_data_Triacyl.glycerols) # adding the groupings to the TG data specifically

# now ask R for averages for each group to compare. Use aggregate which splits data into specified groups, applies a function and then gives back a new data frame with the results of the function

TG_group_averages <- aggregate(. ~ groupings, data = TG_grouping, FUN = mean, na.rm = TRUE) # Asks for 'all other columns' against column 1 (groupings) in TG_grouping to be processed for their mean. This leaves 19 columns, one for each group. 

colnames(TG_group_averages)[1] <- "Group" # renaming the column to 'Group'

# I now want to compare these results

# firstly, I want to make a bar chart for all of them individually


################################################################################
########### WORKING ON A LOOP FOR ALL LIPID FAMILIES ###########################
################################################################################

raw_data_names <- ls(pattern = "^raw_data_") # makes sure we're analysing everything that has data in the final data frame

raw_data_names <- raw_data_names[sapply(raw_data_names, function(x) {
  df <- get(x)
  is.data.frame(df) && all(sapply(df, is.numeric))
})]

if (!dir.exists("outputs")) {
  dir.create("outputs")
} # makes a folder called 'outputs'

lipid_families_folder <- file.path("outputs", "lipid_families")

if (!dir.exists(lipid_families_folder)) {
  dir.create(lipid_families_folder)
} # makes an extra folder in 'outputs' called 'lipid families'


for (name in raw_data_names) {
  cat("Creating folder for:", name, "\n")
  lipid_family <- sub("^raw_data_", "", name)
  folder_path <- file.path("outputs/lipid_families", lipid_family)
  if (!dir.exists(folder_path)) {
    dir.create(folder_path)
  }
} # creates a folder in 'outputs/lipid_families' for all information relating to each of the lipid families




#### FOR LOOP FOR COMPARING THE AVERAGES OF ALL LIPIDS WITHIN A FAMILY #########

for (df_name in raw_data_names) {
  
  lipid_family <- sub("^raw_data_", "", df_name)
  folder_path <- file.path("outputs", "lipid_families", lipid_family)
  
  ## setting message
  cat("Processing:", df_name, "\n")
  
  ## Just in case folder doesn't exist, you can optionally check or skip
  if (!dir.exists(folder_path)) {
    warning(paste("Folder does not exist:", folder_path))
    next  # skip to next iteration if folder missing
  }
  
  df <- get(df_name)
  
  # Calculate average values for each lipid
  avg_values <- colMeans(df, na.rm = TRUE)
  
  plot_df <- data.frame(
    Lipid = names(avg_values),
    Average = as.numeric(avg_values)
  )
  
  p <- ggplot(plot_df, aes(x = reorder(Lipid, Average), y = Average, fill = Average)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_gradient(low = "pink", high = "red") +
    labs(title = paste("Average Lipid Values -", lipid_family),
         x = "Lipid",
         y = "Average Value") +
    theme_minimal()
  
  # Save the plot inside the existing folder
  plot_file <- file.path(folder_path, paste0(lipid_family, "_barplot.png"))
  ggsave(filename = plot_file, plot = p, width = 8, height = 4)
}














































