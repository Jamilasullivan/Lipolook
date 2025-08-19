################################################################################
############################ VARIABLES TO FILL #################################
################################################################################

working_directory <- "C:/Users/jamsu/OneDrive - Cardiff University/University/Masters/Big Data Biology/Modules/Dissertation/R Programme/Dissertation/test_data"

raw_data_file_name <- "raw_data_1.csv"

lipids_tested_file_name <- "lipids_tested_1.csv"

adjustment_method <- "fdr" # options are "holm" = Holm–Bonferroni
                        #             "hochber" = Hochberg
                        #             "hommel" = Hommel
                        #             "bonferroni" = Bonferroni
                        #             "BH" or "fdr" = Benjamini–Hochberg
                        #             "BY" = Benjamini–Yekutieli
                        #             "none" = No adjustment

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

all_raw_data <- ls(pattern = "^raw_data_")

for (name in all_raw_data) {
  df <- get(name)
  if (is.data.frame(df) && ncol(df) == 0) {
    message("Removing empty data frame: ", name)
    rm(list = name, envir = .GlobalEnv)
    all_raw_data <- setdiff(all_raw_data, name)  # remove from vector too
  }
}

raw_data_names <- Filter(function(x) {
  obj <- get(x)
  
  # Must be a data frame
  if (!is.data.frame(obj)) {
    message("Excluded for not being a data frame: ", x)
    return(FALSE)
  }
  
  # Must be fully numeric
  if (!all(sapply(obj, is.numeric))) {
    message("Excluded for non-numeric columns: ", x)
    return(FALSE)
  }
  
  TRUE
}, all_raw_data)

raw_data_names

################################################################################
########### CREATING FOLDERS FOR ALL LIPID FAMILIES AND CATEGORIES #############
################################################################################

category_mapping <- read.csv("lipid_categories_1.csv", stringsAsFactors = FALSE)

category_mapping$Family_clean <- make.names(category_mapping$Family)
category_mapping$Category_clean <- make.names(category_mapping$Category)

top_level_dir <- file.path("outputs", "lipid_categories")

if (!dir.exists(top_level_dir)) {
  dir.create(top_level_dir, recursive = TRUE, showWarnings = FALSE)
}

count <- 1

for (name in raw_data_names) {
  lipid_family <- make.names(sub("^raw_data_", "", name))
  category <- category_mapping$Category_clean[category_mapping$Family_clean == lipid_family][1]
  
  if (is.na(category) || length(category) == 0) {
    message("No category found for family: ", lipid_family)
    next
  }
  
  # Build path: outputs/lipid_categories/category/family
  folder_path <- file.path(top_level_dir, category, lipid_family)
  
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE, showWarnings = FALSE)
    message(count, ". Folder created: ", folder_path)
    count <- count + 1
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
######################## ADDING GROUPS BACK TO DFS #############################
################################################################################

groups <- raw_data[,3]

for (name in raw_data_names) {
  cat("Adding 'groups' column to:", name, "\n")
  df <- get(name)
  df <- cbind(groups, df)
  assign(name, df)
} 

################################################################################
######################## CALCULATING GROUP AVERGAGES ###########################
################################################################################

top_level_dir <- file.path("outputs", "lipid_categories")

for (name in raw_data_names) {
  # Extract lipid family
  lipid_family <- make.names(sub("^raw_data_", "", name))
  
  # Look up category
  category <- category_mapping$Category_clean[category_mapping$Family_clean == lipid_family][1]
  
  if (is.na(category) || length(category) == 0) {
    message("No category found for family: ", lipid_family)
    next
  }
  
  # Use existing folder path
  folder_path <- file.path(top_level_dir, category, lipid_family)
  
  cat("Saving averages for:", name, "to folder:", folder_path, "\n")
  
  df_avg <- aggregate(. ~ groups, data = get(name), FUN = mean)
  csv_file <- file.path(folder_path, paste0(lipid_family, "_avg.csv"))
  write.csv(df_avg, file = csv_file, row.names = FALSE)
}

################################################################################
################################# ANOVA ########################################
################################################################################

## for individual lipid families ###############################################

for (name in raw_data_names) {
  cat("Running ANOVA for:", name, "\n")
  df <- get(name)
  lipid_columns <- setdiff(names(df), "groups")
  p_values <- numeric(length(lipid_columns))
  names(p_values) <- lipid_columns
  
  for (lipid in lipid_columns) {
    formula <- as.formula(paste(lipid, "~ groups"))
    aov_result <- aov(formula, data = df)
    summary_result <- summary(aov_result)
    p_values[lipid] <- summary_result[[1]][["Pr(>F)"]][1]
  }
  
  p_adjusted <- p.adjust(p_values, method = adjustment_method)

  output_df <- data.frame(
    lipid = lipid_columns,
    p_value = p_values,
    p_adjusted = p_adjusted
  )
  output_df$significance <- cut(output_df$p_adjusted,
                                breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                                labels = c("***", "**", "*", "NO"),
                                right = TRUE
  )
  output_df <- output_df[order(output_df$p_adjusted, decreasing = FALSE), ]
  lipid_family <- sub("^raw_data_", "", name)
  folder_path <- file.path("outputs", "lipid_families", lipid_family)
  
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
  }
  
  csv_file <- file.path(folder_path, paste0(lipid_family, "_anova_pvalues.csv"))
  write.csv(output_df, file = csv_file, row.names = FALSE)
  cat("Saved results to:", csv_file, "\n")
}

## for all lipids combined #####################################################

raw_data_lipids$Group <- as.factor(raw_data_lipids$Group)
lipid_columns <- names(raw_data_lipids)[-1]
p_values <- sapply(lipid_columns, function(lipid) {
  formula <- as.formula(paste(lipid, "~ groups"))
  summary(aov(formula, data = raw_data_lipids))[[1]][["Pr(>F)"]][1]
})
p_adjusted <- p.adjust(p_values, method = adjustment_method)
output_df <- data.frame(
  lipid = lipid_columns,
  p_value = p_values,
  p_adjusted = p_adjusted,
  significance = ifelse(p_adjusted <= 0.001, "***",
                        ifelse(p_adjusted <= 0.01, "**",
                               ifelse(p_adjusted <= 0.05, "*", "NO")))
)
output_df <- output_df[order(output_df$p_adjusted, decreasing = FALSE), ]
print(output_df)
write.csv(output_df, "outputs/total_lipids/total_lipid_anova_results.csv", row.names = FALSE)

num_vstrong <- sum(output_df$significance == "***")
num_strong <- sum(output_df$significance == "**")
num_significant <- sum(output_df$significance == "*")
num_not_signif <- sum(output_df$significance == "NO")

cat("Number of very highly significant lipids (***):", num_vstrong, "\n")
cat("Number of highly significant lipids (**):", num_strong, "\n")
cat("Number of significant lipids (*):", num_significant, "\n")
cat("Number of insignificant lipids (NO):", num_not_signif, "\n")

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






































































































