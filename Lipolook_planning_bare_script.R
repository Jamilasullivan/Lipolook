################################################################################
########################## 1. VARIABLES TO FILL ################################
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

groups <- ## what is the current name of your column containing groupings in raw data?
  
control_group <- "A" ## the name of your control group in your raw data file. what is your control group called?
  
glm_variables <- c("","","") ## include the names of columns with metadata to be included

################################################################################
########################## 2. PACKAGES TO INSTALL ##############################
################################################################################

packages <- c("moments", "tidyr", "dplyr", "ggplot2", "stats")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

lapply(packages, library, character.only = TRUE)
if (!require(dunn.test)) install.packages("dunn.test")

library(moments)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stats)
library(dunn.test)

################################################################################
######################## 3. SETTING WORK DIRECTOTY #############################
################################################################################

setwd(working_directory)

################################################################################
###################### 4. CREATING DIRECTORY STRUCTURE #########################
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
######################### 5. INITIAL DATA TIDYING ##############################
################################################################################

raw_data <- read.csv(raw_data_file_name, header = T)
metadata <- raw_data[ , grep("^X", names(raw_data))]
colnames(metadata) <- as.character(unlist(metadata[1, ]))
metadata <- metadata[-1,]
groups <- metadata[,"Group"]
raw_data <- raw_data[ , !grepl("^X", names(raw_data))]
concentration <- raw_data[1,1]
raw_data <- raw_data[-1,]
str(raw_data)
raw_data[] <- lapply(raw_data, function(x) as.numeric(as.character(x)))
str(raw_data)
sum(is.na(raw_data[,ncol(raw_data)]))

################################################################################
######################### 6. RAW DATA MANIPULATION #############################
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
####################### 7. LIPIDS TESTED MANIPULATION ##########################
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
###################### 8. CHECKING FOR MISMATCHED COLUMNS ######################
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
####################### 9. SUBSET DATA BY LIPID FAMILY #########################
################################################################################

sanitize_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)  # replace non-alphanumeric with underscore
  x <- gsub("_+", "_", x)             # collapse multiple underscores
  x <- gsub("^_|_$", "", x)           # trim leading/trailing underscores
  x
}

lipids_tested$Family_clean <- sanitize_name(lipids_tested$family)
lipid_families <- split(lipids_tested$lipid, lipids_tested$Family_clean)
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
    all_raw_data <- setdiff(all_raw_data, name)  
  }
}

raw_data_names <- Filter(function(x) {
  obj <- get(x)
  
  if (!is.data.frame(obj)) {
    message("Excluded for not being a data frame: ", x)
    return(FALSE)
  }
  
  if (!all(sapply(obj, is.numeric))) {
    message("Excluded for non-numeric columns: ", x)
    return(FALSE)
  }

  TRUE
}, all_raw_data)

#1
raw_data_CPEs <- raw_data_Ceramide_phospho_ethanolamines
rm(raw_data_Ceramide_phospho_ethanolamines)
#2
raw_data_Chol_CEs <- raw_data_Cholesterol_cholesteryl_esters
rm(raw_data_Cholesterol_cholesteryl_esters)
#3
raw_data_LPEs <- raw_data_Lysophosphatidy_ethanolamines
rm(raw_data_Lysophosphatidy_ethanolamines)

raw_data_names <- c(raw_data_names, c("raw_data_CPEs", "raw_data_Chol_CEs", "raw_data_LPEs"))
raw_data_names <- raw_data_names[!raw_data_names %in% c("raw_data_Ceramide_phospho_ethanolamines", "raw_data_Cholesterol_cholesteryl_esters", "raw_data_Lysophosphatidy_ethanolamines")]

raw_data_names

################################################################################
######### 10. CREATING FOLDERS FOR ALL LIPID FAMILIES AND CATEGORIES ###########
################################################################################

category_mapping <- read.csv("lipid_categories_1.csv", stringsAsFactors = FALSE)

category_mapping$Family_clean   <- sanitize_name(category_mapping$Family)
category_mapping$Category_clean <- sanitize_name(category_mapping$Category)

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
###################### 11. ADDING GROUPS BACK TO DFS ###########################
################################################################################

for (name in raw_data_names) {
  cat("Adding 'groups' column to:", name, "\n")
  df <- get(name)
  df <- cbind(groups, df)
  assign(name, df)
} # only run once!!

for (name in raw_data_names) {
  lipid_family <- make.names(sub("^raw_data_", "", name))
  category <- category_mapping$Category_clean[category_mapping$Family_clean == lipid_family][1]
  
  if (is.na(category) || length(category) == 0) {
    message("No category found for family: ", lipid_family)
    next
  }
  
  folder_path <- file.path(top_level_dir, category, lipid_family)
  cat("Saving raw data for:", name, "\n")
  df <- get(name)
  write.csv(df, file = file.path(folder_path, paste0(lipid_family, "_raw_data.csv")), row.names = FALSE)
}

################################################################################
###################### 12. CALCULATING GROUP AVERGAGES #########################
################################################################################

top_level_dir <- file.path("outputs", "lipid_categories")

for (name in raw_data_names) {
  lipid_family <- make.names(sub("^raw_data_", "", name))
  category <- category_mapping$Category_clean[category_mapping$Family_clean == lipid_family][1]
  
  if (is.na(category) || length(category) == 0) {
    message("No category found for family: ", lipid_family)
    next
  }
  
  folder_path <- file.path(top_level_dir, category, lipid_family)
  cat("Saving averages for:", name, "to folder:", folder_path, "\n")
  df_avg <- aggregate(. ~ groups, data = get(name), FUN = mean)
  csv_file <- file.path(folder_path, paste0(lipid_family, "_avg.csv"))
  write.csv(df_avg, file = csv_file, row.names = FALSE)
}

################################################################################
####################### 13. HISTOGRAMS FOR NORMALITY ###########################
################################################################################

top_level_dir <- file.path("outputs", "lipid_categories")

for (name in raw_data_names) {
  lipid_family <- sub("^raw_data_", "", name)
  category <- category_mapping$Category_clean[category_mapping$Family_clean == lipid_family][1]
  
  if (is.na(category) || length(category) == 0) {
    message("No category found for family: ", lipid_family)
    next
  }
  
  folder_path <- file.path(top_level_dir, category, lipid_family)
  
  hist_folder <- file.path(folder_path, "histograms")
  if (!dir.exists(hist_folder)) {
    dir.create(hist_folder, recursive = TRUE, showWarnings = FALSE)
    message("Created histograms subfolder for: ", lipid_family)
  }
  
  cat("Processing histograms for:", name, "\n")
  df <- get(name)
  
  for (col in setdiff(names(df), "groups")) {
    png(file.path(hist_folder, paste0(col, "_hist.png")))
    hist(df[[col]],
         main = paste("Histogram of", col),
         xlab = col,
         breaks = 30)
    dev.off()
  }
}

################################################################################
############### 14. NORMALITY OUTPUT FOR ALL SAMPLES COMBINED ##################
################################################################################

## across all samples for each lipid ###########################################

top_level_dir <- file.path("outputs", "lipid_categories")

for (name in raw_data_names) {
  lipid_family <- sub("^raw_data_", "", name)
  category <- category_mapping$Category_clean[category_mapping$Family_clean == lipid_family][1]
  
  if (is.na(category) || length(category) == 0) {
    message("No category found for family: ", lipid_family)
    next
  }
  
  folder_path <- file.path(top_level_dir, category, lipid_family)  
  df <- get(name)
  lipid_columns <- names(df)[-1]
  distribution_summary <- data.frame(
    lipid = lipid_columns,
    skewness = numeric(length(lipid_columns)),
    kurtosis = numeric(length(lipid_columns)),
    shapiro_p = numeric(length(lipid_columns)),
    normality = character(length(lipid_columns)),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(lipid_columns)) {
    lipid <- lipid_columns[i]
    values <- df[[lipid]]  
    distribution_summary$skewness[i] <- skewness(values)
    distribution_summary$kurtosis[i] <- kurtosis(values)
    
    if (length(values) >= 3 & length(values) <= 5000) {
      shapiro_res <- shapiro.test(values)
      distribution_summary$shapiro_p[i] <- shapiro_res$p.value
      distribution_summary$normality[i] <- ifelse(shapiro_res$p.value > 0.05, "Normal", "Non-normal")
    } else {
      distribution_summary$shapiro_p[i] <- NA
      distribution_summary$normality[i] <- "NA"
    }
  }
  
  distribution_summary$shapiro_p_adj <- p.adjust(distribution_summary$shapiro_p, method = adjustment_method)
  distribution_summary$normality <- ifelse(distribution_summary$shapiro_p_adj > 0.05, "Normal", "Non-normal")
  distribution_summary <- distribution_summary[, c("lipid", "skewness", "kurtosis", "shapiro_p", "shapiro_p_adj", "normality")]
  cat("Saving distribution summary for:", name, "\n")
  
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE, showWarnings = FALSE)
  }
  
  dis_csv_file <- file.path(folder_path, paste0(lipid_family, "_distribution.csv"))
  write.csv(distribution_summary, file = dis_csv_file, row.names = FALSE)
}

## number of all lipids that reached normality with shapiro-wilk test (plus correction for multiple testing)

top_level_dir <- file.path("outputs", "lipid_categories")
all_distribution_csv_files <- list.files(top_level_dir, pattern = "_distribution\\.csv$", full.names = TRUE, recursive = TRUE)
all_lipids_distribution <- list()

for (file in all_distribution_csv_files) {
  df <- read.csv(file, header = TRUE)
  summary_df <- df[, c(1, ncol(df))]
  summary_df$family <- basename(file)
  all_lipids_distribution[[length(all_lipids_distribution) + 1]] <- summary_df
}

combined_summary <- bind_rows(all_lipids_distribution)
combined_summary <- combined_summary[-3]
write.csv(combined_summary, file = file.path(top_level_dir,"..","total_lipids", "combined_lipid_normality.csv"), row.names = FALSE)
distribution_counts <- table(combined_summary$normality)
distribution_counts_df <- as.data.frame(distribution_counts)
distribution_total <- sum(distribution_counts)
distribution_counts_df$percent <- round(100 * distribution_counts_df$Freq / distribution_total, 1)
distribution_counts_df
write.csv(distribution_counts_df, file = file.path(top_level_dir,"..","total_lipids", "combined_lipid_normality_summary.csv"), row.names = FALSE)

################################################################################
################# 15. NORMALITY OUTPUT FOR INDIVIDUAL GROUPS ###################
################################################################################

top_level_dir <- file.path("outputs", "lipid_categories")

for (name in raw_data_names) {
  lipid_family <- sub("^raw_data_", "", name)
  category <- category_mapping$Category_clean[category_mapping$Family_clean == lipid_family][1]
  
  if (is.na(category) || length(category) == 0) {
    message("No category found for family: ", lipid_family)
    next
  }
  
  folder_path <- file.path(top_level_dir, category, lipid_family)  
  df <- get(name)
  lipid_columns <- names(df)[-1]
  group_col <- names(df)[1]
  groups <- unique(df[[group_col]])
  
  per_group_results <- list()
  
  for (i in seq_along(lipid_columns)) {
    lipid <- lipid_columns[i]
    
    group_p <- c()
    group_normality <- c()
    
    for (g in groups) {
      vals_g <- df[df[[group_col]] == g, lipid]
      if (length(vals_g) >= 3 & length(vals_g) <= 5000) {
        if (length(unique(vals_g)) > 1) {
          p <- shapiro.test(vals_g)$p.value
          group_p <- c(group_p, p)
          group_normality <- c(group_normality, ifelse(p > 0.05, "Normal", "Non-normal"))
        } else {
          group_p <- c(group_p, NA)
          group_normality <- c(group_normality, "All identical")
        }
      } else {
        group_p <- c(group_p, NA)
        group_normality <- c(group_normality, NA)
      }
    }
    
    # adjust all p-values for this lipid (not inside the loop)
    group_p_adj <- p.adjust(group_p, method = adjustment_method)
    group_normality <- ifelse(group_p_adj > 0.05, "Normal", "Non-normal")
    
    per_group_results[[lipid]] <- data.frame(
      group = groups, 
      p_value = group_p, 
      p_value_adj = group_p_adj, 
      normality = group_normality
    )
  }
  
  # save all per-group results
  for (lipid in names(per_group_results)) {
    write.csv(
      per_group_results[[lipid]], 
      file = file.path(folder_path, paste0(lipid, "_distribution_groups.csv")), 
      row.names = FALSE
    )
    message("Saving distribution per group for: ", lipid_family, " -> ", lipid)
  }
}

## combining normality of all lipids across groups #############################

top_level_dir <- file.path("outputs", "lipid_categories")

all_files <- list.files(top_level_dir, pattern = "_distribution_groups\\.csv$", full.names = TRUE, recursive = T)

combined_df <- NULL

for (file in all_files) {
  df <- read.csv(file)
  
  # Assume first column is group, last column is normality
  lipid_name <- sub("_distribution_group\\.csv$", "", basename(file))
  df_lipid <- df[, c(1, ncol(df))]
  colnames(df_lipid) <- c("Group", lipid_name)  # rename normality column to lipid name
  
  if (is.null(combined_df)) {
    combined_df <- df_lipid
  } else {
    combined_df <- merge(combined_df, df_lipid, by = "Group", all = TRUE)
  }
}

write.csv(combined_df, file = file.path(top_level_dir,"..","total_lipids", "combined_per_group_normality.csv"), row.names = FALSE)

normality_per_group <- read.csv("outputs/total_lipids/combined_per_group_normality.csv", row.names = 1)
all_values <- unlist(normality_per_group)
table(all_values)
total <- length(all_values)
counts <- table(all_values)
percentages <- round(100 * counts / total, 1)
write.csv(percentages,file = file.path(top_level_dir,"..","total_lipids", "combined_per_group_normality_percentages.csv"), row.names = FALSE)

for (cat in names(counts)) {
  cat(cat, ":", counts[cat], "(", percentages[cat], "%)\n")
}

################################################################################
##### 16. LOG TRANSFORMATION AND NORMALITY OUTPUTS FOR ALL SAMPLES COMBINED ####
################################################################################

## For all lipids combined #####################################################

top_level_dir <- file.path("outputs", "lipid_categories")

for (name in raw_data_names) {
  lipid_family <- sub("^raw_data_", "", name)
  category <- category_mapping$Category_clean[category_mapping$Family_clean == lipid_family][1]
  
  if (is.na(category) || length(category) == 0) {
    message("No category found for family: ", lipid_family)
    next
  }
  
  folder_path <- file.path(top_level_dir, category, lipid_family)  
  df <- get(name)
  lipid_columns <- names(df)[-1]
  
  distribution_summary <- data.frame(
    lipid = lipid_columns,
    skewness = numeric(length(lipid_columns)),
    kurtosis = numeric(length(lipid_columns)),
    shapiro_p = numeric(length(lipid_columns)),
    normality = character(length(lipid_columns)),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(lipid_columns)) {
    lipid <- lipid_columns[i]
    
    # log-transform the values
    values <- log1p(df[[lipid]])  
    
    distribution_summary$skewness[i] <- skewness(values)
    distribution_summary$kurtosis[i] <- kurtosis(values)
    
    if (length(values) >= 3 & length(values) <= 5000) {
      shapiro_res <- shapiro.test(values)
      distribution_summary$shapiro_p[i] <- shapiro_res$p.value
      distribution_summary$normality[i] <- ifelse(shapiro_res$p.value > 0.05, "Normal", "Non-normal")
    } else {
      distribution_summary$shapiro_p[i] <- NA
      distribution_summary$normality[i] <- "NA"
    }
  }
  
  distribution_summary$shapiro_p_adj <- p.adjust(distribution_summary$shapiro_p, method = adjustment_method)
  distribution_summary$normality <- ifelse(distribution_summary$shapiro_p_adj > 0.05, "Normal", "Non-normal")
  distribution_summary <- distribution_summary[, c("lipid", "skewness", "kurtosis", "shapiro_p", "shapiro_p_adj", "normality")]
  
  cat("Saving log-transformed distribution summary for:", name, "\n")
  
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE, showWarnings = FALSE)
  }
  
  dis_csv_file <- file.path(folder_path, paste0(lipid_family, "_distribution_log.csv"))
  write.csv(distribution_summary, file = dis_csv_file, row.names = FALSE)
}

################################################################################
###### 17. LOG TRANSFORMATION AND NORMALITY OUTPUTS FOR INDIVIDUAL GROUPS ######
################################################################################

top_level_dir <- file.path("outputs", "lipid_categories")

for (name in raw_data_names) {
  lipid_family <- sub("^raw_data_", "", name)
  category <- category_mapping$Category_clean[category_mapping$Family_clean == lipid_family][1]
  
  if (is.na(category) || length(category) == 0) {
    message("No category found for family: ", lipid_family)
    next
  }
  
  folder_path <- file.path(top_level_dir, category, lipid_family)  
  df <- get(name)
  lipid_columns <- names(df)[-1]
  group_col <- names(df)[1]
  groups <- unique(df[[group_col]])
  
  per_group_results <- list()
  
  for (i in seq_along(lipid_columns)) {
    lipid <- lipid_columns[i]
    
    group_p <- c()
    group_normality <- c()
    
    for (g in groups) {
      vals_g <- df[df[[group_col]] == g, lipid]
      
      # log-transform values (safe with zeros)
      vals_g <- log1p(vals_g)  # replace with log(vals_g) if all > 0
      
      if (length(vals_g) >= 3 & length(vals_g) <= 5000) {
        if (length(unique(vals_g)) > 1) {
          p <- shapiro.test(vals_g)$p.value
          group_p <- c(group_p, p)
          group_normality <- c(group_normality, ifelse(p > 0.05, "Normal", "Non-normal"))
        } else {
          group_p <- c(group_p, NA)
          group_normality <- c(group_normality, "All identical")
        }
      } else {
        group_p <- c(group_p, NA)
        group_normality <- c(group_normality, NA)
      }
    }
    
    # adjust all p-values for this lipid (not inside the loop)
    group_p_adj <- p.adjust(group_p, method = adjustment_method)
    group_normality <- ifelse(group_p_adj > 0.05, "Normal", "Non-normal")
    
    per_group_results[[lipid]] <- data.frame(
      group = groups, 
      p_value = group_p, 
      p_value_adj = group_p_adj, 
      normality = group_normality
    )
  }
  
  # save all per-group results
  for (lipid in names(per_group_results)) {
    write.csv(
      per_group_results[[lipid]], 
      file = file.path(folder_path, paste0(lipid, "_log_distribution_groups.csv")), 
      row.names = FALSE
    )
    message("Saving *log-transformed* distribution per group for: ", lipid_family, " -> ", lipid)
  }
}

## combining normality of all lipids across groups #############################

top_level_dir <- file.path("outputs", "lipid_categories")

all_files <- list.files(top_level_dir, pattern = "_log_distribution_groups\\.csv$", full.names = TRUE, recursive = T)

combined_df <- NULL

for (file in all_files) {
  df <- read.csv(file)
  
  # Assume first column is group, last column is normality
  lipid_name <- sub("_log_distribution_groups\\.csv$", "", basename(file))
  df_lipid <- df[, c(1, ncol(df))]
  colnames(df_lipid) <- c("Group", lipid_name)  # rename normality column to lipid name
  
  if (is.null(combined_df)) {
    combined_df <- df_lipid
  } else {
    combined_df <- merge(combined_df, df_lipid, by = "Group", all = TRUE)
  }
}

# save combined normality results
write.csv(combined_df, file = file.path(top_level_dir,"..","total_lipids", "combined_log_per_group_normality.csv"), row.names = FALSE)

normality_per_group <- read.csv("outputs/total_lipids/combined_log_per_group_normality.csv", row.names = 1)
all_values <- unlist(normality_per_group)
table(all_values)
total <- length(all_values)
counts <- table(all_values)
percentages <- round(100 * counts / total, 1)
write.csv(percentages,file = file.path(top_level_dir,"..","total_lipids", "combined_log_per_group_normality_percentages.csv"), row.names = FALSE)

for (cat in names(counts)) {
  cat(cat, ":", counts[cat], "(", percentages[cat], "%)\n")
}

################################################################################
############### 18. MANN-WHITNEY U TEST - TWO INDEPENDENT GROUPS ###############
################################################################################

top_level_dir <- file.path("outputs", "lipid_categories")

count <- 1  # counter for families

for (name in raw_data_names) {
  lipid_family <- sub("^raw_data_", "", name)
  category <- category_mapping$Category_clean[category_mapping$Family_clean == lipid_family][1]
  
  if (is.na(category) || length(category) == 0) {
    message(count, ". Skipping (no category): ", lipid_family)
    count <- count + 1
    next
  }
  
  folder_path <- file.path(top_level_dir, category, lipid_family)
  if (!dir.exists(folder_path)) dir.create(folder_path, recursive = TRUE)
  
  df <- get(name)
  
  lipid_columns <- names(df)[-1]   # exclude 'groups'
  group_col <- names(df)[1]
  
  mw_results <- data.frame(
    lipid = lipid_columns,
    p_value = NA,
    p_value_adj = NA,
    significance = NA,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(lipid_columns)) {
    lipid <- lipid_columns[i]
    control_vals <- df[df[[group_col]] == control_group, lipid]
    other_groups <- setdiff(unique(df[[group_col]]), control_group)
    
    pvals <- c()
    for (g in other_groups) {
      test_vals <- df[df[[group_col]] == g, lipid]
      if (length(control_vals) > 0 & length(test_vals) > 0) {
        test_res <- wilcox.test(control_vals, test_vals)
        pvals <- c(pvals, test_res$p.value)
      } else {
        pvals <- c(pvals, NA)
      }
    }
    
    if (all(is.na(pvals))) {
      mw_results$p_value[i] <- NA
    } else {
      mw_results$p_value[i] <- min(pvals, na.rm = TRUE)
    }
  }
  
  # adjust p-values
  mw_results$p_value_adj <- p.adjust(mw_results$p_value, method = adjustment_method)
  
  # assign significance stars
  mw_results$significance <- sapply(mw_results$p_value_adj, function(p) {
    if (is.na(p)) {
      "not significant"
    } else if (p < 0.001) {
      "***"
    } else if (p < 0.01) {
      "**"
    } else if (p < 0.05) {
      "*"
    } else {
      "not significant"
    }
  })
  
  write.csv(
    mw_results, 
    file = file.path(folder_path, paste0(lipid_family, "_mannwhitney.csv")), 
    row.names = FALSE
  )
  
  message(count, ". Mann–Whitney results saved for family: ", lipid_family)
  count <- count + 1
}

################################################################################
############### 19. KRUSKAL-WALLIS H TEST - >2 INDEPENDENT GROUPS ##############
################################################################################

top_level_dir <- file.path("outputs", "lipid_categories")

count <- 1  # counter for families

for (name in raw_data_names) {
  lipid_family <- sub("^raw_data_", "", name)
  category <- category_mapping$Category_clean[category_mapping$Family_clean == lipid_family][1]
  
  if (is.na(category) || length(category) == 0) {
    message(count, ". Skipping (no category): ", lipid_family)
    count <- count + 1
    next
  }
  
  folder_path <- file.path(top_level_dir, category, lipid_family)
  if (!dir.exists(folder_path)) dir.create(folder_path, recursive = TRUE)
  
  df <- get(name)
  
  lipid_columns <- names(df)[-1]   # exclude 'groups' column
  group_col <- names(df)[1]
  
  kw_results <- data.frame(
    lipid = lipid_columns,
    p_value = NA,
    p_value_adj = NA,
    significance = NA,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(lipid_columns)) {
    lipid <- lipid_columns[i]
    
    # remove NA values
    df_subset <- df[!is.na(df[[lipid]]), c(group_col, lipid)]
    
    # run Kruskal-Wallis if at least 2 groups exist
    if (length(unique(df_subset[[group_col]])) > 1) {
      test_res <- kruskal.test(df_subset[[lipid]] ~ df_subset[[group_col]])
      kw_results$p_value[i] <- test_res$p.value
    }
  }
  
  # adjust p-values
  kw_results$p_value_adj <- p.adjust(kw_results$p_value, method = adjustment_method)
  
  # assign significance stars
  kw_results$significance <- sapply(kw_results$p_value_adj, function(p) {
    if (is.na(p)) {
      "not significant"
    } else if (p < 0.001) {
      "***"
    } else if (p < 0.01) {
      "**"
    } else if (p < 0.05) {
      "*"
    } else {
      "not significant"
    }
  })
  
  write.csv(
    kw_results, 
    file = file.path(folder_path, paste0(lipid_family, "_kruskalwallis.csv")), 
    row.names = FALSE
  )
  
  message(count, ". Kruskal–Wallis results saved for family: ", lipid_family)
  count <- count + 1
}

################################################################################
######################## 12. DUNN'S TEST - POST HOC ############################
################################################################################

top_level_dir <- file.path("outputs", "lipid_categories")
count <- 1  # counter for families

for (name in raw_data_names) {
  lipid_family <- sub("^raw_data_", "", name)
  category <- category_mapping$Category_clean[category_mapping$Family_clean == lipid_family][1]
  
  if (is.na(category) || length(category) == 0) {
    message(count, ". Skipping (no category): ", lipid_family)
    count <- count + 1
    next
  }
  
  folder_path <- file.path(top_level_dir, category, lipid_family)
  if (!dir.exists(folder_path)) dir.create(folder_path, recursive = TRUE)
  
  df <- get(name)
  lipid_columns <- names(df)[-1]   # exclude 'groups' column
  group_col <- names(df)[1]
  
  for (i in seq_along(lipid_columns)) {
    lipid <- lipid_columns[i]
    df_subset <- df[!is.na(df[[lipid]]), c(group_col, lipid)]
    
    # Skip lipids with all identical values
    if (length(unique(df_subset[[lipid]])) <= 1) {
      message("Skipping lipid ", lipid, " (all values identical)")
      next
    }
    
    # Count non-NA observations per group
    group_counts <- table(df_subset[[group_col]])
    
    # Only run Dunn if at least 2 groups have ≥1 value
    if (sum(group_counts > 0) >= 2) {
      # Run Dunn quietly
      dunn_res <- suppressMessages(
        dunn.test(df_subset[[lipid]], df_subset[[group_col]], method = "bonferroni", alpha = 0.05, altp = TRUE)
      )
      
      # Only create results if valid
      if (!is.null(dunn_res$comparisons) && length(dunn_res$comparisons) > 0 &&
          !is.null(dunn_res$P) && length(dunn_res$P) == length(dunn_res$comparisons) &&
          !is.null(dunn_res$P.adjusted) && length(dunn_res$P.adjusted) == length(dunn_res$comparisons)) {
        
        results_df <- data.frame(
          comparison = dunn_res$comparisons,
          p_value = dunn_res$P,
          p_value_adj = dunn_res$P.adjusted,
          significance = sapply(dunn_res$P.adjusted, function(p) {
            if (is.na(p)) {
              "not significant"
            } else if (p < 0.001) {
              "***"
            } else if (p < 0.01) {
              "**"
            } else if (p < 0.05) {
              "*"
            } else {
              "not significant"
            }
          }),
          stringsAsFactors = FALSE
        )
        
        write.csv(
          results_df,
          file = file.path(folder_path, paste0(lipid_family, "_", lipid, "_dunn.csv")),
          row.names = FALSE
        )
      } else {
        message("Skipping lipid ", lipid, " (Dunn test returned no valid comparisons)")
      }
      
    } else {
      message("Skipping lipid ", lipid, " (less than 2 groups with data)")
    }
  }
  
  message(count, ". Dunn's post hoc results saved for family: ", lipid_family)
  count <- count + 1
}

################################################################################
#################### 20. SPEARMAN'S RHO - CORRELATION ##########################
################################################################################



################################################################################
##################### 21. KENDALL'S TAU - CORRELATION ##########################
################################################################################



################################################################################
####################### 22. GENERALISED LINEAR MODEL ###########################
################################################################################

















































