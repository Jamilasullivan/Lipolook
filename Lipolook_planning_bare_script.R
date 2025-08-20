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

groups <- ## what is the current name of your column containing groupings in raw data?
  
control_group <- ## what is the name of your control group in raw data?

################################################################################
############################ PACKAGES TO INSTALL ###############################
################################################################################

##install.packages("moments")
library(moments)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stats)

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
  
  # Build path: outputs/lipid_categories/category/family
  folder_path <- file.path(top_level_dir, category, lipid_family)
  
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE, showWarnings = FALSE)
    message(count, ". Folder created: ", folder_path)
    count <- count + 1
  }
}

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
  
  folder_path <- file.path(top_level_dir, category, lipid_family)
  
  cat("Saving averages for:", name, "to folder:", folder_path, "\n")
  
  df_avg <- aggregate(. ~ groups, data = get(name), FUN = mean)
  csv_file <- file.path(folder_path, paste0(lipid_family, "_avg.csv"))
  write.csv(df_avg, file = csv_file, row.names = FALSE)
}

################################################################################
######################### HISTOGRAMS FOR NORMALITY #############################
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
################# NORMALITY OUTPUT FOR ALL SAMPLES COMBINED ####################
################################################################################

## across all samples for each lipid

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
    values <- df[[lipid]]  # across all samples
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
################### NORMALITY OUTPUT FOR INDIVIDUAL GROUPS #####################
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
  lipid_columns <- names(df)[-1]  # assuming first column is group/sample
  group_col <- names(df)[1]
  groups <- unique(df[[group_col]])
  
  distribution_summary <- data.frame(
    lipid = lipid_columns,
    skewness_overall = numeric(length(lipid_columns)),
    kurtosis_overall = numeric(length(lipid_columns)),
    shapiro_p_overall = numeric(length(lipid_columns)),
    shapiro_p_adj_overall = numeric(length(lipid_columns)),
    normality_overall = character(length(lipid_columns)),
    stringsAsFactors = FALSE
  )
  
  per_group_results <- list()
  
  for (i in seq_along(lipid_columns)) {
    lipid <- lipid_columns[i]
    values <- df[[lipid]]
    
    # Overall stats
    distribution_summary$skewness_overall[i] <- skewness(values)
    distribution_summary$kurtosis_overall[i] <- kurtosis(values)
    
    if (length(values) >= 3 & length(values) <= 5000 && length(unique(values)) > 1) {
      distribution_summary$shapiro_p_overall[i] <- shapiro.test(values)$p.value
    } else {
      distribution_summary$shapiro_p_overall[i] <- NA
    }
    
    # Per-group normality
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
    per_group_results[[lipid]] <- data.frame(group = groups, p_value = group_p, normality = group_normality)
  }
  
  # Adjust overall p-values and decide normality
  distribution_summary$shapiro_p_adj_overall <- p.adjust(distribution_summary$shapiro_p_overall, method = adjustment_method)
  distribution_summary$normality_overall <- ifelse(distribution_summary$shapiro_p_adj_overall > 0.05, "Normal", "Non-normal")
  
  # Reorder columns
  distribution_summary <- distribution_summary[, c("lipid", "skewness_overall", "kurtosis_overall",
                                                   "shapiro_p_overall", "shapiro_p_adj_overall",
                                                   "normality_overall")]
  
  # Save overall summary
  if (!dir.exists(folder_path)) dir.create(folder_path, recursive = TRUE, showWarnings = FALSE)
  write.csv(distribution_summary, file = file.path(folder_path, paste0(lipid_family, "_distribution.csv")), row.names = FALSE)
  message("Saving")
  
  
  
  # Save per-group summaries per lipid
  for (lipid in names(per_group_results)) {
    write.csv(per_group_results[[lipid]], file = file.path(folder_path, paste0(lipid, "distribution_group.csv")), row.names = FALSE)
    message("Saving distribution per group for: ", lipid_family)
  }
}

## combining normality of all lipids across groups

top_level_dir <- file.path("outputs", "lipid_categories")

all_files <- list.files(top_level_dir, pattern = ".distribution_group\\.csv$", full.names = TRUE, recursive = T)

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

# Save combined CSV
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
############################ LOG TRANSFORMATION ################################
################################################################################










































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
