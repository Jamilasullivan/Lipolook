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

groups <- "Group" ## what is the current name of your column containing groupings in raw data?
  
control_group <- "E" ## the name of your control group in your raw data file. what is your control group called?
  
glm_variables <- c("","","") ## include the names of columns with metadata to be included

################################################################################
########################## 2. PACKAGES TO INSTALL ##############################
################################################################################

# List all required packages
packages <- c(
  "moments", "tidyr", "dplyr", "ggplot2", "stats",
  "dunn.test", "rmarkdown", "corrplot", "pheatmap"
)

# Install any that are missing
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Load them all
lapply(packages, library, character.only = TRUE)

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

raw_data_lipids <- raw_data
duplicated_columns <- duplicated(as.list(raw_data_lipids))
any(duplicated_columns)
duplicated_column_names <- names(raw_data_lipids)[duplicated_columns]
print(duplicated_column_names)
duplicated_column_names <- data.frame(duplicated_columns = duplicated_column_names)
write.csv(duplicated_column_names, "outputs/error_files/duplicated_columns.csv", row.names = FALSE)

metadata <- raw_data[ , !(grepl("\\.$", names(raw_data)) | names(raw_data) == "Cholesterol")]
groups <- metadata[,"Group"]
raw_data_lipids <- raw_data[ , grepl("\\.$", names(raw_data)) | names(raw_data) == "Cholesterol"]
str(raw_data_lipids)
raw_data_lipids[] <- lapply(raw_data_lipids, function(x) as.numeric(as.character(x)))
str(raw_data_lipids)
sum(is.na(raw_data_lipids[,ncol(raw_data_lipids)]))

################################################################################
####################### 6. LIPIDS TESTED MANIPULATION ##########################
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
###################### 7. CHECKING FOR MISMATCHED COLUMNS ######################
######################## I.E. COLUMNS TESTED THAT THERE ########################
########################## WAS NO RESULTANT DATA FOR ###########################
################################################################################

raw_data_columns <- colnames(raw_data_lipids[ncol(raw_data_lipids)])
lipids <- lipids_tested$lipid
unmatched_cols <- setdiff(raw_data_columns, lipids)
number_unmatched <- length(unmatched_cols)
unmatched_cols <- data.frame(Unmatched_columns = unmatched_cols)
write.csv(unmatched_cols, "outputs/error_files/unmatched_columns.csv", row.names = FALSE) # this should be empty because the mismatched columns now end up in the metadata object

################################################################################
####################### 8. SUBSET DATA BY LIPID FAMILY #########################
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
raw_data_names <- raw_data_names[!raw_data_names %in% c("raw_data_Ceramide_phospho_ethanolamines", "raw_data_Cholesterol_cholesteryl_esters", "raw_data_Lysophosphatidy_ethanolamines")] # this section of code is redundant for the results and this should be considered going forwards. The input should be standardized before being given to the client to save this problem in the future.

raw_data_names <- raw_data_names[raw_data_names != "raw_data_lipids"]
raw_data_names

################################################################################
######### 9. CREATING FOLDERS FOR ALL LIPID FAMILIES AND CATEGORIES ############
################################################################################

write.csv(
  data.frame(Family = names(lipid_families)),
  "lipid_categories_1.csv",
  row.names = FALSE
) ## THIS NEEDS TO BE FILLED IN BY THE CLIENT AND SAVED AS 'complete_lipid_categories.csv'

category_mapping <- read.csv("complete_lipid_categories.csv", stringsAsFactors = FALSE)

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
###################### 10. ADDING GROUPS BACK TO DFS ###########################
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
###################### 11. CALCULATING GROUP AVERGAGES #########################
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
####################### 12. HISTOGRAMS FOR NORMALITY ###########################
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
################# 13. FOREST PLOTS - GROUPS WITHIN A FAMILIY ###################
################################################################################

#### SAVING TOTALS FOR EACH FAMILY #############################################

counter <- 1
total_families <- length(raw_data_names)

for (name in raw_data_names) {
  # Display progress
  lipid_family <- sub("^raw_data_", "", name)
  message("Processing family ", counter, " of ", total_families, ": ", lipid_family)
  
  # Get the data frame
  df <- get(name)
  group_col <- names(df)[1]
  lipid_columns <- names(df)[-1]
  
  # Convert lipid columns to numeric safely
  df[, lipid_columns] <- lapply(df[, lipid_columns, drop = FALSE], function(x) as.numeric(as.character(x)))
  
  # Compute total per sample
  total_df <- data.frame(
    group = df[[group_col]],
    total = rowSums(df[, lipid_columns, drop = FALSE], na.rm = TRUE)
  )
  
  # Assign the new data frame in the environment
  total_name <- paste0("total_", lipid_family)
  assign(total_name, total_df)
  
  # Increment counter
  counter <- counter + 1
}

##### FOREST PLOTS (ABSOLUTE DIFFERENCE) #######################################

top_level_dir <- file.path("outputs", "lipid_categories")
counter <- 1
total_families <- length(raw_data_names)

for (name in raw_data_names) {
  lipid_family <- sub("^raw_data_", "", name)
  message("Processing family ", counter, " of ", total_families, ": ", lipid_family)
  category <- category_mapping$Category_clean[category_mapping$Family_clean == lipid_family][1]
  
  if (is.na(category) || length(category) == 0) {
    message("No category found for family: ", lipid_family)
    counter <- counter + 1
    next
  }
  
  folder_path <- file.path(top_level_dir, category, lipid_family)
  if (!dir.exists(folder_path)) dir.create(folder_path, recursive = TRUE)
  
  total_name <- paste0("total_", lipid_family)
  if (!exists(total_name)) {
    message("No total data frame found for: ", lipid_family)
    counter <- counter + 1
    next
  }
  
  df_total <- get(total_name)  # must have columns: group, total
  group_col <- "group"
  
  df_total$total <- as.numeric(df_total$total)
  
  summary_df <- df_total %>%
    group_by(.data[[group_col]]) %>%
    summarize(
      mean_total = mean(total, na.rm = TRUE),
      sd_total = sd(total, na.rm = TRUE),
      n = n(),
      se_total = sd_total / sqrt(n),
      .groups = "drop"
    )
  
  if (!control_group %in% summary_df[[group_col]]) {
    message("Control group not found for ", lipid_family, " — skipping plot")
    counter <- counter + 1
    next
  }
  
  control_mean <- summary_df$mean_total[summary_df[[group_col]] == control_group]
  summary_df <- summary_df %>%
    mutate(diff_from_control = mean_total - control_mean)

  plot_df <- summary_df %>% filter(.data[[group_col]] != control_group)

  max_diff <- max(abs(plot_df$diff_from_control + plot_df$se_total), 
                  abs(plot_df$diff_from_control - plot_df$se_total), na.rm = TRUE)

  plot <- ggplot(plot_df, aes(y = .data[[group_col]], x = diff_from_control)) +
    geom_point(size = 3, color = "#990101") +
    geom_errorbarh(aes(
      xmin = diff_from_control - se_total,
      xmax = diff_from_control + se_total
    ), height = 0.1, color = "#990101") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#747373") +
    labs(
      title = paste("Forest plot (relative to Control):", lipid_family),
      x = "Difference from Control (sum of lipids)",
      y = "Group"
    ) +
    theme_bw() +
    coord_cartesian(xlim = c(-max_diff, max_diff))
  
  ggsave(
    filename = file.path(folder_path, paste0(lipid_family, "_forest_plot.png")),
    plot = plot,
    width = 6,
    height = 8
  )
  
  message("Saved forest plot for: ", lipid_family)
  counter <- counter + 1
}

#### FOREST PLOTS (LOG2 FOLD CHANGE) ###########################################

top_level_dir <- file.path("outputs", "lipid_categories")
counter <- 1
total_families <- length(raw_data_names)

for (name in raw_data_names) {
  lipid_family <- sub("^raw_data_", "", name)
  message("Processing family ", counter, " of ", total_families, ": ", lipid_family)
  
  category <- category_mapping$Category_clean[category_mapping$Family_clean == lipid_family][1]
  if (is.na(category) || length(category) == 0) {
    message("No category found for family: ", lipid_family)
    counter <- counter + 1
    next
  }
  
  folder_path <- file.path(top_level_dir, category, lipid_family)
  if (!dir.exists(folder_path)) dir.create(folder_path, recursive = TRUE)
  
  total_name <- paste0("total_", lipid_family)
  if (!exists(total_name)) {
    message("No total data frame found for: ", lipid_family)
    counter <- counter + 1
    next
  }
  
  df_total <- get(total_name)  # must have columns: group, total
  df_total$total <- as.numeric(df_total$total)
  group_col <- "group"
  
  # Summary per group
  summary_df <- df_total %>%
    group_by(.data[[group_col]]) %>%
    summarize(
      mean_total = mean(total, na.rm = TRUE),
      sd_total = sd(total, na.rm = TRUE),
      n = n(),
      se_total = sd_total / sqrt(n),
      .groups = "drop"
    )
  
  # Skip if control group not present
  if (!control_group %in% summary_df[[group_col]]) {
    message("Control group not found for ", lipid_family, " — skipping plot")
    counter <- counter + 1
    next
  }
  
  # Compute log2 fold change relative to control
  control_mean <- summary_df$mean_total[summary_df[[group_col]] == control_group]
  summary_df <- summary_df %>%
    mutate(
      fold_change = mean_total / control_mean,
      log2_fold_change = log2(fold_change),
      log2_se = se_total / (control_mean * log(2))  # approximate SE on log2 scale
    )
  
  # Exclude control group from plotting
  plot_df <- summary_df %>% filter(.data[[group_col]] != control_group)
  
  # Determine symmetric x-axis limits
  max_abs_log2 <- max(abs(plot_df$log2_fold_change + plot_df$log2_se),
                      abs(plot_df$log2_fold_change - plot_df$log2_se),
                      na.rm = TRUE)
  
  # Forest plot
  plot <- ggplot(plot_df, aes(y = .data[[group_col]], x = log2_fold_change)) +
    geom_point(size = 3, color = "#05016F") +
    geom_errorbarh(aes(
      xmin = log2_fold_change - log2_se,
      xmax = log2_fold_change + log2_se
    ), height = 0.1, color = "#05016F") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#747373") +
    labs(
      title = paste("Forest plot (log2 fold change vs Control):", lipid_family),
      x = "Log2 Fold Change from Control",
      y = "Group"
    ) +
    theme_bw() +
    coord_cartesian(xlim = c(-max_abs_log2, max_abs_log2))
  
  # Save plot
  ggsave(
    filename = file.path(folder_path, paste0(lipid_family, "_log2_fold_change_forest_plot.png")),
    plot = plot,
    width = 6,
    height = 8
  )
  
  message("Saved log2 fold change forest plot for: ", lipid_family)
  counter <- counter + 1
}

################################################################################
################ 14. FOREST PLOTS - GROUPS WITHIN A CATEGORY ###################
################################################################################

#### SAVING TOTALS FOR EACH CATEGORY ###########################################

top_level_dir <- file.path("outputs", "lipid_categories")
unique_categories <- unique(category_mapping$Category_clean)
counter <- 1
total_categories <- length(unique_categories)

category_dfs <- list()  

for (category in unique_categories) {
  message("Processing category ", counter, " of ", total_categories, ": ", category)
  families_in_cat <- category_mapping$Family_clean[category_mapping$Category_clean == category]
  
  family_dfs <- list()
  for (family in families_in_cat) {
    total_name <- paste0("total_", family)
    if (exists(total_name)) {
      family_dfs[[family]] <- get(total_name)
    }
  }
  
  if (length(family_dfs) == 0) {
    message("No family data found for category: ", category)
    counter <- counter + 1
    next
  }
  
  totals_only <- lapply(family_dfs, function(df) df$total)
  
  if (length(totals_only) == 1) {
    category_df <- data.frame(
      group = family_dfs[[1]]$group
    )
    category_df[[paste0("total_", names(family_dfs))]] <- totals_only[[1]]
    
    category_df$total <- totals_only[[1]]
    
  } else {
    category_df <- data.frame(
      group = family_dfs[[1]]$group,
      do.call(cbind, totals_only)
    )
    names(category_df)[-1] <- paste0("total_", names(family_dfs))
    category_df$total <- rowSums(category_df[,-1], na.rm = TRUE)
  }
  
  assign(paste0("category_", category), category_df)
  category_dfs[[category]] <- category_df
  
  message("Category data frame created with ", nrow(category_df), " rows and ", ncol(category_df), " columns (including total column)")
  
  counter <- counter + 1
}

#### FOREST PLOTS (ABSOLUTE DIFFERENCE FROM MEAN) ##############################

top_level_dir <- file.path("outputs", "lipid_categories")
counter <- 1
total_categories <- length(category_dfs)

for (category in names(category_dfs)) {
  message("Processing category ", counter, " of ", total_categories, ": ", category)
  
  df_category <- category_dfs[[category]]
  
  df_category$total <- as.numeric(df_category$total)
  
  summary_df <- df_category %>%
    group_by(group) %>%
    summarize(
      mean_total = mean(total, na.rm = TRUE),
      sd_total = sd(total, na.rm = TRUE),
      n = n(),
      se_total = sd_total / sqrt(n),
      .groups = "drop"
    )
  
  if (!control_group %in% summary_df$group) {
    message("Control group not found for category ", category, " — skipping plot")
    counter <- counter + 1
    next
  }
  
  control_mean <- summary_df$mean_total[summary_df$group == control_group]
  summary_df <- summary_df %>%
    mutate(diff_from_control = mean_total - control_mean)

  plot_df <- summary_df %>% filter(group != control_group)
  
  if (nrow(plot_df) == 0) {
    message("No groups to plot for category ", category, " — skipping plot")
    counter <- counter + 1
    next
  }

  max_diff <- max(abs(plot_df$diff_from_control + plot_df$se_total), 
                  abs(plot_df$diff_from_control - plot_df$se_total), na.rm = TRUE)

  plot <- ggplot(plot_df, aes(y = group, x = diff_from_control)) +
    geom_point(size = 3, color = "#990101") +
    geom_errorbarh(aes(
      xmin = diff_from_control - se_total,
      xmax = diff_from_control + se_total
    ), height = 0.1, color = "#990101") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#747373") +
    labs(
      title = paste("Forest plot (relative to Control):", category),
      x = "Difference from Control (sum of lipids)",
      y = "Group"
    ) +
    theme_bw() +
    coord_cartesian(xlim = c(-max_diff, max_diff))

  folder_path <- file.path(top_level_dir, category)
  if (!dir.exists(folder_path)) dir.create(folder_path, recursive = TRUE)

  ggsave(
    filename = file.path(folder_path, paste0(category, "_absolute_change_forest_plot.png")),
    plot = plot,
    width = 6,
    height = 8
  )
  
  message("Saved forest plot for category: ", category)
  counter <- counter + 1
}

#### FOREST PLOTS (LOG2 FOLD CHANGE DIFFERENCE FROM MEAN) ######################

top_level_dir <- file.path("outputs", "lipid_categories")
counter <- 1
total_categories <- length(category_dfs)  # category_dfs from previous step

for (category in names(category_dfs)) {
  message("Processing category ", counter, " of ", total_categories, ": ", category)
  
  # Get category data frame
  df_category <- category_dfs[[category]]  # must have 'group' and 'total'
  
  # Ensure numeric
  df_category$total <- as.numeric(df_category$total)
  
  # Compute summary per group
  summary_df <- df_category %>%
    group_by(group) %>%
    summarize(
      mean_total = mean(total, na.rm = TRUE),
      sd_total = sd(total, na.rm = TRUE),
      n = n(),
      se_total = sd_total / sqrt(n),
      .groups = "drop"
    )
  
  # Skip if control group is missing
  if (!control_group %in% summary_df$group) {
    message("Control group not found for category ", category, " — skipping plot")
    counter <- counter + 1
    next
  }
  
  # Compute fold change relative to control
  control_mean <- summary_df$mean_total[summary_df$group == control_group]
  summary_df <- summary_df %>%
    mutate(
      fold_change = mean_total / control_mean,
      log2_fold_change = log2(fold_change)
    )
  
  # Exclude control group from plotting
  plot_df <- summary_df %>% filter(group != control_group)
  
  if (nrow(plot_df) == 0) {
    message("No groups to plot for category ", category, " — skipping plot")
    counter <- counter + 1
    next
  }
  
  # Determine symmetric x-axis limits using SEs from plot_df only
  max_abs_log2 <- max(
    abs(plot_df$log2_fold_change + (plot_df$se_total / (control_mean * log(2)))),
    abs(plot_df$log2_fold_change - (plot_df$se_total / (control_mean * log(2)))),
    na.rm = TRUE
  )
  
  # Forest plot using log2 fold change
  plot <- ggplot(plot_df, aes(y = group, x = log2_fold_change)) +
    geom_point(size = 3, color = "#990101") +
    geom_errorbarh(aes(
      xmin = log2_fold_change - (se_total / (control_mean * log(2))),
      xmax = log2_fold_change + (se_total / (control_mean * log(2)))
    ), height = 0.1, color = "#990101") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#747373") +
    labs(
      title = paste("Category-level Forest Plot (log2 fold change):", category),
      x = "Log2 Fold Change from Control",
      y = "Group"
    ) +
    theme_bw() +
    coord_cartesian(xlim = c(-max_abs_log2, max_abs_log2))
  
  # Create category folder if it doesn't exist
  folder_path <- file.path(top_level_dir, category)
  if (!dir.exists(folder_path)) dir.create(folder_path, recursive = TRUE)
  
  # Save the forest plot
  ggsave(
    filename = file.path(folder_path, paste0(category, "_log2_change_forest_plot.png")),
    plot = plot,
    width = 6,
    height = 8
  )
  
  message("Saved forest plot for category: ", category)
  counter <- counter + 1
}

################################################################################
############### 15. NORMALITY OUTPUT FOR ALL SAMPLES COMBINED ##################
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
################# 16. NORMALITY OUTPUT FOR INDIVIDUAL GROUPS ###################
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
    
    group_p_adj <- p.adjust(group_p, method = adjustment_method)
    group_normality <- ifelse(group_p_adj > 0.05, "Normal", "Non-normal")
    
    per_group_results[[lipid]] <- data.frame(
      group = groups, 
      p_value = group_p, 
      p_value_adj = group_p_adj, 
      normality = group_normality
    )
  }
  
  for (lipid in names(per_group_results)) {
    write.csv(
      per_group_results[[lipid]], 
      file = file.path(folder_path, paste0(lipid, "_distribution_groups_raw.csv")), 
      row.names = FALSE
    )
    message("Saving RAW distribution per group for: ", lipid_family, " -> ", lipid)
  }
}

#### conbining values ##########################################################

all_files <- list.files(
  top_level_dir, 
  pattern = "_distribution_groups_raw\\.csv$", 
  full.names = TRUE, 
  recursive = TRUE
)

combined_df <- NULL

for (file in all_files) {
  df <- read.csv(file)
  
  lipid_name <- sub("_distribution_groups_raw\\.csv$", "", basename(file))
  df_lipid <- df[, c(1, ncol(df))]
  colnames(df_lipid) <- c("Group", lipid_name) 
  
  if (is.null(combined_df)) {
    combined_df <- df_lipid
  } else {
    combined_df <- merge(combined_df, df_lipid, by = "Group", all = TRUE)
  }
}

write.csv(
  combined_df, 
  file = file.path(top_level_dir,"..","total_lipids", "combined_per_group_normality_raw.csv"), 
  row.names = FALSE
)

normality_per_group <- read.csv("outputs/total_lipids/combined_per_group_normality_raw.csv", row.names = 1)
all_values <- unlist(normality_per_group)
counts <- table(all_values)
total <- length(all_values)
percentages <- round(100 * counts / total, 1)

write.csv(
  percentages,
  file = file.path(top_level_dir,"..","total_lipids", "combined_per_group_normality_raw_percentages.csv"),
  row.names = TRUE
)

for (cat in names(counts)) {
  cat(cat, ":", counts[cat], "(", percentages[cat], "%)\n")
}

################################################################################
##### 17. LOG TRANSFORMATION AND NORMALITY OUTPUTS FOR ALL SAMPLES COMBINED ####
################################################################################

## For all lipids combined #####################################################

top_level_dir <- file.path("outputs", "lipid_categories")

all_lipids_distribution <- list()

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
    values <- log1p(df[[lipid]])
    
    distribution_summary$skewness[i] <- skewness(values, na.rm = TRUE)
    distribution_summary$kurtosis[i] <- kurtosis(values, na.rm = TRUE)
    
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
  
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE, showWarnings = FALSE)
  }
  
  dis_csv_file <- file.path(folder_path, paste0(lipid_family, "_distribution_log.csv"))
  write.csv(distribution_summary, file = dis_csv_file, row.names = FALSE)
  cat("Saved distribution summary for:", lipid_family, "\n")

  distribution_summary$family <- lipid_family
  distribution_summary$category <- category
  all_lipids_distribution[[length(all_lipids_distribution) + 1]] <- distribution_summary
}

combined_summary <- bind_rows(all_lipids_distribution)

combined_file <- file.path("outputs", "total_lipids", "combined_lipid_normality_log.csv")
write.csv(combined_summary, file = combined_file, row.names = FALSE)

distribution_counts <- table(combined_summary$normality)
distribution_counts_df <- as.data.frame(distribution_counts)
distribution_total <- sum(distribution_counts)
distribution_counts_df$percent <- round(100 * distribution_counts_df$Freq / distribution_total, 1)

summary_file <- file.path("outputs", "total_lipids", "combined_lipid_normality_summary_log.csv")
write.csv(distribution_counts_df, file = summary_file, row.names = FALSE)

cat("Combined summary saved with", distribution_total, "lipids.\n")
print(distribution_counts_df)

################################################################################
###### 18. LOG TRANSFORMATION AND NORMALITY OUTPUTS FOR INDIVIDUAL GROUPS ######
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
############### 19. MANN-WHITNEY U TEST - TWO INDEPENDENT GROUPS ###############
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
############### 20. KRUSKAL-WALLIS H TEST (>2 INDEPENDENT GROUPS) ##############
########################## AND DUNN'S POST HOC TEST ############################
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
  
  dunn_results_all <- data.frame(
    lipid = character(),
    comparison = character(),
    p_value = numeric(),
    p_value_adj = numeric(),
    significance = character(),
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
      
      # Run Dunn's test if Kruskal-Wallis p < 0.05
      if (test_res$p.value < 0.05) {
        # Run Dunn's test
        dunn_res <- dunn.test(df_subset[[lipid]], g = df_subset[[group_col]], method = "none", kw = FALSE, altp = FALSE)
        
        if (length(dunn_res$P) > 0) {
          # manually create pairwise comparison names
          group_levels <- unique(df_subset[[group_col]])
          comparisons <- combn(group_levels, 2, FUN = function(x) paste(x[1], "vs", x[2]))
          
          pvals <- dunn_res$P
          pvals_adj <- p.adjust(pvals, method = adjustment_method)
          
          significance <- sapply(pvals_adj, function(p) {
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
          
          dunn_results_all <- rbind(dunn_results_all,
                                    data.frame(
                                      lipid = lipid,
                                      comparison = comparisons,
                                      p_value = pvals,
                                      p_value_adj = pvals_adj,
                                      significance = significance,
                                      stringsAsFactors = FALSE
                                    ))
        }
      }
    }
  }
  # write Kruskal-Wallis results
  write.csv(
    kw_results, 
    file = file.path(folder_path, paste0(lipid_family, "_kruskalwallis.csv")), 
    row.names = FALSE
  )
  
  # write Dunn test results if any
  if (nrow(dunn_results_all) > 0) {
    write.csv(
      dunn_results_all, 
      file = file.path(folder_path, paste0(lipid_family, "_dunn.csv")), 
      row.names = FALSE
    )
  }
  
  message(count, ". Kruskal–Wallis (+Dunn if significant) results saved for family: ", lipid_family)
  count <- count + 1
}

################################################################################
########################## SPEARMAN CORRELATION TEST ###########################
################################################################################

#### FOR ALL LIPIDS COMBINED ###################################################

cor_mat <- cor(raw_data_lipids, use = "pairwise.complete.obs", method = "spearman")

range(cor_mat, na.rm = TRUE)

output_folder <- "outputs/total_lipids"
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

csv_file <- file.path(output_folder, "lipid_spearman_correlation.csv")
write.csv(cor_mat, file = csv_file, row.names = TRUE)

message("Spearman correlation CSV saved to: ", csv_file)

#### FOR EACH LIPID CATEGORY ###################################################

top_level_dir <- file.path("outputs", "lipid_categories")
if (!dir.exists(top_level_dir)) dir.create(top_level_dir, recursive = TRUE)

count <- 1

for (category in unique(category_mapping$Category_clean)) {

  families_in_cat <- category_mapping$Family_clean[category_mapping$Category_clean == category]

  lipid_species <- unlist(lapply(families_in_cat, function(fam) {
    lipid_families[[fam]] 
  }))

  lipid_species <- lipid_species[lipid_species %in% colnames(cor_mat)]
  
  if (length(lipid_species) < 2) {
    message("Not enough lipids for category ", category, " — skipping")
    next
  }

  cor_sub <- cor_mat[lipid_species, lipid_species]
  
  folder_path <- file.path(top_level_dir, category)
  if (!dir.exists(folder_path)) dir.create(folder_path, recursive = TRUE)
  
  csv_file <- file.path(folder_path, paste0(category, "_spearman_correlation.csv"))
  write.csv(cor_sub, file = csv_file, row.names = TRUE)
  
  message(count, ". Saved Spearman correlation CSV for category: ", category)
  count <- count + 1
}

################################################################################
############################ EXTRA VISUALISATIONS ##############################
################################################################################

################################################################################
####################### AVERAGES BAR CHART FOR CATEGORIES ######################
################################################################################

#### GETTING AVERAGES PER CATEGORY #############################################

top_level_dir <- file.path("outputs", "lipid_categories")

all_objs <- ls()

pattern_end <- paste0(category_mapping$Category_clean, collapse = "|")

lipid_categories <- grep(paste0("^category_.*(", pattern_end, ")$"), 
                      all_objs, 
                      value = TRUE)

all_avg <- data.frame() 

for (name in lipid_categories) {
  
  df <- get(name)
  cat_total <- df[, c(1, ncol(df))]
  colnames(cat_total) <- c("group", "total")
  cat_avg <- aggregate(total ~ group, data = cat_total, FUN = mean)
  cat_avg$category <- name
  all_avg <- rbind(all_avg, cat_avg)
}

all_avg$category <- factor(all_avg$category)
all_avg$group <- factor(all_avg$group)

category_colors <- c(
  "category_fatty_acyls" = "#D56401",
  "category_glycerolipids" = "#7F0101",
  "category_glycerophospholipids" = "#017B7F",
  "category_sphingolipids" = "#057F01",
  "category_sterol_lipids" = "#7F0163"
)

cat_avg <- ggplot(all_avg, aes(x = group, y = total, fill = category)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Average by Group Across Lipid Categories",
       x = "Group",
       y = "Average Value") +
  theme_classic() +
  scale_fill_manual(values = category_colors)

ggsave(file.path(top_level_dir, "..", "total_lipids/category_averages.png"),
       plot = cat_avg,
       width = 10, height = 6, dpi = 300)

#### SAME BUT LOG TRANSFORMED ##################################################

all_avg <- data.frame()  

for (name in lipid_categories) {
  
  df <- get(name)
  cat_total <- df[, c(1, ncol(df))]
  colnames(cat_total) <- c("group", "total")
  cat_total$total <- log10(cat_total$total + 1e-6)
  cat_avg <- aggregate(total ~ group, data = cat_total, FUN = mean)
  cat_avg$category <- name
  all_avg <- rbind(all_avg, cat_avg)
}

all_avg$category <- factor(all_avg$category)
all_avg$group <- factor(all_avg$group)

cat_avg_log <- ggplot(all_avg, aes(x = group, y = total, fill = category)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Average by Group Across Lipid Categories",
       x = "Group",
       y = "Average Value") +
  theme_classic() +
  scale_fill_manual(values = category_colors)

ggsave(file.path(top_level_dir, "..", "total_lipids/category_log_averages.png"),
       plot = cat_avg_log,
       width = 10, height = 6, dpi = 300)

################################################################################
################### CORRELATION MATRIX FOR EACH CATEGORY #######################
################################################################################

category_folders <- list.dirs(top_level_dir, recursive = FALSE)

for (category_folder in category_folders) {
  
  csv_file <- list.files(category_folder, pattern = "_spearman_correlation\\.csv$", full.names = TRUE)
  
  if (length(csv_file) == 0) {
    message("No CSV found in folder: ", category_folder, " — skipping")
    next
  }
  
  cor_mat <- as.matrix(read.csv(csv_file, row.names = 1, check.names = FALSE))
  
  if (ncol(cor_mat) < 2) {
    message("Not enough lipids in correlation matrix for folder: ", category_folder, " — skipping")
    next
  }
  
  # Reorder rows/cols alphabetically
  ord <- order(colnames(cor_mat))
  cor_mat <- cor_mat[ord, ord]
  
  heatmap_file <- file.path(category_folder, paste0(basename(category_folder), "_correlation_heatmap.png"))
  
  png(filename = heatmap_file, width = 3000, height = 3000, res = 300)
  
  pheatmap(
    cor_mat,
    color = colorRampPalette(c("blue", "white", "red"))(200),
    cluster_rows = FALSE,   # keep alphabetical order
    cluster_cols = FALSE,   # keep alphabetical order
    show_rownames = TRUE,
    show_colnames = TRUE,
    main = paste("Spearman Correlation:", basename(category_folder))
  )
  
  dev.off()
  message("Saved heatmap for category: ", basename(category_folder))
}

#### FOR ALL LIPIDS COMBINED ###################################################

cor_mat <- cor(raw_data_lipids, use = "pairwise.complete.obs", method = "spearman")

# Alphabetical order for complete matrix
ord <- order(colnames(cor_mat))
cor_mat <- cor_mat[ord, ord]

write.csv(cor_mat, file = "outputs/total_lipids/complete_lipid_correlation.csv", row.names = TRUE)

print(range(cor_mat, na.rm = TRUE))

col_palette <- colorRampPalette(c("blue", "white", "red"))(200)

png("outputs/total_lipids/complete_lipid_correlation_heatmap.png",
    width = 4000, height = 4000, res = 300)

pheatmap(
  cor_mat,
  color = col_palette,
  cluster_rows = FALSE,   # alphabetical
  cluster_cols = FALSE,   # alphabetical
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "Complete Spearman Correlation (All Lipids)"
)

dev.off()



