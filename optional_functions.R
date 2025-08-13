#### add groups column #########################################################

add_groups_column <- function(raw_data_names, groups) {
  for (name in raw_data_names) {
    cat("Adding 'groups' column to:", name, "\n")
    df <- get(name)              # Get the data frame by name
    df <- cbind(groups, df)      # Add the groups column
    assign(name, df, envir = .GlobalEnv)  # Save it back into global environment
  }
}

add_groups_column(raw_data_names, groups) # test working

#### calculating group averages ################################################

calculate_group_averages <- function(raw_data_names, groups){
  for (name in raw_data_names) {
    lipid_family <- sub("^raw_data_", "", name)
    folder_path <- file.path("outputs", "lipid_families", lipid_family)
    cat("Processing group averages for:", name, "\n")
    df <- get(name)
    df_avg <- aggregate(. ~ groups, data = df, FUN = mean)
    new_name <- paste0(name, "_avg")
    assign(new_name, df_avg)
    csv_file <- file.path(folder_path, paste0(lipid_family, "_avg.csv"))
    write.csv(df_avg, file = csv_file, row.names = FALSE)
  }
}

calculate_group_averages(raw_data_names, groups)

##### ANOVAS ###################################################################

## for individual lipid families ###############################################

anova_per_lipid_family <- function(raw_data_names, groups){
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
}

anova_per_lipid_family <- function(raw_data_names, groups)

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



































