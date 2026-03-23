run_ttest <- function(region_df, group1_pattern, group2_pattern) {
  library(dplyr)
  library(tidyr)
  library(stringr)
  
  # Step 1: Convert region_df to long format with mz, replicate, intensity
  df_long <- region_df %>%
    tibble::rownames_to_column("mz") %>%
    pivot_longer(-mz, names_to = "sample", values_to = "intensity") %>%
    mutate(replicate = sub("@.*", "", sample))
  
  # Step 2: Average intensity across replicates with same mz+replicate
  df_avg <- df_long %>%
    group_by(mz, replicate) %>%
    summarise(mean_intensity = mean(intensity, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = replicate, values_from = mean_intensity) %>%
    tibble::column_to_rownames("mz")
  
  # Step 3: Identify columns belonging to each group using pattern matching
  group1_cols <- grep(group1_pattern, colnames(df_avg), value = TRUE)
  group2_cols <- grep(group2_pattern, colnames(df_avg), value = TRUE)
  
  # Step 4: Select group columns + restore mz as a column
  df_avg_select <- df_avg %>%
    dplyr::select(all_of(c(group1_cols, group2_cols)))
  
  df_avg_select$mz <- rownames(df_avg_select)
  
  # Step 5: Convert to long format for t-test
  df_long2 <- df_avg_select %>%
    pivot_longer(cols = all_of(c(group1_cols, group2_cols)),
                 names_to = "replicate",
                 values_to = "intensity") %>%
    mutate(
      group = case_when(
        str_detect(replicate, group1_pattern) ~ group1_pattern,
        str_detect(replicate, group2_pattern) ~ group2_pattern
      )
    )
  
  # Step 6: Perform t-test across all m/z
  ttest_results <- df_long2 %>%
    filter(group %in% c(group1_pattern, group2_pattern)) %>%
    group_by(mz) %>%
    summarise(
      pvalue = tryCatch(
        t.test(intensity ~ group)$p.value,
        error = function(e) NA
      ),
      !!paste0("mean_", group1_pattern) :=
        mean(intensity[group == group1_pattern], na.rm = TRUE),
      !!paste0("mean_", group2_pattern) :=
        mean(intensity[group == group2_pattern], na.rm = TRUE),
      diff =
        mean(intensity[group == group2_pattern], na.rm = TRUE) -
        mean(intensity[group == group1_pattern], na.rm = TRUE),
      .groups = "drop"
    )
  ttest_results <- ttest_results %>% mutate(change = ifelse(diff < 0, "Decreased", "Increased"))
  
  return(ttest_results)
}
