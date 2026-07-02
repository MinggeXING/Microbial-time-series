# Get environment stage data ----

rm(list = ls())

# Load required packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
})

# 1. Environment data preprocessing ----
# 1.1 Filter odd time point samples from original environment data ----

# Read original environment data
env_original <- read.csv("data/540/env.csv", check.names = FALSE)

# Filter samples at odd time points (T1, T3, T5, ..., T29)
env_filtered <- env_original %>%
  filter(str_detect(Sample_ID, "^T(1|3|5|7|9|11|13|15|17|19|21|23|25|27|29)[AMNM].*_\\d+$"))

message(sprintf("Retained %d environment samples after filtering", nrow(env_filtered)))

# Save filtered environment data
write.csv(env_filtered, "data/stage_data/env_filtered_stage_equ.csv", row.names = FALSE)

# 1.2 Reorder rows and rename groups ----

# Re-read filtered environment data
env_data <- read.csv("data/stage_data/env_filtered_stage_equ.csv", check.names = FALSE)

# ~~Data sorting preparation ----

# Define time period priority groups (consistent with ASV table)
time_groups <- list(
  c(1, 3, 5),    # First priority group - stage 1
  c(7, 9, 11),   # Second priority group - stage 2
  c(13, 15, 17), # Third priority group - stage 3
  c(19, 21, 23), # Fourth priority group - stage 4
  c(25, 27, 29)  # Fifth priority group - stage 5
)

# Create time priority mapping table
time_priority <- stack(setNames(rep(1:5, lengths(time_groups)), unlist(time_groups))) %>%
  setNames(c("group_priority", "time")) %>%
  mutate(time = as.numeric(as.character(time)))

# Define treatment group priority (strictly AM0→NM0→AM2→NM2→AM5→NM5→AM15→NM15)
treatment_order <- c("AM0", "NM0", "AM2", "NM2", "AM5", "NM5", "AM15", "NM15")

# ~~Data sorting ----

# Extract Sample_ID and parse sorting elements
env_data <- env_data %>%
  mutate(
    # Extract time point (number after T)
    time = as.numeric(str_extract(Sample_ID, "(?<=T)\\d+")),
    # Extract treatment type (letters after T and number)
    treat = str_extract(Sample_ID, "(?<=T\\d{1,2})[A-Z]+"),
    # Extract group number (number after letters, until underscore)
    group = as.numeric(str_extract(Sample_ID, "(?<=[A-Z]{2})\\d+")),
    # Extract replicate number (number after underscore)
    rep = as.numeric(str_extract(Sample_ID, "(?<=_)\\d+$")),
    treatment = paste0(treat, group)  # Create AM0/AM2 combinations
  ) %>%
  left_join(time_priority, by = "time") %>%
  mutate(
    treatment_priority = match(treatment, treatment_order) # Assign treatment group priority
  )

# Check parsing results
cat("Preview of parsed data:\n")
print(head(env_data %>% select(Sample_ID, time, treat, group, rep, treatment, group_priority, treatment_priority)))

# Perform multi-level sorting
env_sorted <- env_data %>%
  arrange(
    group_priority,    # Level 1: time period group priority
    treatment_priority, # Level 2: treatment group priority
    time,              # Level 3: time point within group
    rep                # Level 4: replicate number
  )

# ~~Rename Sample_ID ----

# Create mapping from time point to stage number
time_to_stage <- list(
  "1" = "1", "3" = "1", "5" = "1",
  "7" = "2", "9" = "2", "11" = "2",
  "13" = "3", "15" = "3", "17" = "3",
  "19" = "4", "21" = "4", "23" = "4",
  "25" = "5", "27" = "5", "29" = "5"
)

# Precise replacement function
replace_time_to_stage <- function(sample_id) {
  # Extract time point number - use simpler method
  time_point <- str_extract(sample_id, "T(\\d+)") %>% str_remove("T")
  
  # Get corresponding stage number
  stage <- time_to_stage[[time_point]]
  
  if(is.null(stage)) {
    warning(paste("No mapping found for time point", time_point))
    return(sample_id)
  }
  
  # Replace T with S and stage number
  str_replace(sample_id, "T\\d+", paste0("S", stage))
}

# Apply precise replacement
new_sample_ids <- sapply(env_sorted$Sample_ID, replace_time_to_stage)

# ~~Handle duplicate Sample_ID issues ----

# Define treatment types
treatment_types <- c("AM0", "NM0", "AM2", "NM2", "AM5", "NM5", "AM15", "NM15")

# Generate all possible group combinations
all_groups <- paste0(rep(c("S1", "S2", "S3", "S4", "S5"), each = length(treatment_types)), treatment_types)

# Initialize new Sample_ID vector
final_sample_ids <- character(length(new_sample_ids))

# Continuous numbering for each group
for(group in all_groups){
  # Find all samples in current group
  group_samples <- which(str_detect(new_sample_ids, paste0("^", group)))
  
  if(length(group_samples) > 0){
    # Generate continuous numbers (1 to n)
    new_suffix <- seq_along(group_samples)
    # Build new Sample_ID
    final_sample_ids[group_samples] <- paste0(group, "_", new_suffix)
  }
}

# ~~Apply final Sample_ID and validate ----

# Apply new Sample_ID
env_sorted$Sample_ID <- final_sample_ids

# Remove temporary columns
env_final <- env_sorted %>%
  select(-time, -treat, -group, -rep, -treatment, -group_priority, -treatment_priority)

# Validate results
cat("\nProcessed environment sample IDs:\n")
print(head(env_final$Sample_ID, 20))

# Check for duplicate Sample_ID
if(any(duplicated(env_final$Sample_ID))){
  warning("Duplicate Sample_ID found: ", paste(env_final$Sample_ID[duplicated(env_final$Sample_ID)], collapse = ", "))
} else {
  cat("Sample_ID uniqueness validation passed, no duplicate sample IDs\n")
}

# Check final data dimensions
message("Final environment data dimensions:")
message("- Samples: ", nrow(env_final))
message("- Variables: ", ncol(env_final))

# Save final environment data
write.csv(env_final, "data/stage_data/env_stage.csv", row.names = FALSE)

# Output statistics
cat("\nEnvironment data preprocessing completed!\n")
cat("Input samples:", nrow(env_original), "\n")
cat("Output samples:", nrow(env_final), "\n")
cat("Sample ID format example:", paste(head(env_final$Sample_ID, 3), collapse = ", "), "\n")

rm(list = ls())

##################################################################################
# 2. Correlation analysis between functional genera and environmental factors for stage data (grouped by Cd gradient) ----
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)

# Create output directories
dir.create("out/microbe/corr_heatmap", recursive = TRUE, showWarnings = FALSE)
dir.create("data/stage_data", recursive = TRUE, showWarnings = FALSE)

genus_stage_abs_target <- read.csv("data/stage_data/target_genera_abs_stage.csv", header = TRUE, row.names = 1)
stage_metadata <- read.csv("data/stage_data/stage_metadata.csv", header = TRUE)
env <- read.csv("data/stage_data/env_stage.csv", header = TRUE, row.names = 1)
env <- env[, 3:16]  # Remove first two non-environment columns

target_genera <- c(
  "Massilia", "Streptomyces", "Clostridium", "Flavisolibacter", "Gemmatimonas", 
  "Lysobacter", "Nocardioides", "Paenibacillus","Bacillus", "Pseudomonas"
)

# Define function to remove zero columns
rem_0 <- function(data, num = 0) {
  data[, colSums(data != num) > 0]
}

# Merge genus and environment data -----
common_samples <- intersect(rownames(genus_stage_abs_target), rownames(env))
genus_data <- genus_stage_abs_target[common_samples, target_genera, drop = FALSE]
env_data <- env[common_samples, ]

# Get Cd information from stage_metadata
cd_info <- stage_metadata[match(common_samples, stage_metadata$Sample_ID), "Cd", drop = FALSE]
rownames(cd_info) <- common_samples

# Define Cd gradients
cd_gradients <- c(0, 2, 5, 15)

# Create list to store all correlation results
all_cor_results <- list()

# Generate correlation heatmap for each Cd gradient
for (cd_level in cd_gradients) {
  cat("Processing Cd gradient:", cd_level, "\n")
  
  # Filter samples for current Cd gradient
  cd_samples <- rownames(cd_info)[cd_info$Cd == cd_level]
  
  if (length(cd_samples) < 3) {
    cat("Cd gradient", cd_level, "has insufficient samples, skipping\n")
    next
  }
  
  # Extract genus and environment data for current gradient
  genus_cd <- genus_data[cd_samples, , drop = FALSE]
  env_cd <- env_data[cd_samples, , drop = FALSE]
  
  # Remove all-zero columns
  genus_cd <- rem_0(genus_cd, num = 0)
  env_cd <- rem_0(env_cd, num = 0)
  
  # Select numeric environmental factors
  numeric_env <- env_cd[, sapply(env_cd, is.numeric), drop = FALSE]
  
  # Check if data is valid
  if (ncol(genus_cd) == 0 || ncol(numeric_env) == 0) {
    cat("Cd gradient", cd_level, "has invalid data, skipping\n")
    next
  }
  
  # Calculate correlation matrix between genera and environmental factors
  cor_matrix <- matrix(NA, nrow = ncol(genus_cd), ncol = ncol(numeric_env))
  p_value_matrix <- matrix(NA, nrow = ncol(genus_cd), ncol = ncol(numeric_env))
  
  rownames(cor_matrix) <- colnames(genus_cd)
  colnames(cor_matrix) <- colnames(numeric_env)
  rownames(p_value_matrix) <- colnames(genus_cd)
  colnames(p_value_matrix) <- colnames(numeric_env)
  
  # Calculate Spearman correlation for each genus and environmental factor
  for (i in 1:ncol(genus_cd)) {
    for (j in 1:ncol(numeric_env)) {
      cor_test <- cor.test(genus_cd[, i], numeric_env[, j], 
                           method = "spearman", 
                           exact = FALSE)
      cor_matrix[i, j] <- cor_test$estimate
      p_value_matrix[i, j] <- cor_test$p.value
    }
  }
  
  # Convert correlation matrix to long format for ggplot2
  cor_melted <- melt(cor_matrix)
  p_melted <- melt(p_value_matrix)
  
  # Merge correlation and p-values
  plot_data <- cbind(cor_melted, p_value = p_melted$value)
  colnames(plot_data) <- c("Genus", "Env_Factor", "Correlation", "p_value")
  
  # Add Cd gradient information
  plot_data$Cd_level <- cd_level
  
  # Add significance markers
  plot_data$signif <- ""
  plot_data$signif[plot_data$p_value < 0.001] <- "***"
  plot_data$signif[plot_data$p_value >= 0.001 & plot_data$p_value < 0.01] <- "**"
  plot_data$signif[plot_data$p_value >= 0.01 & plot_data$p_value < 0.05] <- "*"
  
  # Store results
  all_cor_results[[as.character(cd_level)]] <- plot_data
  
  # Plot heatmap using ggplot2
  p <- ggplot(plot_data, aes(x = Env_Factor, y = Genus, fill = Correlation)) +
    geom_tile(color = "white", linewidth = 0.5) +
    # Correlation values
    geom_text(aes(label = sprintf("%.2f", Correlation)), 
              size = 5, color = "black") +
    # Significance markers
    geom_text(aes(label = signif), 
              size = 6, color = "black", fontface = "bold",
              nudge_y = 0.25) +  
    scale_fill_gradient2(low = "#64DDCE", high = "#EA6A47", mid = "white", 
                         midpoint = 0, limit = c(-1, 1), 
                         name = "Correlation\n") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 13, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 11),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank()
    ) +
    labs(
      title = paste("Cd_", cd_level),
      x = "", 
      y = "Genera"
    ) +
    coord_fixed(ratio = 1)
  
  # Save plots to specified path
  ggsave(paste0("out/microbe/corr_heatmap/correlation_heatmap_Cd", cd_level, ".pdf"), p, width = 12, height = 8)
  ggsave(paste0("out/microbe/corr_heatmap/correlation_heatmap_Cd", cd_level, ".jpg"), p, width = 12, height = 8, dpi = 300)
  
}

# Combine correlation results from all Cd gradients ----
if (length(all_cor_results) > 0) {
  # Combine all data
  combined_cor <- do.call(rbind, all_cor_results)
  
  # Reorder columns
  combined_cor <- combined_cor %>%
    select(Cd_level, Genus, Env_Factor, Correlation, p_value, signif)
  
  # Save combined data to specified path
  write.csv(
    combined_cor, 
    "data/stage_data/function_bacteria_env_corr.csv", 
    row.names = FALSE
  )
  
  
} else {
  cat("Insufficient data to generate results\n")
}


################################################################################
# 3. Data correlation and causality calculation ----

# Create output directories
dir.create("out/microbe/causality_heatmap", recursive = TRUE, showWarnings = FALSE)
dir.create("data/540", recursive = TRUE, showWarnings = FALSE)

# Load required packages
library(rEDM)
library(boot)
library(Kendall)
library(ggplot2)
library(reshape2)
library(dplyr)
library(patchwork)

# Read data
genus_abs_target <- read.csv("data/540/target_genera_abundance.csv", header = TRUE, row.names = 1)
metadata <- read.csv("data/540/metadata.csv", header = TRUE)
env <- read.csv("data/540/env.csv", header = TRUE, row.names = 1)
env <- env[, 4:18]  # Remove first three non-environment columns

# Define target genus order
target_genera_ordered <- c(
  "Massilia", "Streptomyces", "Clostridium", "Flavisolibacter", "Gemmatimonas", 
  "Lysobacter", "Nocardioides", "Paenibacillus", "Bacillus", "Pseudomonas",
  "Sphingomonas", "Microvirga", "Nitrospira"
)

# Define environmental factor order
specified_env_order <- c(
  "AP", "NH4_N", "NO3_N", "NO2_N", "pH", "soil_tem", "soil_moi", 
  "Pn", "Tr", "WUE", "env_tem_ave", "env_moi_ave"
)  

# Check if environmental factor column names exist
cat("Checking environmental factor columns...\n")
available_env <- colnames(env)
cat("Available environmental factors:\n")
print(available_env)

# Keep only existing environmental factors
env_order_final <- specified_env_order[specified_env_order %in% available_env]
if (length(env_order_final) < length(specified_env_order)) {
  missing <- setdiff(specified_env_order, available_env)
  cat("Warning: Some specified environmental factors are missing:\n")
  print(missing)
  cat("Using available factors only.\n")
}

# Reorder environmental factor data frame
env <- env[, env_order_final, drop = FALSE]

# Define environmental factors that need reverse calculation (Genus → Environment)
reverse_env_factors <- c("AP", "NH4_N", "NO3_N", "NO2_N")

# Define environmental factors for forward calculation (Environment → Genus)
forward_env_factors <- setdiff(env_order_final, reverse_env_factors)

cat("\nEnvironmental factor directions:\n")
cat("Reverse direction (Genus → Env):", paste(reverse_env_factors, collapse = ", "), "\n")
cat("Forward direction (Env → Genus):", paste(forward_env_factors, collapse = ", "), "\n")

# Define function to remove zero columns
rem_0 <- function(data, num = 0) {
  data[, colSums(data != num) > 0]
}

# Merge genus data and environment data -----
common_samples <- intersect(rownames(genus_abs_target), rownames(env))
genus_data <- genus_abs_target[common_samples, target_genera_ordered, drop = FALSE]
env_data <- env[common_samples, , drop = FALSE]

# Get Cd, Mycorrhizal and Time information from metadata
meta_info <- metadata[match(common_samples, metadata$Sample_ID), 
                      c("Cd", "Mycorrhizal", "Time", "New_time"), drop = FALSE]
rownames(meta_info) <- common_samples

# Filter AM group samples
am_samples <- rownames(meta_info)[meta_info$Mycorrhizal == "AM"]
cat("Total AM samples:", length(am_samples), "\n")

if (length(am_samples) == 0) {
  stop("No AM samples found! Please check metadata.")
}

# Use only AM group data
genus_data_am <- genus_data[am_samples, target_genera_ordered, drop = FALSE]
env_data_am <- env_data[am_samples, env_order_final, drop = FALSE]
meta_info_am <- meta_info[am_samples, , drop = FALSE]

# Define Cd gradients
cd_gradients <- c(0, 2, 5, 15)

################################################################################
# Time series validation function ----
validate_time_series <- function(time_vector) {
  # Check if numeric
  if (!is.numeric(time_vector)) {
    time_vector <- as.numeric(as.character(time_vector))
    if (any(is.na(time_vector))) {
      return(list(
        is_numeric = FALSE,
        is_increasing = FALSE,
        has_na = TRUE,
        has_duplicates = NA,
        length = length(time_vector),
        min_time = NA,
        max_time = NA,
        message = "Time vector contains non-numeric values"
      ))
    }
  }
  
  # Check for missing values
  has_na <- any(is.na(time_vector))
  
  # Check if monotonically increasing
  is_increasing <- all(diff(time_vector) >= 0)
  
  # Check for duplicate values
  has_duplicates <- any(duplicated(time_vector))
  
  return(list(
    is_numeric = TRUE,
    is_increasing = is_increasing,
    has_na = has_na,
    has_duplicates = has_duplicates,
    length = length(time_vector),
    min_time = min(time_vector, na.rm = TRUE),
    max_time = max(time_vector, na.rm = TRUE)
  ))
}

# Validate time information
cat("\n=== Validating time series for AM samples ===\n")
for (cd in cd_gradients) {
  cd_samples <- rownames(meta_info_am)[meta_info_am$Cd == cd]
  if (length(cd_samples) > 0) {
    # Use New_time for sorting
    time_vec <- meta_info_am[cd_samples, "New_time"]
    validation <- validate_time_series(time_vec)
    
    cat("\nCd =", cd, ":\n")
    cat("  Samples:", length(cd_samples), "\n")
    cat("  Numeric:", validation$is_numeric, "\n")
    cat("  Increasing:", validation$is_increasing, "\n")
    cat("  Has NA:", validation$has_na, "\n")
    cat("  Has duplicates:", validation$has_duplicates, "\n")
    cat("  New_time range:", validation$min_time, "-", validation$max_time, "\n")
  }
}

################################################################################
# CCM analysis function - bidirectional ----
ccm.fast.demo.improved.bidirectional <- function(ds, Epair = TRUE, cri = 'rmse', Emax = 10, 
                                                 alpha = 0.05, min_lib_size = 3, 
                                                 bootstrap_reps = 1000, 
                                                 n_genus = NULL, n_env = NULL,
                                                 reverse_indices = NULL) {
  
  cat("Starting BIDIRECTIONAL CCM analysis with", nrow(ds), "time points and", ncol(ds), "variables\n")
  
  if (cri == 'rho') {
    jfun <- match.fun('which.max')
  } else {
    jfun <- match.fun('which.min')
  }
  
  # Z-score standardization
  ds <- as.matrix(apply(ds, 2, scale))
  np <- nrow(ds) # time series length
  ns <- ncol(ds) # number of variables
  
  # Check if time series length is sufficient
  if (np < 10) {
    cat("WARNING: Time series length (", np, ") is short for CCM analysis\n")
  }
  
  # Library sizes - gradually increasing from minimum to maximum
  lib.s <- c(seq(10, min(50, np), by = 10), np)
  lib.s <- unique(lib.s[lib.s <= np])
  
  # Critical rho value based on specified alpha level
  crirho <- qt(1 - alpha/2, np - 2) / sqrt(np - 2 + qt(1 - alpha/2, np - 2)^2)
  
  # Initialize result matrices
  ccm.rho <- ccm.sig <- ccm.kendall.p <- matrix(0, ns, ns)
  ccm.boot.p <- matrix(1, ns, ns) # initialize bootstrap p-value to 1 (not significant)
  
  # If numbers of genera and environmental factors are specified, we can optimize the loop
  if (!is.null(n_genus) && !is.null(n_env) && (n_genus + n_env) == ns) {
    cat("Optimized BIDIRECTIONAL CCM calculation for", n_genus, "genera and", n_env, "environmental factors\n")
    
    # Analyze forward direction: Environment → Genus
    for (env_idx in 1:n_env) {
      i <- n_genus + env_idx  # environmental factor index (in the second half of the data matrix)
      env_name <- colnames(ds)[i]
      
      t.begin <- proc.time()
      
      # Check the direction of this environmental factor
      if (!is.null(reverse_indices) && env_idx %in% reverse_indices) {
        # Reverse: Genus → Environment
        cat("  Processing reverse direction: Genus →", env_name, "\n")
        
        for (j in 1:n_genus) {  # j: genus index
          
          # Select optimal E value
          ccm.E <- NULL  
          for (E.t in 2:min(Emax, np-2)) {
            ccm.E <- rbind(ccm.E, ccm(cbind(x = ds[, j], y = ds[, i]), E = E.t, tp = -1,
                                      lib_column = "x", target_column = "y", 
                                      lib_sizes = np, random_libs = FALSE, silent = TRUE))
          }
          
          if (is.null(ccm.E) || nrow(ccm.E) == 0) {
            next
          }
          
          Eop <- ccm.E[jfun(ccm.E[, cri]), 'E']
          
          # Perform CCM analysis - genus j → environmental factor i
          ccm.out <- ccm(cbind(x = ds[, j], y = ds[, i]), E = Eop, tp = 0, 
                         lib_column = "x", target_column = "y", 
                         lib_sizes = lib.s, random_libs = FALSE, silent = TRUE)
          
          # Aggregate results
          ccm.seq <- aggregate(ccm.out[, 'rho'], list(ccm.out[, 'lib_size']), mean, na.rm = TRUE)
          ccm.seq <- ccm.seq[!(is.na(ccm.seq[, 2]) | is.infinite(ccm.seq[, 2])), ]
          
          if (nrow(ccm.seq) == 0) {
            next
          }
          
          ccm.seq[ccm.seq[, 2] < 0, 2] <- 0
          termrho <- ccm.seq[nrow(ccm.seq), 2]  # terminal rho value
          
          # Calculate bootstrap p-value
          if (bootstrap_reps > 0) {
            boot_rhos <- numeric(bootstrap_reps)
            successful_boots <- 0
            
            for (b in 1:bootstrap_reps) {
              # Randomly shuffle y variable (environmental factor), destroying causality but preserving time structure
              y_shuffled <- sample(ds[, i])
              
              ccm.boot <- ccm(cbind(x = ds[, j], y = y_shuffled), E = Eop, tp = 0,
                              lib_column = "x", target_column = "y", 
                              lib_sizes = np, random_libs = FALSE, silent = TRUE)
              
              if (!is.null(ccm.boot) && !all(is.na(ccm.boot$rho))) {
                boot_rhos[b] <- mean(ccm.boot$rho, na.rm = TRUE)
                successful_boots <- successful_boots + 1
              }
            }
            
            if (successful_boots > 0) {
              boot_rhos <- boot_rhos[1:successful_boots]
              ccm.boot.p[j, i] <- sum(boot_rhos >= termrho) / successful_boots
            }
          }
          
          # Mann-Kendall test
          if (nrow(ccm.seq) >= min_lib_size) {
            kend <- MannKendall(ccm.seq[, 2])
            ccm.kendall.p[j, i] <- ifelse(is.null(kend$sl[1]), 1, kend$sl[1])
            
            # Dual significance test
            ccm.sig[j, i] <- (kend$tau[1] > 0) * 
              (kend$sl[1] < alpha) * 
              (termrho > crirho) *
              (ccm.boot.p[j, i] < alpha)
          } else {
            ccm.kendall.p[j, i] <- 1
            ccm.sig[j, i] <- 0
          }
          ccm.rho[j, i] <- termrho
        }
      } else {
        # Forward: Environment → Genus
        cat("  Processing forward direction:", env_name, "→ Genus\n")
        
        for (j in 1:n_genus) {     # genus (in the first half of the data matrix)
          
          # Select optimal E value
          ccm.E <- NULL  
          for (E.t in 2:min(Emax, np-2)) {
            ccm.E <- rbind(ccm.E, ccm(cbind(x = ds[, i], y = ds[, j]), E = E.t, tp = -1,
                                      lib_column = "x", target_column = "y", 
                                      lib_sizes = np, random_libs = FALSE, silent = TRUE))
          }
          
          if (is.null(ccm.E) || nrow(ccm.E) == 0) {
            next
          }
          
          Eop <- ccm.E[jfun(ccm.E[, cri]), 'E']
          
          # Perform CCM analysis - environmental factor i → genus j
          ccm.out <- ccm(cbind(x = ds[, i], y = ds[, j]), E = Eop, tp = 0, 
                         lib_column = "x", target_column = "y", 
                         lib_sizes = lib.s, random_libs = FALSE, silent = TRUE)
          
          # Aggregate results
          ccm.seq <- aggregate(ccm.out[, 'rho'], list(ccm.out[, 'lib_size']), mean, na.rm = TRUE)
          ccm.seq <- ccm.seq[!(is.na(ccm.seq[, 2]) | is.infinite(ccm.seq[, 2])), ]
          
          if (nrow(ccm.seq) == 0) {
            next
          }
          
          ccm.seq[ccm.seq[, 2] < 0, 2] <- 0
          termrho <- ccm.seq[nrow(ccm.seq), 2]  # terminal rho value
          
          # Calculate bootstrap p-value
          if (bootstrap_reps > 0) {
            boot_rhos <- numeric(bootstrap_reps)
            successful_boots <- 0
            
            for (b in 1:bootstrap_reps) {
              # Randomly shuffle y variable (genus), destroying causality but preserving time structure
              y_shuffled <- sample(ds[, j])
              
              ccm.boot <- ccm(cbind(x = ds[, i], y = y_shuffled), E = Eop, tp = 0,
                              lib_column = "x", target_column = "y", 
                              lib_sizes = np, random_libs = FALSE, silent = TRUE)
              
              if (!is.null(ccm.boot) && !all(is.na(ccm.boot$rho))) {
                boot_rhos[b] <- mean(ccm.boot$rho, na.rm = TRUE)
                successful_boots <- successful_boots + 1
              }
            }
            
            if (successful_boots > 0) {
              boot_rhos <- boot_rhos[1:successful_boots]
              ccm.boot.p[i, j] <- sum(boot_rhos >= termrho) / successful_boots
            }
          }
          
          # Mann-Kendall test
          if (nrow(ccm.seq) >= min_lib_size) {
            kend <- MannKendall(ccm.seq[, 2])
            ccm.kendall.p[i, j] <- ifelse(is.null(kend$sl[1]), 1, kend$sl[1])
            
            # Dual significance test
            ccm.sig[i, j] <- (kend$tau[1] > 0) * 
              (kend$sl[1] < alpha) * 
              (termrho > crirho) *
              (ccm.boot.p[i, j] < alpha)
          } else {
            ccm.kendall.p[i, j] <- 1
            ccm.sig[i, j] <- 0
          }
          ccm.rho[i, j] <- termrho
        }
      }
      time.used <- proc.time() - t.begin 
      cat("  Environmental factor", env_name, "CCM completed:", time.used[3], "sec\n")
    }
  } else {
    cat("error")
  }
  
  return(list(ccm.rho = ccm.rho, 
              ccm.sig = ccm.sig, 
              ccm.kendall.p = ccm.kendall.p,
              ccm.boot.p = ccm.boot.p))
}

################################################################################
# Bootstrap correlation analysis function ----
bootstrap_correlation <- function(x, y, nboot = 1000) {
  # Calculate observed correlation
  obs_cor <- cor(x, y, method = "spearman")
  
  # Check if there is sufficient valid data
  if (is.na(obs_cor) || length(x) < 5 || length(y) < 5) {
    return(list(correlation = NA,
                p_value = 1,
                ci_lower = NA,
                ci_upper = NA))
  }
  
  # Bootstrap function
  boot_func <- function(data, indices) {
    cor(data[indices, 1], data[indices, 2], method = "spearman")
  }
  
  # Perform bootstrap
  boot_data <- cbind(x, y)
  boot_data <- boot_data[complete.cases(boot_data), , drop = FALSE]
  
  if (nrow(boot_data) < 5) {
    return(list(correlation = obs_cor,
                p_value = 1,
                ci_lower = NA,
                ci_upper = NA))
  }
  
  tryCatch({
    boot_results <- boot(boot_data, boot_func, R = nboot)
    
    # Calculate confidence intervals
    ci <- tryCatch({
      boot.ci(boot_results, type = "bca")
    }, error = function(e) {
      boot.ci(boot_results, type = "perc")
    })
    
    # Calculate p-value
    if (!is.null(ci) && "bca" %in% names(ci)) {
      p_value <- ifelse(ci$bca[4] > 0 & ci$bca[5] > 0 | 
                          ci$bca[4] < 0 & ci$bca[5] < 0, 0, 1)
      ci_lower <- ci$bca[4]
      ci_upper <- ci$bca[5]
    } else if (!is.null(ci) && "percent" %in% names(ci)) {
      p_value <- ifelse(ci$percent[4] > 0 & ci$percent[5] > 0 | 
                          ci$percent[4] < 0 & ci$percent[5] < 0, 0, 1)
      ci_lower <- ci$percent[4]
      ci_upper <- ci$percent[5]
    } else {
      p_value <- 1
      ci_lower <- NA
      ci_upper <- NA
    }
    
    return(list(correlation = obs_cor,
                p_value = p_value,
                ci_lower = ci_lower,
                ci_upper = ci_upper))
  }, error = function(e) {
    return(list(correlation = obs_cor,
                p_value = 1,
                ci_lower = NA,
                ci_upper = NA))
  })
}

################################################################################
# Main function for analyzing by Cd group - using AM group data, sorted by New_time (1~90) ----
analyze_by_cd_am_with_time <- function(genus_data, env_data, meta_info, cd_gradients) {
  
  # Store all results
  all_cor_results <- list()
  all_ccm_results <- list()
  
  for (cd in cd_gradients) {
    cat("\n" , strrep("=", 50), "\n")
    cat("Analyzing AM group - Cd =", cd, "\n")
    cat(strrep("=", 50), "\n")
    
    # Filter samples for current Cd concentration
    cd_samples <- rownames(meta_info)[meta_info$Cd == cd]
    
    if (length(cd_samples) < 10) {  # Increase minimum sample requirement
      cat("  Not enough AM samples for Cd =", cd, "(", length(cd_samples), "samples), skipping...\n")
      cat("  Minimum 10 samples required for reliable CCM analysis\n")
      next
    }
    
    # Key step: sort by New_time (1~90)
    cd_meta <- meta_info[cd_samples, ]
    
    # Ensure time information is numeric
    cd_meta$New_time <- as.numeric(as.character(cd_meta$New_time))
    
    # Handle missing time information
    if (any(is.na(cd_meta$New_time))) {
      cat("  WARNING: Some samples have missing New_time values\n")
      # Use Time as fallback
      if (all(is.na(cd_meta$New_time)) && !all(is.na(cd_meta$Time))) {
        cat("  Using Time for sorting\n")
        cd_meta$New_time <- as.numeric(as.character(cd_meta$Time))
      } else {
        cat("  No valid time information, using sample order\n")
        cd_meta$New_time <- 1:nrow(cd_meta)
      }
    }
    
    # Sort by New_time
    cd_meta <- cd_meta[order(cd_meta$New_time), ]
    cd_samples_sorted <- rownames(cd_meta)
    
    cat("  Number of samples:", length(cd_samples_sorted), "\n")
    cat("  New_time range:", min(cd_meta$New_time, na.rm = TRUE), "to", max(cd_meta$New_time, na.rm = TRUE), "\n")
    cat("  New_time points (first 10):", paste(head(cd_meta$New_time, 10), collapse = ", "), "\n")
    if (length(cd_meta$New_time) > 10) {
      cat("  New_time points (last 10):", paste(tail(cd_meta$New_time, 10), collapse = ", "), "\n")
    }
    
    # Check time series quality
    time_validation <- validate_time_series(cd_meta$New_time)
    if (!time_validation$is_increasing) {
      cat("  WARNING: New_time series is not strictly increasing!\n")
    }
    if (time_validation$has_duplicates) {
      cat("  WARNING: New_time series has duplicate time points!\n")
    }
    
    # Extract data sorted by time
    cd_genus <- genus_data[cd_samples_sorted, target_genera_ordered, drop = FALSE]
    cd_env <- env_data[cd_samples_sorted, env_order_final, drop = FALSE]
    
    # Remove all-zero columns
    cd_genus <- rem_0(cd_genus)
    cd_env <- rem_0(cd_env)
    
    # Ensure correct column order
    cd_genus <- cd_genus[, intersect(target_genera_ordered, colnames(cd_genus)), drop = FALSE]
    cd_env <- cd_env[, intersect(env_order_final, colnames(cd_env)), drop = FALSE]
    
    if (ncol(cd_genus) == 0 | ncol(cd_env) == 0) {
      cat("  No valid data for AM group Cd =", cd, ", skipping...\n")
      next
    }
    
    # Check time series length
    n_time_points <- nrow(cd_genus)
    cat("  Time series length:", n_time_points, "\n")
    cat("  Number of bacterial genera:", ncol(cd_genus), "\n")
    cat("  Number of environmental factors:", ncol(cd_env), "\n")
    
    if (n_time_points < 10) {
      cat("  WARNING: Time series too short for reliable CCM analysis (min 10 points recommended)\n")
    }
    
    # Z-score standardization (independent standardization within each Cd group)
    cat("  Performing Z-score standardization...\n")
    cd_genus_scaled <- as.data.frame(scale(cd_genus))
    cd_env_scaled <- as.data.frame(scale(cd_env))
    
    # Handle NA values after standardization
    cd_genus_scaled[is.na(cd_genus_scaled)] <- 0
    cd_env_scaled[is.na(cd_env_scaled)] <- 0
    
    # Ensure correct column order
    cd_genus_scaled <- cd_genus_scaled[, colnames(cd_genus), drop = FALSE]
    cd_env_scaled <- cd_env_scaled[, colnames(cd_env), drop = FALSE]
    
    # Store current Cd group results
    cor_matrix <- matrix(NA, nrow = ncol(cd_genus_scaled), ncol = ncol(cd_env_scaled))
    cor_p_matrix <- matrix(NA, nrow = ncol(cd_genus_scaled), ncol = ncol(cd_env_scaled))
    rownames(cor_matrix) <- colnames(cd_genus_scaled)
    colnames(cor_matrix) <- colnames(cd_env_scaled)
    rownames(cor_p_matrix) <- colnames(cd_genus_scaled)
    colnames(cor_p_matrix) <- colnames(cd_env_scaled)
    
    # Perform bootstrap correlation analysis
    cat("  Performing bootstrap correlation analysis...\n")
    correlation_count <- 0
    for (i in 1:ncol(cd_genus_scaled)) {
      for (j in 1:ncol(cd_env_scaled)) {
        genus_vec <- cd_genus_scaled[, i]
        env_vec <- cd_env_scaled[, j]
        
        # Remove NA values
        valid_idx <- !is.na(genus_vec) & !is.na(env_vec)
        if (sum(valid_idx) >= 10) {  # Increase minimum valid point requirement
          cor_result <- bootstrap_correlation(genus_vec[valid_idx], env_vec[valid_idx])
          cor_matrix[i, j] <- cor_result$correlation
          cor_p_matrix[i, j] <- cor_result$p_value
          correlation_count <- correlation_count + 1
        }
      }
    }
    cat("  Completed", correlation_count, "correlation calculations\n")
    
    # Perform CCM causality analysis
    cat("  Performing BIDIRECTIONAL CCM causality analysis...\n")
    
    # Combine genus and environment data for CCM - note order: genera first, environmental factors second
    combined_data <- cbind(cd_genus_scaled, cd_env_scaled)
    
    # Determine indices of reverse environmental factors
    reverse_indices <- which(colnames(cd_env_scaled) %in% reverse_env_factors)
    
    cat("  Reverse CCM indices (Genus → Env):", reverse_indices, "\n")
    cat("  Reverse CCM factors:", colnames(cd_env_scaled)[reverse_indices], "\n")
    
    # Verify data is time-ordered
    cat("  Verifying time order...\n")
    cat("    First 5 New_time points:", head(cd_meta$New_time, 5), "\n")
    cat("    Last 5 New_time points:", tail(cd_meta$New_time, 5), "\n")
    
    # Use new bidirectional CCM function
    n_genus <- ncol(cd_genus_scaled)
    n_env <- ncol(cd_env_scaled)
    
    cat("  Starting bidirectional CCM with", nrow(combined_data), "time points\n")
    ccm_result <- ccm.fast.demo.improved.bidirectional(combined_data, 
                                                       bootstrap_reps = 200,  # reduce to speed up
                                                       n_genus = n_genus,
                                                       n_env = n_env,
                                                       reverse_indices = reverse_indices,
                                                       Emax = min(8, nrow(combined_data)-2))  # limit maximum embedding dimension
    
    # Store results
    all_cor_results[[paste0("Cd_", cd)]] <- list(
      correlation = cor_matrix,
      p_value = cor_p_matrix,
      time_info = cd_meta$New_time,  # save New_time for validation
      sample_order = cd_samples_sorted,
      n_samples = n_time_points
    )
    
    all_ccm_results[[paste0("Cd_", cd)]] <- ccm_result
    
    # Save current Cd group results
    saveRDS(list(
      correlation = all_cor_results[[paste0("Cd_", cd)]],
      causality = all_ccm_results[[paste0("Cd_", cd)]],
      time_order = cd_meta$New_time,
      samples = cd_samples_sorted,
      n_genus = n_genus,
      n_env = n_env,
      reverse_factors = reverse_env_factors
    ), file = paste0("out/microbe/causality_heatmap/am_cd_", cd, "_results.rds"))
    
    cat("  Analysis completed for Cd =", cd, "\n")
  }
  
  return(list(correlation = all_cor_results, causality = all_ccm_results))
}

################################################################################
# Execute analysis ----
cat("\n" , strrep("*", 60), "\n")
cat("Starting AM group BIDIRECTIONAL CCM analysis with New_time ordering\n")
cat(strrep("*", 60), "\n\n")

results_am <- analyze_by_cd_am_with_time(genus_data_am, env_data_am, meta_info_am, cd_gradients)

# Verify time order in results
cat("\n" , strrep("=", 60), "\n")
cat("Verifying New_time order in results\n")
cat(strrep("=", 60), "\n")

for (cd in cd_gradients) {
  cd_name <- paste0("Cd_", cd)
  if (cd_name %in% names(results_am$correlation)) {
    time_info <- results_am$correlation[[cd_name]]$time_info
    n_samples <- results_am$correlation[[cd_name]]$n_samples
    cat("\nCd =", cd, "\n")
    cat("  Samples:", n_samples, "\n")
    cat("  New_time points summary:\n")
    print(summary(time_info))
    cat("  Is increasing:", all(diff(time_info) >= 0), "\n")
  }
}

cat("\nAM group bidirectional analysis completed with New_time ordering!\n")

################################################################################
# Results visualization and summary ----

# Dimension extraction function - extract CCM results by direction ----
extract_ccm_by_direction <- function(ccm_result, genus_names, env_names, reverse_env_factors) {
  # Get total dimensions of CCM matrix
  total_vars <- nrow(ccm_result$ccm.rho)
  
  # Number of genera and environmental factors
  n_genus <- length(genus_names)
  n_env <- length(env_names)
  
  # Verify dimension match
  if (n_genus + n_env != total_vars) {
    cat("Warning: Dimension mismatch!\n")
    cat("Total variables in CCM:", total_vars, "\n")
    cat("Genus count:", n_genus, "\n")
    cat("Environment count:", n_env, "\n")
    
    # Adjust dimensions: take the smaller value
    n_genus <- min(n_genus, total_vars)
    n_env <- total_vars - n_genus
    cat("Adjusted - Genus:", n_genus, "Environment:", n_env, "\n")
  }
  
  # Initialize result matrices
  ccm_rho_final <- matrix(NA, nrow = n_env, ncol = n_genus)
  ccm_sig_final <- matrix(0, nrow = n_env, ncol = n_genus)
  
  rownames(ccm_rho_final) <- env_names
  colnames(ccm_rho_final) <- genus_names
  rownames(ccm_sig_final) <- env_names
  colnames(ccm_sig_final) <- genus_names
  
  # Extract results
  for (env_idx in 1:n_env) {
    env_name <- env_names[env_idx]
    i <- n_genus + env_idx  # environmental factor row index in CCM matrix
    
    for (genus_idx in 1:n_genus) {
      genus_name <- genus_names[genus_idx]
      j <- genus_idx  # genus column index in CCM matrix
      
      if (env_name %in% reverse_env_factors) {
        # Reverse: Genus → Environmental factor (stored in matrix at genus row, environment column)
        ccm_rho_final[env_idx, genus_idx] <- ccm_result$ccm.rho[j, i]
        ccm_sig_final[env_idx, genus_idx] <- ccm_result$ccm.sig[j, i]
      } else {
        # Forward: Environmental factor → Genus (stored in matrix at environment row, genus column)
        ccm_rho_final[env_idx, genus_idx] <- ccm_result$ccm.rho[i, j]
        ccm_sig_final[env_idx, genus_idx] <- ccm_result$ccm.sig[i, j]
      }
    }
  }
  
  # Create a direction marker matrix (for color coding)
  direction_matrix <- matrix("Env→Genus", nrow = n_env, ncol = n_genus)
  rownames(direction_matrix) <- env_names
  colnames(direction_matrix) <- genus_names
  
  for (env_idx in 1:n_env) {
    env_name <- env_names[env_idx]
    if (env_name %in% reverse_env_factors) {
      direction_matrix[env_idx, ] <- "Genus→Env"  # genus points to environment
    }
  }
  
  # Create a color coding matrix
  color_matrix <- matrix("darkblue", nrow = n_env, ncol = n_genus)  # default blue: Env→Genus
  rownames(color_matrix) <- env_names
  colnames(color_matrix) <- genus_names
  
  for (env_idx in 1:n_env) {
    env_name <- env_names[env_idx]
    if (env_name %in% reverse_env_factors) {
      color_matrix[env_idx, ] <- "darkred"  # red: Genus→Env
    }
  }
  
  return(list(ccm_rho = ccm_rho_final, 
              ccm_sig = ccm_sig_final,
              direction = direction_matrix,
              direction_color = color_matrix))
}

################################################################################
# Visualization function: plot CCM heatmap with direction-colored significance boxes ----

# ~1. Plot correlation heatmap ----
plot_correlation_heatmap <- function(cor_matrix, cor_p_matrix, title, fixed_limits = TRUE) {
  
  # Convert data to long format
  cor_melted <- melt(cor_matrix)
  p_melted <- melt(cor_p_matrix)
  
  # Merge data
  plot_data <- data.frame(
    Genus = cor_melted$Var1,
    Environment = cor_melted$Var2,
    Correlation = cor_melted$value,
    p_value = p_melted$value
  )
  
  # Add significance markers
  plot_data$signif <- ifelse(plot_data$p_value < 0.05, "sig", "ns")
  
  # Create heatmap
  p <- ggplot(plot_data, aes(x = Environment, y = Genus, fill = Correlation)) +
    geom_tile(color = "white", linewidth = 0.5) +
    # Add significance boxes only on significant correlations
    geom_tile(data = subset(plot_data, signif == "sig"),
              aes(color = signif), fill = NA, linewidth = 1.2) +
    scale_fill_gradient2(
      low = "#3B7ABA", 
      mid = "white",
      high = "#C04C5A",
      midpoint = 0,
      limits = if(fixed_limits) c(-1, 1) else NULL,
      name = "Correlation (r)"
    ) +
    scale_color_manual(
      values = c("sig" = "black"),
      labels = c("sig" = "p < 0.05"),
      name = "Significance"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 30, face = "plain"),
      axis.text.y = element_text(size = 30, face = "plain"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_blank(),
      legend.title = element_text(
        size = 26, 
        face = "plain",
        margin = margin(b = 15)
      ),
      legend.text = element_text(size = 26),
      legend.key.height = unit(0.8, "cm"),
      legend.key.width = unit(0.8, "cm"),
      panel.grid = element_blank(),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
    ) +
    coord_fixed()
  
  return(p)
}

# ~2. Plot CCM heatmap with direction-colored significance boxes ----
plot_ccm_heatmap_with_direction_colors <- function(ccm_rho_matrix, ccm_sig_matrix, direction_color_matrix, title) {
  
  # Convert data to long format
  rho_melted <- melt(ccm_rho_matrix)
  sig_melted <- melt(ccm_sig_matrix)
  color_melted <- melt(direction_color_matrix)
  
  # Merge data
  plot_data <- data.frame(
    Genus = rho_melted$Var2,  # note: in matrix, environment is row, genus is column
    Environment = rho_melted$Var1,
    CCM_rho = rho_melted$value,
    CCM_sig = sig_melted$value,
    Direction_Color = color_melted$value  # direction color coding
  )
  
  # Ensure CCM rho values are in 0-1 range
  plot_data$CCM_rho <- ifelse(plot_data$CCM_rho < 0, 0, plot_data$CCM_rho)
  plot_data$CCM_rho <- ifelse(plot_data$CCM_rho > 1, 1, plot_data$CCM_rho)
  
  # Add significance direction
  plot_data$signif_direction <- ifelse(
    plot_data$CCM_sig == 1, 
    plot_data$Direction_Color,  # if significant, use direction color
    "not_sig"                   # if not significant
  )
  
  # Create heatmap
  p <- ggplot(plot_data, aes(x = Environment, y = Genus, fill = CCM_rho)) +
    geom_tile(color = "white", linewidth = 0.5) +
    # Add significance boxes with direction-dependent colors
    geom_tile(data = subset(plot_data, signif_direction %in% c( "darkred", "darkblue" )),
              aes(color = signif_direction), fill = NA, linewidth = 1.5) +
    scale_fill_gradient(
      low = "#FFF8E1", 
      high = "#D2691E",
      limits = c(0, 1),
      name = "CCM (ρ)"
    ) +
    scale_color_manual(
      values = c(
        "darkred" = "#C04C5A",       # red: Genus → Env
        "darkblue" = "#3B7ABA",      # blue: Env → Genus
        "not_sig" = "transparent"    # not significant: transparent (not shown)
      ),
      breaks = c("darkred", "darkblue"), 
      labels = c(
        "darkred" = "Genus → Env ", 
        "darkblue" = "Env → Genus "
      ),
      name = "Causality Direction",
      guide = guide_legend(
        override.aes = list(fill = NA, color = c("#C04C5A", "#3B7ABA" )),  
        keywidth = unit(1.5, "cm"),
        keyheight = unit(0.8, "cm")
      )
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 30, face = "plain"),
      axis.text.y = element_text(size = 30, face = "plain"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_blank(),
      legend.title = element_text(
        size = 26, 
        face = "plain",
        margin = margin(b = 15)
      ),
      legend.text = element_text(size = 26),
      legend.key.height = unit(0.8, "cm"),
      legend.key.width = unit(0.8, "cm"),
      panel.grid = element_blank(),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
    ) +
    coord_fixed()
  
  return(p)
}

# ~3. Function to combine two heatmaps ----
plot_combined_separated_heatmaps_bidirectional_colors <- function(cor_matrix, cor_p_matrix, 
                                                                  ccm_rho_matrix, ccm_sig_matrix, 
                                                                  direction_color_matrix, title_suffix, 
                                                                  output_dir, panel_label = NULL) {
  
  # Plot correlation heatmap
  p_corr <- plot_correlation_heatmap(
    cor_matrix = cor_matrix,
    cor_p_matrix = cor_p_matrix,
    title = title_suffix,
    fixed_limits = TRUE
  )
  
  # Plot CCM heatmap with direction color markers
  p_ccm <- plot_ccm_heatmap_with_direction_colors(
    ccm_rho_matrix = ccm_rho_matrix,
    ccm_sig_matrix = ccm_sig_matrix,
    direction_color_matrix = direction_color_matrix,
    title = title_suffix
  )
  
  # Combine two heatmaps
  p_combined <- p_corr + p_ccm +
    plot_layout(ncol = 2, widths = c(1, 1), guides = "collect") +
    plot_annotation(
      title = paste("Cd ", title_suffix),
      theme = theme(
        plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
        plot.margin = margin(1, 1, 1, 1, "cm")
      )
    )
  
  # If panel label provided, add to top-left corner
  if (!is.null(panel_label)) {
    library(grid)
    p_combined <- p_combined +
      annotation_custom(
        grob = grid::textGrob(
          label = panel_label,
          hjust = 0, vjust = 0,
          gp = grid::gpar(
            fontsize = 40,
            fontface = "bold",
            col = "black"
          )
        ),
        xmin = -Inf, xmax = Inf,
        ymin = Inf, ymax = Inf
      )
  }
  
  # Save combined heatmap
  output_path <- file.path(output_dir, paste0("bidirectional_color_heatmaps_Cd_", title_suffix, ".png"))
  ggsave(output_path, p_combined, width = 22, height = 9.8, dpi = 300)
  
  output_path2 <- file.path(output_dir, paste0("bidirectional_color_heatmaps_Cd_", title_suffix, ".pdf"))
  ggsave(output_path2, p_combined, width = 22, height = 9.8, dpi = 300)
  
  # Also save each heatmap separately
  ggsave(file.path(output_dir, paste0("correlation_heatmap_Cd_", title_suffix, ".png")),
         p_corr, width = 12, height = 8, dpi = 300)
  ggsave(file.path(output_dir, paste0("ccm_bidirectional_color_heatmap_Cd_", title_suffix, ".png")),
         p_ccm, width = 12, height = 8, dpi = 300)
  ggsave(file.path(output_dir, paste0("correlation_heatmap_Cd_", title_suffix, ".pdf")),
         p_corr, width = 12, height = 8, dpi = 300)
  ggsave(file.path(output_dir, paste0("ccm_bidirectional_color_heatmap_Cd_", title_suffix, ".pdf")),
         p_ccm, width = 12, height = 8, dpi = 300)
  
  cat("  Saved bidirectional color heatmaps for Cd =", title_suffix, "\n")
  
  return(list(correlation_plot = p_corr, ccm_plot = p_ccm, combined_plot = p_combined))
}

################################################################################
# Generate separate heatmaps for each Cd group and combine ----
cat("\n" , strrep("*", 60), "\n")
cat("Generating bidirectional color heatmaps for AM group results\n")
cat(strrep("*", 60), "\n\n")

# Store all plots
all_plots <- list()
all_combined_plots <- list()

# Define panel labels
panel_labels <- c("A", "B", "C", "D")
cd_labels <- data.frame(
  cd = cd_gradients,
  label = panel_labels[1:length(cd_gradients)]
)

for (i in 1:length(cd_gradients)) {
  cd <- cd_gradients[i]
  cd_name <- paste0("Cd_", cd)
  
  if (cd_name %in% names(results_am$correlation)) {
    
    cat("Processing Cd =", cd, "...\n")
    
    # Get current Cd group data
    cor_matrix <- results_am$correlation[[cd_name]]$correlation
    cor_p_matrix <- results_am$correlation[[cd_name]]$p_value
    ccm_result <- results_am$causality[[cd_name]]
    
    # Ensure genus and environment factor order is correct
    genus_names <- target_genera_ordered[target_genera_ordered %in% rownames(cor_matrix)]
    env_names <- env_order_final[env_order_final %in% colnames(cor_matrix)]
    
    # Reorder matrices
    cor_matrix <- cor_matrix[genus_names, env_names, drop = FALSE]
    cor_p_matrix <- cor_p_matrix[genus_names, env_names, drop = FALSE]
    
    # Use new function to extract CCM results (including direction and color coding)
    ccm_subset <- extract_ccm_by_direction(ccm_result, genus_names, env_names, reverse_env_factors)
    
    # Get corresponding panel label
    panel_label <- cd_labels$label[cd_labels$cd == cd]
    
    # Plot separate heatmaps and combine (add panel label)
    plots <- plot_combined_separated_heatmaps_bidirectional_colors(
      cor_matrix = cor_matrix,
      cor_p_matrix = cor_p_matrix,
      ccm_rho_matrix = ccm_subset$ccm_rho,
      ccm_sig_matrix = ccm_subset$ccm_sig,
      direction_color_matrix = ccm_subset$direction_color,
      title_suffix = cd,
      output_dir = "out/microbe/causality_heatmap",
      panel_label = panel_label
    )
    
    # Store plots
    all_plots[[cd_name]] <- plots
    all_combined_plots[[cd_name]] <- plots$combined_plot
    
    cat("  Generated bidirectional color heatmaps for Cd =", cd, "with label", panel_label, "\n")
    cat("  Color coding: Blue = Env → Genus, Red = Genus → Env\n")
  }
}

################################################################################
# Combine all Cd group plots into one figure ----
cat("\n" , strrep("*", 60), "\n")
cat("Combining all Cd gradient plots into one figure\n")
cat(strrep("*", 60), "\n\n")

if (length(all_combined_plots) > 0) {
  
  # Extract all combined plots
  combined_plot_list <- lapply(names(all_combined_plots), function(cd_name) {
    return(all_combined_plots[[cd_name]])
  })
  
  # Calculate layout
  n_plots <- length(combined_plot_list)
  
  if (n_plots == 4) {
    # 4 plots: 2 rows x 2 columns
    layout_matrix <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2, byrow = TRUE)
  } else if (n_plots == 3) {
    # 3 plots: 3 rows x 1 column
    layout_matrix <- matrix(c(1, 2, 3), nrow = 3, ncol = 1, byrow = TRUE)
  } else if (n_plots == 2) {
    # 2 plots: 2 rows x 1 column
    layout_matrix <- matrix(c(1, 2), nrow = 2, ncol = 1, byrow = TRUE)
  } else {
    # 1 plot or more: use default layout
    layout_matrix <- NULL
  }
  
  # Create final figure
  if (n_plots == 4) {
    # Use wrap_plots with specified layout
    p_final <- wrap_plots(
      combined_plot_list,
      nrow = 2,
      ncol = 2,
      guides = "collect"
    ) +
      plot_annotation(
        title = "",
        theme = theme(
          plot.title = element_text(
            size = 35,
            face = "bold",
            hjust = 0.5,
            margin = margin(b = 20)
          ),
          plot.margin = margin(2, 2, 2, 2, "cm")
        )
      )
    
    # Save final figure
    output_path_final <- "out/microbe/causality_heatmap/all_Cd_gradients_bidirectional_color_combined.png"
    ggsave(
      output_path_final,
      p_final,
      width = 45,
      height = 20,
      dpi = 300,
      bg = "white"
    )
    
    output_path_final2 <- "out/microbe/causality_heatmap/all_Cd_gradients_bidirectional_color_combined.pdf"
    ggsave(
      output_path_final2,
      p_final,
      width = 45,
      height = 20,
      dpi = 300,
      bg = "white"
    )
    
    cat("  Combined figure saved: ", output_path_final, "\n")
    cat("  Figure dimensions: 45 x 20 inches\n")
    
  } else {
    # For other numbers of plots, use more flexible layout
    p_final <- wrap_plots(
      combined_plot_list,
      guides = "collect"
    ) +
      plot_annotation(
        title = "",
        theme = theme(
          plot.title = element_text(
            size = 35,
            face = "bold",
            hjust = 0.5,
            margin = margin(b = 20)
          ),
          plot.margin = margin(2, 2, 2, 2, "cm")
        )
      )
    
    # Adjust size based on number of plots
    if (n_plots == 3) {
      fig_width <- 35
      fig_height <- 30
    } else if (n_plots == 2) {
      fig_width <- 25
      fig_height <- 18
    } else {
      fig_width <- 25
      fig_height <- 10
    }
    
    # Save final figure
    output_path_final <- "out/microbe/causality_heatmap/all_Cd_gradients_bidirectional_color_combined.png"
    ggsave(
      output_path_final,
      p_final,
      width = fig_width,
      height = fig_height,
      dpi = 300,
      bg = "white"
    )
    
    cat("  Combined figure saved: ", output_path_final, "\n")
    cat("  Figure dimensions: ", fig_width, "x", fig_height, "inches\n")
  }
}


# Save workspace
save.image(file = "out/microbe/causality_heatmap/am_genera_env_ccm_corr.RData")