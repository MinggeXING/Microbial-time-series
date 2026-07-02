# Load required packages
library(tidyverse)
library(pheatmap)
library(ggplot2)
library(patchwork) 
library(ggpubr)
library(grid)


# Plot genus heatmap for stage data ----
################################################################################
# Read data
genus_stage_abs_target <- read.csv("data/stage_data/target_genera_abs_stage.csv", header = TRUE)
stage_metadata <- read.csv("data/stage_data/stage_metadata.csv", header = TRUE)

genus_stage_abs_meta <- genus_stage_abs_target %>%
  left_join(stage_metadata, by = "Sample_ID")

# Filter raw abundance columns ----
raw_genus_cols <- names(genus_stage_abs_target)[!str_detect(names(genus_stage_abs_target), "^log_")]
raw_genus_cols <- raw_genus_cols[raw_genus_cols != "Sample_ID"]

# Extract raw abundance data
genus_stage_raw <- genus_stage_abs_meta %>%
  select(Sample_ID, all_of(raw_genus_cols), Cd, Mycorrhizal, stage)

# Check data structure
cat("Total samples:", nrow(genus_stage_raw), "\n")
cat("AM group samples:", sum(genus_stage_raw$Mycorrhizal == "AM"), "\n")
cat("NM group samples:", sum(genus_stage_raw$Mycorrhizal == "NM"), "\n")
cat("Cd gradients:", paste(sort(unique(genus_stage_raw$Cd)), collapse = ", "), "\n")
cat("Stages:", paste(sort(unique(genus_stage_raw$stage)), collapse = ", "), "\n")

# Filter AM and NM groups
am_stage_data <- genus_stage_raw %>%
  filter(Mycorrhizal == "AM")

nm_stage_data <- genus_stage_raw %>%
  filter(Mycorrhizal == "NM")

# Check sample sizes after filtering
cat("\nSample sizes after filtering:\n")
cat("AM group samples:", nrow(am_stage_data), "\n")
cat("NM group samples:", nrow(nm_stage_data), "\n")
cat("AM group stages:", paste(sort(unique(am_stage_data$stage)), collapse = ", "), "\n")
cat("NM group stages:", paste(sort(unique(nm_stage_data$stage)), collapse = ", "), "\n")

# Get all Cd gradients
cd_levels <- unique(genus_stage_abs_meta$Cd)
cd_levels <- cd_levels[!is.na(cd_levels)]
cat("\nAll Cd gradients:", paste(cd_levels, collapse = ", "), "\n")

# Create output directory
stage_output_dir <- "out/microbe/abs_heatmap_stage"
if (!dir.exists(stage_output_dir)) {
  dir.create(stage_output_dir, recursive = TRUE)
}

# Define function to plot heatmap (stage version) ----
plot_heatmap_by_cd_stage <- function(cd_value, genus_order = NULL) {
  
  # Filter data for current Cd gradient
  am_cd <- am_stage_data %>% filter(Cd == cd_value)
  nm_cd <- nm_stage_data %>% filter(Cd == cd_value)
  
  cat("\n=== Cd", cd_value, "gradient analysis ===\n")
  cat("AM group samples:", nrow(am_cd), "\n")
  cat("NM group samples:", nrow(nm_cd), "\n")
  
  # Get common stages
  common_stages <- intersect(am_cd$stage, nm_cd$stage)
  cat("Common stages:", paste(sort(common_stages), collapse = ", "), "\n")
  
  # Initialize result matrices
  diff_matrix <- matrix(NA, nrow = length(common_stages), ncol = length(raw_genus_cols))
  pvalue_matrix <- matrix(NA, nrow = length(common_stages), ncol = length(raw_genus_cols))
  rownames(diff_matrix) <- as.character(sort(common_stages))
  rownames(pvalue_matrix) <- as.character(sort(common_stages))
  colnames(diff_matrix) <- raw_genus_cols
  colnames(pvalue_matrix) <- raw_genus_cols
  
  # Calculate AM-NM difference and t-test p-value for each genus at each stage
  for(stage_point in common_stages) {
    # Get AM and NM data for current stage
    am_stage <- am_cd %>% filter(stage == stage_point)
    nm_stage <- nm_cd %>% filter(stage == stage_point)
    
    # Calculate mean for each genus
    am_means <- colMeans(am_stage[, raw_genus_cols], na.rm = TRUE)
    nm_means <- colMeans(nm_stage[, raw_genus_cols], na.rm = TRUE)
    
    # Calculate difference
    diff_values <- am_means - nm_means
    
    # Store in matrix
    diff_matrix[as.character(stage_point), ] <- diff_values
    
    # Perform t-test for each genus
    for(genus in raw_genus_cols) {
      am_values <- am_stage[[genus]]
      nm_values <- nm_stage[[genus]]
      
      # Check if there is enough data for t-test
      if(length(am_values) >= 2 && length(nm_values) >= 2 && 
         !all(is.na(am_values)) && !all(is.na(nm_values))) {
        # Perform t-test
        t_test_result <- try(t.test(am_values, nm_values, var.equal = TRUE), silent = TRUE)
        if(!inherits(t_test_result, "try-error")) {
          pvalue_matrix[as.character(stage_point), genus] <- t_test_result$p.value
        } else {
          pvalue_matrix[as.character(stage_point), genus] <- NA
        }
      } else {
        pvalue_matrix[as.character(stage_point), genus] <- NA
      }
    }
  }
  
  # Check for NA values
  if(any(is.na(diff_matrix))) {
    cat("Warning: NA values exist, replacing with 0\n")
    diff_matrix[is.na(diff_matrix)] <- 0
  }
  
  # Z-score normalization of the difference matrix (by column)
  zscore_matrix <- apply(diff_matrix, 2, function(x) {
    if(sd(x) == 0) {
      return(rep(0, length(x)))  # If standard deviation is 0, return vector of 0s
    } else {
      return(scale(x))
    }
  })
  rownames(zscore_matrix) <- rownames(diff_matrix)
  
  cat("Normalized matrix dimensions:", dim(zscore_matrix), "\n")
  cat("Normalized data range:", round(range(zscore_matrix, na.rm = TRUE), 3), "\n")
  
  # Convert matrix to long format data frame
  heatmap_data <- as.data.frame(zscore_matrix) %>%
    rownames_to_column("stage_point") %>%
    pivot_longer(cols = -stage_point, names_to = "Genus", values_to = "degree")
  
  # Convert p-value matrix to long format and add significance markers
  pvalue_data <- as.data.frame(pvalue_matrix) %>%
    rownames_to_column("stage_point") %>%
    pivot_longer(cols = -stage_point, names_to = "Genus", values_to = "pvalue")
  
  # Add significance markers
  pvalue_data <- pvalue_data %>%
    mutate(
      significance = case_when(
        pvalue < 0.001 ~ "***",
        pvalue < 0.01 ~ "**",
        pvalue < 0.05 ~ "*",
        TRUE ~ ""
      )
    )
  
  # Merge significance markers into heatmap data
  heatmap_data <- heatmap_data %>%
    left_join(pvalue_data, by = c("stage_point", "Genus"))
  
  # Calculate stage sequence features for ordering
  # Convert stage to ordered factor, maintain stage order
  stage_levels <- as.character(sort(common_stages))
  
  # Calculate stage sequence features for each genus
  genus_patterns <- data.frame(Genus = colnames(zscore_matrix))
  
  for(genus in colnames(zscore_matrix)) {
    values <- zscore_matrix[, genus]
    
    # Calculate stage trend (using stage order rather than numeric values)
    # Assign numeric weights to stages (in order of appearance)
    stage_weights <- 1:length(stage_levels)
    
    # Calculate trend (linear regression slope)
    if(length(unique(stage_weights)) > 1) {
      trend_model <- lm(values ~ stage_weights)
      genus_patterns[genus_patterns$Genus == genus, "trend"] <- coef(trend_model)[2]
    } else {
      genus_patterns[genus_patterns$Genus == genus, "trend"] <- 0
    }
    
    # Calculate early and late period means
    early_period <- 1:floor(length(stage_levels)/2)  # first half stages
    late_period <- (length(stage_levels) - floor(length(stage_levels)/2) + 1):length(stage_levels)  # second half stages
    
    genus_patterns[genus_patterns$Genus == genus, "early_mean"] <- mean(values[early_period], na.rm = TRUE)
    genus_patterns[genus_patterns$Genus == genus, "late_mean"] <- mean(values[late_period], na.rm = TRUE)
    
    # Calculate change pattern: early high vs early low
    genus_patterns[genus_patterns$Genus == genus, "pattern_score"] <- 
      genus_patterns[genus_patterns$Genus == genus, "early_mean"] - 
      genus_patterns[genus_patterns$Genus == genus, "late_mean"]
  }
  
  # Order genera based on change pattern
  # If genus_order is provided, use that order
  # Otherwise, sort by pattern_score descending
  if(is.null(genus_order)) {
    genus_patterns <- genus_patterns[order(-genus_patterns$pattern_score), ]
    genus_order <- genus_patterns$Genus
  } else {
    # Use provided genus order
    genus_order <- genus_order
    cat("Using provided genus order\n")
  }
  
  # Set genus order
  heatmap_data$Genus <- factor(heatmap_data$Genus, levels = genus_order)
  
  # Set stage order reversed (stage1 at top, stage5 at bottom)
  heatmap_data$stage_point <- factor(heatmap_data$stage_point, levels = rev(stage_levels))
  
  # Calculate color midpoint (based on data range)
  data_range <- range(heatmap_data$degree, na.rm = TRUE)
  midpoint_value <- (data_range[1] + data_range[2]) / 2
  
  # Calculate significance test statistics
  total_tests <- sum(!is.na(pvalue_data$pvalue))
  significant_tests <- sum(pvalue_data$pvalue < 0.05, na.rm = TRUE)
  cat("Total t-tests performed:", total_tests, "\n")
  cat("Significant results (p < 0.05):", significant_tests, "\n")
  
  # Plot heatmap
  heatmap_plot <- ggplot(heatmap_data, aes(x = Genus, 
                                           y = stage_point, 
                                           fill = degree)) +
    geom_tile(color = "white", size = 0.5) +
    # Add significance markers
    geom_text(aes(label = significance), 
              size = 18,  
              color = "black",
              fontface = "bold",
              vjust = 0.8) +
    scale_fill_gradient2(
      low = "#254751",
      mid = "white",
      high = "#DB6951",
      midpoint = midpoint_value,
      name = "Abundance\n(AM - NM)",
      guide = guide_colorbar(
        barwidth = unit(1.5, "cm"),
        barheight = unit(8, "cm")
      )
    ) +
    labs(
      title = paste("Cd ", cd_value),
      x = "",
      y = ""
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 55, face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 39),
      axis.text.y = element_text(size = 38),
      axis.title.y  = element_text(size = 39),
      legend.position = "right",
      panel.grid = element_blank(),
      legend.title = element_text(size = 37, face = "plain"),
      legend.text = element_text(size = 34),
      legend.box.background = element_blank()
    )
  
  # Save images
  png_file <- file.path(stage_output_dir, paste0("abs_heatmap_stage_cd", cd_value, ".png"))
  pdf_file <- file.path(stage_output_dir, paste0("abs_heatmap_stage_cd", cd_value, ".pdf"))
  
  
  ggsave(png_file, heatmap_plot, width = 26, height = 14, dpi = 300)
  ggsave(pdf_file, heatmap_plot, width = 26, height = 14)
  
  
  print(heatmap_plot)
  
  # Return results for inspection
  return(list(
    diff_matrix = diff_matrix,
    zscore_matrix = zscore_matrix,
    pvalue_matrix = pvalue_matrix,
    heatmap_data = heatmap_data,
    common_stages = common_stages,
    genus_patterns = genus_patterns,
    genus_order = genus_order, 
    plot = heatmap_plot
  ))
}

# First run Cd 0 to get genus order
cat("\n=== Running Cd 0 to get genus order ===\n")
cd0_result <- plot_heatmap_by_cd_stage(cd_value = 0)
cd0_genus_order <- cd0_result$genus_order
results_stage_ggplot[["Cd0"]] <- cd0_result
plot_list[["Cd0"]] <- cd0_result$plot

# Then run other Cd gradients using Cd 0's genus order
other_cd_levels <- cd_levels[cd_levels != 0]

for(cd_val in other_cd_levels) {
  cat("\n=== Running Cd", cd_val, "using Cd 0's genus order ===\n")
  result <- plot_heatmap_by_cd_stage(
    cd_value = cd_val, 
    genus_order = cd0_genus_order
  )
  results_stage_ggplot[[paste0("Cd", cd_val)]] <- result
  plot_list[[paste0("Cd", cd_val)]] <- result$plot
}

################################################################################
# Define function to plot combined heatmap figure ----
plot_heatmap_cd_stage_combined <- function(cd_value, genus_order = NULL, show_axis_text = TRUE) {
  
  # Filter data for current Cd gradient
  am_cd <- am_stage_data %>% filter(Cd == cd_value)
  nm_cd <- nm_stage_data %>% filter(Cd == cd_value)
  
  cat("\n=== Cd", cd_value, "gradient analysis ===\n")
  cat("AM group samples:", nrow(am_cd), "\n")
  cat("NM group samples:", nrow(nm_cd), "\n")
  
  # Get common stages
  common_stages <- intersect(am_cd$stage, nm_cd$stage)
  cat("Common stages:", paste(sort(common_stages), collapse = ", "), "\n")
  
  # Initialize result matrices
  diff_matrix <- matrix(NA, nrow = length(common_stages), ncol = length(raw_genus_cols))
  pvalue_matrix <- matrix(NA, nrow = length(common_stages), ncol = length(raw_genus_cols))
  rownames(diff_matrix) <- as.character(sort(common_stages))
  rownames(pvalue_matrix) <- as.character(sort(common_stages))
  colnames(diff_matrix) <- raw_genus_cols
  colnames(pvalue_matrix) <- raw_genus_cols
  
  # Calculate AM-NM difference and t-test p-value for each genus at each stage
  for(stage_point in common_stages) {
    # Get AM and NM data for current stage
    am_stage <- am_cd %>% filter(stage == stage_point)
    nm_stage <- nm_cd %>% filter(stage == stage_point)
    
    # Calculate mean for each genus
    am_means <- colMeans(am_stage[, raw_genus_cols], na.rm = TRUE)
    nm_means <- colMeans(nm_stage[, raw_genus_cols], na.rm = TRUE)
    
    # Calculate difference
    diff_values <- am_means - nm_means
    
    # Store in matrix
    diff_matrix[as.character(stage_point), ] <- diff_values
    
    # Perform t-test for each genus
    for(genus in raw_genus_cols) {
      am_values <- am_stage[[genus]]
      nm_values <- nm_stage[[genus]]
      
      # Check if there is enough data for t-test
      if(length(am_values) >= 2 && length(nm_values) >= 2 && 
         !all(is.na(am_values)) && !all(is.na(nm_values))) {
        # Perform t-test
        t_test_result <- try(t.test(am_values, nm_values, var.equal = TRUE), silent = TRUE)
        if(!inherits(t_test_result, "try-error")) {
          pvalue_matrix[as.character(stage_point), genus] <- t_test_result$p.value
        } else {
          pvalue_matrix[as.character(stage_point), genus] <- NA
        }
      } else {
        pvalue_matrix[as.character(stage_point), genus] <- NA
      }
    }
  }
  
  # Check for NA values
  if(any(is.na(diff_matrix))) {
    cat("Warning: NA values exist, replacing with 0\n")
    diff_matrix[is.na(diff_matrix)] <- 0
  }
  
  # Z-score normalization of the difference matrix (by column)
  zscore_matrix <- apply(diff_matrix, 2, function(x) {
    if(sd(x) == 0) {
      return(rep(0, length(x)))  # If standard deviation is 0, return vector of 0s
    } else {
      return(scale(x))
    }
  })
  rownames(zscore_matrix) <- rownames(diff_matrix)
  
  cat("Normalized matrix dimensions:", dim(zscore_matrix), "\n")
  cat("Normalized data range:", round(range(zscore_matrix, na.rm = TRUE), 3), "\n")
  
  # Convert matrix to long format data frame
  heatmap_data <- as.data.frame(zscore_matrix) %>%
    rownames_to_column("stage_point") %>%
    pivot_longer(cols = -stage_point, names_to = "Genus", values_to = "degree")
  
  # Convert p-value matrix to long format and add significance markers
  pvalue_data <- as.data.frame(pvalue_matrix) %>%
    rownames_to_column("stage_point") %>%
    pivot_longer(cols = -stage_point, names_to = "Genus", values_to = "pvalue")
  
  # Add significance markers
  pvalue_data <- pvalue_data %>%
    mutate(
      significance = case_when(
        pvalue < 0.001 ~ "***",
        pvalue < 0.01 ~ "**",
        pvalue < 0.05 ~ "*",
        TRUE ~ ""
      )
    )
  
  # Merge significance markers into heatmap data
  heatmap_data <- heatmap_data %>%
    left_join(pvalue_data, by = c("stage_point", "Genus"))
  
  # Calculate stage sequence features for ordering
  # Convert stage to ordered factor, maintain stage order
  stage_levels <- as.character(sort(common_stages))
  
  # Calculate stage sequence features for each genus
  genus_patterns <- data.frame(Genus = colnames(zscore_matrix))
  
  for(genus in colnames(zscore_matrix)) {
    values <- zscore_matrix[, genus]
    
    # Calculate stage trend (using stage order rather than numeric values)
    # Assign numeric weights to stages (in order of appearance)
    stage_weights <- 1:length(stage_levels)
    
    # Calculate trend (linear regression slope)
    if(length(unique(stage_weights)) > 1) {
      trend_model <- lm(values ~ stage_weights)
      genus_patterns[genus_patterns$Genus == genus, "trend"] <- coef(trend_model)[2]
    } else {
      genus_patterns[genus_patterns$Genus == genus, "trend"] <- 0
    }
    
    # Calculate early and late period means
    early_period <- 1:floor(length(stage_levels)/2)  # first half stages
    late_period <- (length(stage_levels) - floor(length(stage_levels)/2) + 1):length(stage_levels)  # second half stages
    
    genus_patterns[genus_patterns$Genus == genus, "early_mean"] <- mean(values[early_period], na.rm = TRUE)
    genus_patterns[genus_patterns$Genus == genus, "late_mean"] <- mean(values[late_period], na.rm = TRUE)
    
    # Calculate change pattern: early high vs early low
    genus_patterns[genus_patterns$Genus == genus, "pattern_score"] <- 
      genus_patterns[genus_patterns$Genus == genus, "early_mean"] - 
      genus_patterns[genus_patterns$Genus == genus, "late_mean"]
  }
  
  # Order genera based on change pattern
  # If genus_order is provided, use that order
  # Otherwise, sort by pattern_score descending
  if(is.null(genus_order)) {
    genus_patterns <- genus_patterns[order(-genus_patterns$pattern_score), ]
    genus_order <- genus_patterns$Genus
  } else {
    # Use provided genus order
    genus_order <- genus_order
    cat("Using provided genus order\n")
  }
  
  # Set genus order
  heatmap_data$Genus <- factor(heatmap_data$Genus, levels = genus_order)
  
  # Set stage order reversed (stage1 at top, stage5 at bottom)
  heatmap_data$stage_point <- factor(heatmap_data$stage_point, levels = rev(stage_levels))
  
  # Calculate color midpoint (based on data range)
  data_range <- range(heatmap_data$degree, na.rm = TRUE)
  midpoint_value <- (data_range[1] + data_range[2]) / 2
  
  # Calculate significance test statistics
  total_tests <- sum(!is.na(pvalue_data$pvalue))
  significant_tests <- sum(pvalue_data$pvalue < 0.05, na.rm = TRUE)
  cat("Total t-tests performed:", total_tests, "\n")
  cat("Significant results (p < 0.05):", significant_tests, "\n")
  
  # Plot heatmap
  heatmap_plot <- ggplot(heatmap_data, aes(x = Genus, 
                                           y = stage_point, 
                                           fill = degree)) +
    geom_tile(color = "white", size = 0.5) +
    # Add significance markers
    geom_text(aes(label = significance), 
              size = 20,  
              color = "black",
              fontface = "bold",
              vjust = 0.8) +
    scale_fill_gradient2(
      low = "#254751",
      mid = "white",
      high = "#DB6951",
      midpoint = midpoint_value,
      name = "Abundance\n(AM - NM)"
    ) +
    labs(
      title = paste("Cd ", cd_value),
      x = "",
      y = ""
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 75, face = "plain", hjust = 0.5),
      axis.text.x = element_blank(),  
      axis.text.y = element_blank(),  
      axis.ticks = element_blank(),   
      axis.title = element_blank(),   
      legend.position = "none",      
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", size = 1),  
      plot.margin = margin(5, 5, 5, 5, "pt")
    )
  
  # Return results for inspection
  return(list(
    diff_matrix = diff_matrix,
    zscore_matrix = zscore_matrix,
    pvalue_matrix = pvalue_matrix,
    heatmap_data = heatmap_data,
    common_stages = common_stages,
    genus_patterns = genus_patterns,
    genus_order = genus_order,  
    plot = heatmap_plot
  ))
}


# Run function for each Cd gradient and collect plot objects ----
results_stage_ggplot <- list()
plot_list <- list()  # Store plot objects

# First run Cd 0 to get genus order
cat("\n=== Running Cd 0 to get genus order ===\n")
cd0_result <- plot_heatmap_cd_stage_combined(cd_value = 0)
cd0_genus_order <- cd0_result$genus_order
results_stage_ggplot[["Cd0"]] <- cd0_result
plot_list[["Cd0"]] <- cd0_result$plot

# Then run other Cd gradients using Cd 0's genus order
other_cd_levels <- cd_levels[cd_levels != 0]

for(cd_val in other_cd_levels) {
  cat("\n=== Running Cd", cd_val, "using Cd 0's genus order ===\n")
  result <- plot_heatmap_cd_stage_combined(
    cd_value = cd_val, 
    genus_order = cd0_genus_order
  )
  results_stage_ggplot[[paste0("Cd", cd_val)]] <- result
  plot_list[[paste0("Cd", cd_val)]] <- result$plot
}

# Add labels to each plot
plot_list$Cd0 <- plot_list$Cd0 + 
  labs(tag = "A") +
  theme(plot.tag = element_text(size = 75, face = "bold", hjust = 0, vjust = 0))

plot_list$Cd2 <- plot_list$Cd2 + 
  labs(tag = "B") +
  theme(plot.tag = element_text(size = 75, face = "bold", hjust = 0, vjust = 0))

plot_list$Cd5 <- plot_list$Cd5 + 
  labs(tag = "C") +
  theme(plot.tag = element_text(size = 75, face = "bold", hjust = 0, vjust = 0))

plot_list$Cd15 <- plot_list$Cd15 + 
  labs(tag = "D") +
  theme(plot.tag = element_text(size = 75, face = "bold", hjust = 0, vjust = 0))

# Modify heatmap themes to add axis labels for different positions
# Cd0 (top-left): show left y-axis labels
plot_list$Cd0 <- plot_list$Cd0 + 
  theme(
    axis.text.y = element_text(size = 75, color = "black"),
    axis.ticks.y = element_line(color = "black")
  )

# Cd2 (top-right): hide y-axis labels
plot_list$Cd2 <- plot_list$Cd2 + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  )

# Cd5 (bottom-left): show left y-axis labels and bottom x-axis labels
plot_list$Cd5 <- plot_list$Cd5 + 
  theme(
    axis.text.y = element_text(size = 75, color = "black"),
    axis.ticks.y = element_line(color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 75, color = "black"),
    axis.ticks.x = element_line(color = "black")
  )

# Cd15 (bottom-right): show bottom x-axis labels
plot_list$Cd15 <- plot_list$Cd15 + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 75, color = "black"),
    axis.ticks.x = element_line(color = "black")
  )

# Create shared legend
# Use one of the plots to extract legend, reduce legend size
legend_plot <- plot_list$Cd0 + 
  theme(legend.position = "right",
        legend.title = element_text(size = 75, face = "plain"),  
        legend.text = element_text(size = 75),  
        legend.key.height = unit(6, "cm"),  
        legend.key.width = unit(3, "cm"))  

# Extract legend
legend <- get_legend(legend_plot)

# Create ggplot object for legend
legend_ggplot <- ggplot() + 
  theme_void() + 
  annotation_custom(legend, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
  theme(plot.margin = margin(0, 0, 0, 0))

# Combine 2x2 heatmap grid using patchwork
heatmap_grid <- (plot_list$Cd0 | plot_list$Cd2) / (plot_list$Cd5 | plot_list$Cd15)

final_plot <- wrap_plots(
  heatmap_grid,
  legend_ggplot,
  ncol = 2,
  widths = c(8, 1)  
)

# Save combined figure
png_file_combined <- file.path(stage_output_dir, "abs_heatmap_stage_combined_2x2.png")
pdf_file_combined <- file.path(stage_output_dir, "abs_heatmap_stage_combined_2x2.pdf")


ggsave(png_file_combined, final_plot, width = 58, height = 32, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(pdf_file_combined, final_plot, width = 58, height = 32, limitsize = FALSE)


print(final_plot)

# Save data for further analysis
save.image(file = paste0(stage_output_dir, "/abs_heatmap_stage_results.RData"))