# Beta diversity analysis for all samples - three-factor version
my_packages <- c('phyloseq', 'tidyverse', 'ggplot2', 'ggpubr', 'igraph', 'vegan', 
                 'tibble', 'gridExtra', 'grid')

lapply(my_packages, library, character.only = TRUE) # Load required packages

# Beta diversity analysis by Cd gradient and stage ----
# Load required packages
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(vegan)
library(gridExtra)
library(grid)

# Define function: analyze and plot for a single Cd level and stage combination
create_nmds_plot <- function(cd_level, stage_level, pl_rarefied) {
  subset_samples <- subset_samples(pl_rarefied, cd == cd_level & stage == stage_level)
  if (nsamples(subset_samples) < 3) {
    p <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "Insufficient data", size = 6) +
      theme_void() +
      labs(title = paste0("Cd", cd_level, "-Stage", stage_level))
    return(p)
  }
  
  # Calculate Bray-Curtis distance
  otu <- as(otu_table(subset_samples), "matrix")
  if(taxa_are_rows(subset_samples)) { 
    otu <- t(otu)
  }
  bray_dist <- vegdist(otu, method = "bray")
  
  metadata <- data.frame(sample_data(subset_samples))
  
  myco_counts <- table(metadata$treat)
  print(paste("Cd", cd_level, "Stage", stage_level, "Myco counts:"))
  print(myco_counts)
  
  if (length(myco_counts) < 2 || any(myco_counts < 2)) {
    p <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, 
               label = paste("Insufficient groups:\nAM:", myco_counts["AM"], "NM:", myco_counts["NM"]), 
               size = 5) +
      theme_void() +
      labs(title = paste0("Cd", cd_level, "-Stage", stage_level))
    return(p)
  }
  
  
  myco_result <- NULL
  r2_value <- NA
  p_value <- NA
  significance <- "ns"
  
  # Perform PERMANOVA test
  tryCatch({
    permanova <- adonis2(
      bray_dist ~ treat,
      data = metadata,
      permutations = 999
    )
    
    print(paste("Cd", cd_level, "Stage", stage_level, "PERMANOVA results:"))
    print(permanova)
    
    # Extract Myco effect results
    if ("treat" %in% rownames(permanova)) {
      myco_result <- as.data.frame(permanova)["treat", ]
    } else if (nrow(permanova) >= 2) {
      # If no "treat" row name, try to extract the first row (excluding Total and Residual)
      myco_result <- as.data.frame(permanova)[1, ]
    }
    
    # Check if Myco results exist
    if (is.null(myco_result) || nrow(myco_result) == 0) {
      p <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No Myco effect found", size = 5) +
        theme_void() +
        labs(title = paste0("Cd", cd_level, "-Stage", stage_level))
      return(p)
    }
    
    # Extract R2 and p-value
    r2_value <- myco_result$R2 * 100
    p_value <- myco_result$`Pr(>F)`
    
    # Add significance markers
    significance <- case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    )
    
  }, error = function(e) {
    print(paste("Error in PERMANOVA for Cd", cd_level, "Stage", stage_level, ":", e$message))
    p <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "PERMANOVA error", size = 5) +
      theme_void() +
      labs(title = paste0("", cd_level, "-", stage_level))
    return(p)
  })
  
  # NMDS analysis
  set.seed(123)
  nmds <- ordinate(
    subset_samples,
    method = "NMDS",
    distance = "bray"
  )
  
  # Get stress value and keep three decimal places
  stress_value <- nmds$stress
  stress_text <- paste("Stress =", round(stress_value, 3))
  
  # Extract NMDS coordinates
  nmds_scores <- as.data.frame(scores(nmds)$sites)
  nmds_scores$Myco <- metadata$treat
  
  # Format R2 text, keep three significant digits
  r2_text <- ifelse(!is.na(r2_value), sprintf("R² = %.1f%%%s", r2_value, significance), "R² = NA")
  
  # Visualize NMDS
  p <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
    # Add confidence ellipses
    stat_ellipse(aes(fill = Myco, color = Myco), 
                 geom = "polygon", alpha = 0.2, level = 0.95) +
    # Add sample points
    geom_point(aes(color = Myco), size = 4, shape = 16) +
    # Set Myco colors
    scale_color_manual(values = c("AM" = "#D0908A", "NM" = "#77AAE4"), name = "Myco") +
    scale_fill_manual(values = c("AM" = "#D0908A", "NM" = "#77AAE4"), name = "Myco") +
    # Add PERMANOVA result
    annotate("text", x = Inf, y = Inf, 
             label = r2_text, 
             hjust = 1.1, vjust = 1.5, size = 12,
             fontface = "italic", color = "darkred") +
    # Add stress value
    annotate("text", x = Inf, y = -Inf, 
             label = stress_text, 
             hjust = 1.1, vjust = -0.5, size = 12,
             fontface = "italic") +
    # Set title and legend
    labs(title = paste0("", cd_level, "-", stage_level),
         x = "NMDS1", y = "NMDS2") +
    theme_bw(base_size = 14) +
    theme(
      panel.border = element_rect(size = 1.4, color = "black"),
      panel.grid = element_blank(),
      axis.title = element_text(size = 30),
      axis.text = element_text(size = 28),
      legend.text = element_text(size = 30),
      legend.title = element_text(size = 30, face = "plain"),
      legend.position = "right",
      legend.box = "vertical",
      legend.key.size = unit(2.5, "lines"),
      plot.title = element_text(hjust = 0.5, size = 30, face = "bold")
    )
  
  return(p)
}

# Get all Cd levels and stage levels
cd_levels <- unique(sample_data(pl_rarefied)$cd)
stage_levels <- unique(sample_data(pl_rarefied)$stage)
cd_levels <- sort(cd_levels)  # Ensure consistent order
stage_levels <- sort(stage_levels)  # Ensure consistent order


# Create plots for each Cd and stage combination
plot_list <- list()
plot_index <- 1

for (cd_level in cd_levels) {
  for (stage_level in stage_levels) {
    plot_list[[plot_index]] <- create_nmds_plot(cd_level, stage_level, pl_rarefied)
    plot_index <- plot_index + 1
  }
}

# Combine plots using ggarrange (4 rows x 5 columns)
combined_plot <- ggarrange(
  plotlist = plot_list,
  ncol = 5, nrow = 4,
  common.legend = TRUE,
  legend = "right"
)

# Display the combined plot
print(combined_plot)

# Save the combined plot
ggsave("out/nmds_by_cd_stage.png", combined_plot, 
       width = 13000, height = 7500, units = "px", dpi = 450, bg ="white")
ggsave("out/nmds_by_cd_stage.pdf", combined_plot, 
       width = 13000, height = 7500, units = "px", dpi = 450, bg ="white")


# Combine all beta diversity analysis result plots ----
all_plot <- ggarrange( final_plot,combined_plot,
                       labels = c("A", "B"),
                       font.label = list(size = 24, face = "bold"),
                       ncol = 1, nrow = 2,
                       heights = c(1, 1.4))

ggsave("out/beta_div.png", all_plot, 
       width = 12000, height = 11000, units = "px", dpi = 450, bg ="white")

ggsave("out/beta_div.pdf", all_plot, 
       width = 12000, height = 11000, units = "px", dpi = 450, bg ="white")

save.image("out/bata.RData")