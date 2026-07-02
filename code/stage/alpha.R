# Alpha diversity analysis code
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(forcats)

# Calculate Alpha diversity indices # time-consuming
alpha_div <- estimate_richness(pl_un_rarefied, measures = c("Chao1", "Shannon"))

# Add sample grouping information
sample_info <- as.data.frame(sample_data(pl_un_rarefied))
alpha_div <- cbind(alpha_div, sample_info)
print(levels(alpha_div$cd))

# Set factor levels
alpha_div$stage <- factor(alpha_div$stage, 
                          levels = c("stage 1", "stage 2", "stage 3", "stage 4", "stage 5"))
alpha_div$treat <- factor(alpha_div$treat, levels = c( "NM", "AM"))


# Define colors
am_color <- alpha("#84170B", 0.9)   
nm_color <- "#9EA5A9"  
am_fill <- alpha("#84170B", 0.6)  
nm_fill <- alpha("#9EA5A9", 0.5)  


# Create significance test function ----
calculate_significance <- function(data, y_var) {
  # Get the order of cd levels
  cd_levels <- levels(data$cd)
  
  # Create all possible combinations of cd and stage
  combinations <- expand.grid(
    cd = cd_levels,
    stage = levels(data$stage),
    stringsAsFactors = FALSE
  )
  
  # Perform t-test for each combination
  results <- apply(combinations, 1, function(row) {
    cd_val <- row["cd"]
    stage_val <- row["stage"]
    
    # Subset data
    subset_data <- data[data$cd == cd_val & data$stage == stage_val, ]
    
    # Ensure each group has at least 2 observations
    if (nrow(subset_data[subset_data$treat == "AM", ]) < 2 || 
        nrow(subset_data[subset_data$treat == "NM", ]) < 2) {
      return(data.frame(
        cd = cd_val,
        stage = stage_val,
        p.value = NA,
        p.label = NA,
        y.position = NA,
        x.position = NA
      ))
    }
    
    # Perform t-test
    t_test <- t.test(as.formula(paste(y_var, "~ treat")), data = subset_data)
    
    # Calculate y-axis position
    y_max <- max(subset_data[[y_var]], na.rm = TRUE)
    y_position <- y_max * 1.03
    
    # Determine p-value label
    p_label <- ifelse(t_test$p.value < 0.001, "***",
                      ifelse(t_test$p.value < 0.01, "**",
                             ifelse(t_test$p.value < 0.05, "*", "")))
    
    # Calculate x-axis position
    # Get the number of stage levels
    stage_levels <- levels(subset_data$stage)
    stage_index <- which(stage_levels == stage_val)
    
    return(data.frame(
      cd = cd_val,
      stage = stage_val,
      p.value = t_test$p.value,
      p.label = p_label,
      y.position = y_position,
      x.position = stage_index  
    ))
  })
  
  # Combine all results
  do.call(rbind, results)
}

# Calculate significance for Chao1 index
stat.test.chao1 <- calculate_significance(alpha_div, "Chao1")

# Calculate significance for Shannon index
stat.test.shannon <- calculate_significance(alpha_div, "Shannon")



# Plot Chao1 index boxplot ----
p_chao1 <- ggplot(alpha_div, aes(x = stage, y = Chao1, fill = treat, color = treat)) +  # add color mapping
  geom_boxplot(
    outlier.shape = NA, 
    position = position_dodge(0.9),
    size = 1.5,  
    alpha = 0.6  
  ) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9),
    size = 2,  
    alpha = 0.8  
  ) +
  facet_wrap(~fct_relevel(cd, "Cd 0", "Cd 2", "Cd 5", "Cd 15"), ncol = 4) +
  # Set fill and border colors separately
  scale_fill_manual(values = c( "NM" = nm_fill, "AM" = am_fill), 
                    labels = c("NM", "AM")) +
  scale_color_manual(values = c("NM" = nm_color, "AM" = am_color), 
                     labels = c("NM", "AM")) +
  # Add significance markers
  geom_text(data = stat.test.chao1[!is.na(stat.test.chao1$p.label) & stat.test.chao1$p.label != "", ],
            aes(x = x.position, y = y.position, label = p.label),
            inherit.aes = FALSE, size = 9, vjust = 0.5) +  
  labs(x = "", y = "Chao1", fill = "", color = "") +  
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 1.1, color = "black", fill = NA),
    axis.text.x = element_text(angle = 80, hjust = 0.5, vjust = 0.5, size = 17),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 22),
    legend.position = "top",
    legend.key.size = unit(1.2, "cm"),
    legend.text = element_text(size = 19),
    strip.background = element_rect(linewidth = 1.1, fill = "lightgrey"),
    strip.text = element_text(size = 19)
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 1.2))
  )

# Plot Shannon index boxplot ----
p_shannon <- ggplot(alpha_div, aes(x = stage, y = Shannon, fill = treat, color = treat)) +  # add color mapping
  geom_boxplot(
    outlier.shape = NA, 
    position = position_dodge(0.9),
    size = 1.5,  
    alpha = 0.6  
  ) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9),
    size = 2, 
    alpha = 0.8  
  ) +
  facet_wrap(~fct_relevel(cd, "Cd 0", "Cd 2", "Cd 5", "Cd 15"), ncol = 4) +
  # Set fill and border colors separately
  scale_fill_manual(values = c( "NM" = nm_fill, "AM" = am_fill), 
                    labels = c("NM", "AM")) +
  scale_color_manual(values = c("NM" = nm_color, "AM" = am_color), 
                     labels = c("NM", "AM")) +
  # Add significance markers
  geom_text(data = stat.test.shannon[!is.na(stat.test.shannon$p.label) & stat.test.shannon$p.label != "", ],
            aes(x = x.position, y = y.position, label = p.label),
            inherit.aes = FALSE, size = 9, vjust = 0.5) + 
  labs(x = "", y = "Shannon", fill = "", color = "") +  
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 1.1, color = "black", fill = NA),
    axis.text.x = element_text(angle = 80, hjust = 0.5, vjust = 0.5, size = 17),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 22),
    legend.position = "top",
    legend.key.size = unit(1.2, "cm"),
    legend.text = element_text(size = 19),
    strip.background = element_rect(size = 1.1, fill = "lightgrey"),
    strip.text = element_text(size = 19)
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 1.2)),
    color = "none"  # Hide the legend for border colors (to avoid duplication)
  ) +
  guides(fill = guide_legend(override.aes = list(size = 1.2)))  # Adjust legend symbol size

# Combine plots
alpha_combined <- ggarrange(p_chao1, p_shannon,
                            ncol = 1, nrow = 2,
                            labels = c("A", "B"),
                            common.legend = TRUE, legend = "top")

# Display plot
print(alpha_combined)

# Save plot
if(save){
  ggsave(filename = "out/alpha_diversity_combined.png", plot = alpha_combined, 
         width = 7000, height = 6500, dpi = 600, units = "px", bg = "white" )
}