# PLOT Genus level----
# Use data.table for processing
library(data.table)
library(ggplot2)
library(dplyr)
library(ggpubr)

# Calculate absolute abundance for Top 15 genera
# Use original phyloseq object
ps.abs <- pl_un_rarefied

# Get ASV table and metadata
ASV_data <- as(otu_table(ps.abs), "matrix")
metadata <- as(sample_data(ps.abs), "matrix")

# Convert to data.table
ASV_dt <- as.data.table(ASV_data, keep.rownames = "ASV")

# Convert to long format
ASV_long <- melt(ASV_dt, id.vars = "ASV", variable.name = "sample", value.name = "Abundance")

# Get taxonomy information
tax <- as.data.table(as.data.frame(tax_table(ps.abs)), keep.rownames = "ASV")

# Merge ASV and taxonomy
df.m2 <- merge(ASV_long, tax, by = "ASV") # time-consuming, after completion df.m2 will have 9 columns

# Clean memory
rm(ASV_dt, ASV_long, ASV_data, tax)

# Merge sample metadata
df.m2 <- merge(df.m2, metadata, by = "sample")

# Calculate total absolute abundance for each Genus
genus_totals <- df.m2[, .(Total = sum(Abundance)), by = Genus]

# Get the top 15 most abundant genera
abs_top15 <- genus_totals[order(-Total)]$Genus[1:15]

# Create Top 15 classification column
df.m2[, Genus.top15 := ifelse(Genus %in% abs_top15, Genus, "Others")] # time-consuming, after completion df.m2 will have 13 columns

# Ensure Others is last
df.m2$Genus.top15 <- factor(df.m2$Genus.top15, levels = c(abs_top15, "Others")) # time-consuming

# Verify results
table(df.m2$Genus.top15)

# Set factor levels for grouping variables
df.m2$Group <- factor(df.m2$Group, 
                      levels = c("S1AM0", "S1NM0", "S1AM2", "S1NM2", "S1AM5", "S1NM5", "S1AM15", "S1NM15",
                                 "S2AM0", "S2NM0", "S2AM2", "S2NM2", "S2AM5", "S2NM5", "S2AM15", "S2NM15",
                                 "S3AM0", "S3NM0", "S3AM2", "S3NM2", "S3AM5", "S3NM5", "S3AM15", "S3NM15",
                                 "S4AM0", "S4NM0", "S4AM2", "S4NM2", "S4AM5", "S4NM5", "S4AM15", "S4NM15",
                                 "S5AM0", "S5NM0", "S5AM2", "S5NM2", "S5AM5", "S5NM5", "S5AM15", "S5NM15")) # time-consuming

df.m2$cd <- factor(df.m2$cd, levels = c("Cd 0", "Cd 2", "Cd 5", "Cd 15"))
df.m2$treat <- factor(df.m2$treat, levels = c("AM", "NM"))
df.m2$stage <- factor(df.m2$stage, levels = c("stage 1", "stage 2", "stage 3", "stage 4", "stage 5"))

# Extract AM and NM data
df.am2 <- df.m2[treat == "AM"] # time-consuming
df.nm2 <- df.m2[treat == "NM"] # time-consuming

# Calculate average absolute abundance per group (average of 9 replicates)
df.am_avg2 <- df.am2[, .(avg_abs_ab = sum(Abundance)/9), by = .(stage, cd, Genus.top15)]
df.nm_avg2 <- df.nm2[, .(avg_abs_ab = sum(Abundance)/9), by = .(stage, cd, Genus.top15)]

# Set factor levels
df.am_avg2$stage <- factor(df.am_avg2$stage, levels = levels(df.m2$stage))
df.am_avg2$cd <- factor(df.am_avg2$cd, levels = levels(df.m2$cd))
df.nm_avg2$stage <- factor(df.nm_avg2$stage, levels = levels(df.m2$stage))
df.nm_avg2$cd <- factor(df.nm_avg2$cd, levels = levels(df.m2$cd))

# Save calculated workspace data----
save.image(file = "abs_genus.RData")

# Calculate mean absolute abundance for each treatment combination (stage × cd × Genus.top15) in AM and NM groups ----
am_absolute_abundance <- df.am_avg2[, .(
  avg_abs_ab = sum(avg_abs_ab)
), by = .(stage, cd, Genus.top15)]

nm_absolute_abundance <- df.nm_avg2[, .(
  avg_abs_ab = sum(avg_abs_ab)
), by = .(stage, cd, Genus.top15)]

# Add treatment type column
am_absolute_abundance[, treat := "AM"]
nm_absolute_abundance[, treat := "NM"]

# Combine AM and NM data
combined_absolute_abundance <- rbindlist(list(am_absolute_abundance, nm_absolute_abundance))

# Reorder columns
setcolorder(combined_absolute_abundance, c("treat", "stage", "cd", "Genus.top15", "avg_abs_ab"))

# Sort by treatment, stage, and Cd concentration
setorder(combined_absolute_abundance, treat, stage, cd)

# Also create a wide format table for easier viewing
# Convert to wide format (each Genus.top15 as a column)
wide_format <- dcast(combined_absolute_abundance, 
                     treat + stage + cd ~ Genus.top15, 
                     value.var = "avg_abs_ab",
                     fill = 0)

# Export wide format data
fwrite(wide_format, "out/absolute_abundance_genus.csv")

# Print result summary
print("Absolute abundance data for AM and NM groups calculated and exported")
print(paste("Total rows processed:", nrow(combined_absolute_abundance)))
print(paste("Including", length(unique(combined_absolute_abundance$Genus.top15)), "genus-level classifications"))
print(paste("AM group:", nrow(am_absolute_abundance), "rows"))
print(paste("NM group:", nrow(nm_absolute_abundance), "rows"))


# Visualization ----
# Calculate maximum y value for AM and NM data
max_am <- max(df.am_avg2[, .(total_ab = sum(avg_abs_ab)), by = .(stage, cd)]$total_ab)
max_nm <- max(df.nm_avg2[, .(total_ab = sum(avg_abs_ab)), by = .(stage, cd)]$total_ab)
max_y <- max(max_am, max_nm) * 1.05  # Add 5% margin

# Define colors - create color scheme for 15 genera
colors_genus <- c(
  "#8dd3c7", "#ffffb3", "#fccde5", "#80b1d3", "#bebada",
  "#fb8072", "#fdb462", "#b3de69", "#d9d9d9", "#bc80bd",
  "#ccebc5", "#ffed6f", "#e5d8bd", "#fddaec", "#f2f2f2",
  "darkgrey"  # Others
)

# Get actual levels of Genus.top15 present
actual_levels <- unique(c(levels(df.am_avg2$Genus.top15), levels(df.nm_avg2$Genus.top15)))
names(colors_genus) <- c(actual_levels[1:15], "Others")

# Plot absolute abundance for AM data
abs_ab_am <- ggplot(df.am_avg2, aes(x = stage, y = avg_abs_ab, fill = Genus.top15)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~cd, ncol = 4) +
  scale_fill_manual(values = colors_genus) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1)) +
  xlab("") +
  ylab("Absolute abundance") +
  # Set uniform y-axis range
  scale_y_continuous(labels = scales::comma, limits = c(0, max_y)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",  
    plot.margin = margin(r = unit(1, "cm")),  
    coord_cartesian(clip = "off"),
    line = element_line(linewidth = 1.5),
    axis.ticks = element_line(linewidth = 1.5),
    panel.border = element_rect(linewidth = 1.5),
    legend.key = element_rect(linewidth = 1.5),
    strip.text = element_text(size = 22, hjust = 0.5, vjust = 0.5, margin = margin(t = 8, b = 8)), 
    strip.background = element_rect(fill = "lightgrey", linewidth = 1.5, linetype = "solid"),
    axis.title.x = element_text(size = 24, colour = "black", margin = margin(t = 10)),  
    axis.title.y = element_text(size = 32, colour = "black", margin = margin(t = 10)),
    axis.text.x = element_text(size = 19, angle = 13),
    axis.text.y = element_text(size = 20, margin = margin(r = 15)),
    legend.title = element_text(size = 26),
    legend.text = element_text(size = 22),  # Reduce legend text size
    legend.key.size = unit(1.7, "cm"),      # Reduce legend symbol size
    legend.spacing.y = unit(0.4, "cm"),     # Reduce vertical legend spacing
    plot.title = element_text(size = 24)
  ) +
  guides(fill = guide_legend(ncol = 1)) +   # Single column legend
  labs(fill = "Top 15 Genera", title = "AM Treatment")

# Plot absolute abundance for NM data
abs_ab_nm <- ggplot(df.nm_avg2, aes(x = stage, y = avg_abs_ab, fill = Genus.top15)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~cd, ncol = 4) +
  scale_fill_manual(values = colors_genus) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1)) +
  xlab("") +
  ylab("Absolute abundance") +
  # Set uniform y-axis range
  scale_y_continuous(labels = scales::comma, limits = c(0, max_y)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",  
    plot.margin = margin(r = unit(1, "cm")),  
    coord_cartesian(clip = "off"),
    line = element_line(linewidth = 1.5),
    axis.ticks = element_line(linewidth = 1.5),
    panel.border = element_rect(linewidth = 1.5),
    legend.key = element_rect(linewidth = 1.5),
    strip.text = element_text(size = 22, hjust = 0.5, vjust = 0.5, margin = margin(t = 8, b = 8)), 
    strip.background = element_rect(fill = "lightgrey", linewidth = 1.5, linetype = "solid"),
    axis.title.x = element_text(size = 24, colour = "black", margin = margin(t = 10)),  
    axis.title.y = element_text(size = 32, colour = "black", margin = margin(t = 10)),
    axis.text.x = element_text(size = 19, angle = 13),
    axis.text.y = element_text(size = 20, margin = margin(r = 15)),
    legend.title = element_text(size = 26),
    legend.text = element_text(size = 22),  # Reduce legend text size
    legend.key.size = unit(1.7, "cm"),      # Reduce legend symbol size
    legend.spacing.y = unit(0.4, "cm"),     # Reduce vertical legend spacing
    plot.title = element_text(size = 24)
  ) +
  guides(fill = guide_legend(ncol = 1)) +   # Single column legend
  labs(fill = "Top 15 Genera", title = "NM Treatment")

# Combine plots using ggarrange
abs_ab_combined <- ggarrange(abs_ab_am, abs_ab_nm,
                             labels = c("A", "B"),
                             font.label = list(size = 24, face = "bold"),
                             ncol = 1, nrow = 2,
                             common.legend = TRUE, legend = "right")

# Save plot----
if(save){
  ggsave(filename = "out/Genera_top15_abs_AM_NM.png", plot = abs_ab_combined, 
         width = 14000, height = 12000, units = "px", dpi = 600, bg = "white")
}