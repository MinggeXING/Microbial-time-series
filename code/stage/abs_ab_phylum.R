#plot absolute abundance result----
# Plot the bar plot for all phyla
# ***Note it can take time to run depending on the computer spec.***
# plot_bar(pl_un_rarefied, "sample_name", fill="Phylum") + geom_bar(stat="identity", color = NA, width = 1)
# the bar plot was exported and edited in AI/CorelDraw to add group names and labels for axes. Relative zotu count was converted to % in AI/CorelDraw. 

# ~ calculate actual relative abundance at the Phylum level ----

# Below code snippet demonstrate how to achieve this.
# The arrange(), rename() and select() are from the dplyr package and spread() is from the tidyr package, both packages are part of the tidyverse.
# The select() and spread() are used to convet the output from long to wide format.
# ref: https://github.com/joey711/phyloseq/issues/1521

# Plot absolute abundance Phylum top 10----
# Use data.table for processing
library(data.table)
library(ggplot2)
library(dplyr)
library(ggpubr)

# ***Note this step can be time consuming ***
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
df.m2 <- merge(ASV_long, tax, by = "ASV") # Merging is time-consuming and requires large memory; wait for memory release from previous task before proceeding.

# Clean memory
rm(ASV_dt, ASV_long, ASV_data, tax)

# Merge sample metadata
df.m2 <- merge(df.m2, metadata, by = "sample")

# Calculate total absolute abundance for each Phylum
phylum_totals <- df.m2[, .(Total = sum(Abundance)), by = Phylum]

# Get the top 10 most abundant phyla
abs_top10 <- phylum_totals[order(-Total)]$Phylum[1:10]

# Create Top 10 classification column
df.m2[, Phylum.top10 := ifelse(Phylum %in% abs_top10, Phylum, "Others")]

# Ensure Others is last
df.m2$Phylum.top10 <- factor(df.m2$Phylum.top10, levels = c(abs_top10, "Others"))

# Verify results
table(df.m2$Phylum.top10)

# Set factor levels for grouping variables
df.m2$Group <- factor(df.m2$Group, 
                      levels = c("S1AM0", "S1NM0", "S1AM2", "S1NM2", "S1AM5", "S1NM5", "S1AM15", "S1NM15",
                                 "S2AM0", "S2NM0", "S2AM2", "S2NM2", "S2AM5", "S2NM5", "S2AM15", "S2NM15",
                                 "S3AM0", "S3NM0", "S3AM2", "S3NM2", "S3AM5", "S3NM5", "S3AM15", "S3NM15",
                                 "S4AM0", "S4NM0", "S4AM2", "S4NM2", "S4AM5", "S4NM5", "S4AM15", "S4NM15",
                                 "S5AM0", "S5NM0", "S5AM2", "S5NM2", "S5AM5", "S5NM5", "S5AM15", "S5NM15"))

df.m2$cd <- factor(df.m2$cd, levels = c("Cd 0", "Cd 2", "Cd 5", "Cd 15"))
df.m2$treat <- factor(df.m2$treat, levels = c("AM", "NM"))
df.m2$stage <- factor(df.m2$stage, levels = c("stage 1", "stage 2", "stage 3", "stage 4", "stage 5"))

# Extract AM and NM data
df.am2 <- df.m2[treat == "AM"]
df.nm2 <- df.m2[treat == "NM"]

# Calculate average absolute abundance per group (average of 9 replicates)
df.am_avg2 <- df.am2[, .(avg_abs_ab = sum(Abundance)/9), by = .(stage, cd, Phylum.top10)]
df.nm_avg2 <- df.nm2[, .(avg_abs_ab = sum(Abundance)/9), by = .(stage, cd, Phylum.top10)]

# Set factor levels
df.am_avg2$stage <- factor(df.am_avg2$stage, levels = levels(df.m2$stage))
df.am_avg2$cd <- factor(df.am_avg2$cd, levels = levels(df.m2$cd))
df.nm_avg2$stage <- factor(df.nm_avg2$stage, levels = levels(df.m2$stage))
df.nm_avg2$cd <- factor(df.nm_avg2$cd, levels = levels(df.m2$cd))

# Calculate mean absolute abundance for each treatment combination (stage × cd × Phylum.top10) in AM and NM groups ----
am_absolute_abundance <- df.am_avg2[, .(
  avg_abs_ab = sum(avg_abs_ab)
), by = .(stage, cd, Phylum.top10)]

nm_absolute_abundance <- df.nm_avg2[, .(
  avg_abs_ab = sum(avg_abs_ab)
), by = .(stage, cd, Phylum.top10)]

# Add treatment type column
am_absolute_abundance[, treat := "AM"]
nm_absolute_abundance[, treat := "NM"]

# Combine AM and NM data
combined_absolute_abundance <- rbindlist(list(am_absolute_abundance, nm_absolute_abundance))

# Reorder columns
setcolorder(combined_absolute_abundance, c("treat", "stage", "cd", "Phylum.top10", "avg_abs_ab"))

# Sort by treatment, stage, and Cd concentration
setorder(combined_absolute_abundance, treat, stage, cd)

# Also create a wide format table for easier viewing
# Convert to wide format (each Phylum.top10 as a column)
wide_format <- dcast(combined_absolute_abundance, 
                     treat + stage + cd ~ Phylum.top10, 
                     value.var = "avg_abs_ab",
                     fill = 0)

# Export wide format data
fwrite(wide_format, "out/absolute_abundance_phylum.csv")

# Print result summary
print("Absolute abundance data for AM and NM groups calculated and exported")
print(paste("Total rows processed:", nrow(combined_absolute_abundance)))
print(paste("Including", length(unique(combined_absolute_abundance$Phylum.top10)), "phylum-level classifications"))
print(paste("AM group:", nrow(am_absolute_abundance), "rows"))
print(paste("NM group:", nrow(nm_absolute_abundance), "rows"))


# Visualization ----
# Define colors (adjust according to your absolute abundance top 10 phyla)
colors_abs <- c(
  "Pseudomonadota" = '#8dd3c7',
  "Bacillota" = '#ffffb3',
  "Bacteroidota" = '#fccde5',
  "Actinomycetota" = '#80b1d3',
  "Chloroflexota" = '#bebada',
  "Myxococcota" = '#fb8072',
  "Planctomycetota" = '#fdb462',
  "Acidobacteriota" = '#b3de69',
  "Verrucomicrobiota" = '#d9d9d9',
  "Gemmatimonadota" = '#bc80bd',
  "Others" = 'darkgrey'
)

# Ensure colors match Phylum.top10 levels
# Get actual levels of Phylum.top10 present
actual_levels <- unique(c(levels(df.am_avg2$Phylum.top10), levels(df.nm_avg2$Phylum.top10)))
colors_abs <- colors_abs[names(colors_abs) %in% actual_levels]

# Calculate maximum y value for AM and NM data
max_am <- max(df.am_avg2[, .(total_ab = sum(avg_abs_ab)), by = .(stage, cd)]$total_ab)
max_nm <- max(df.nm_avg2[, .(total_ab = sum(avg_abs_ab)), by = .(stage, cd)]$total_ab)
max_y <- max(max_am, max_nm) * 1.05  # Add 5% margin

# Plot absolute abundance for AM data ----
abs_ab_am <- ggplot(df.am_avg2, aes(x = stage, y = avg_abs_ab, fill = Phylum.top10)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~cd, ncol = 4) +
  scale_fill_manual(values = colors_abs) +
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
    axis.title.y = element_text(size = 28, colour = "black", margin = margin(t = 10)),
    axis.text.x = element_text(size = 16, angle = 15),
    axis.text.y = element_text(size = 18, margin = margin(r = 15)),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 22),
    legend.key.size = unit(2, "cm"),
    legend.spacing.y = unit(0.4, "cm"),
    plot.title = element_text(size = 22)
  ) +
  labs(fill = "Top 10 Phyla", title = "AM Treatment")

# Plot absolute abundance for NM data ----
abs_ab_nm <- ggplot(df.nm_avg2, aes(x = stage, y = avg_abs_ab, fill = Phylum.top10)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~cd, ncol = 4) +
  scale_fill_manual(values = colors_abs) +
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
    axis.title.y = element_text(size = 28, colour = "black", margin = margin(t = 10)),
    axis.text.x = element_text(size = 16, angle = 15),
    axis.text.y = element_text(size = 18, margin = margin(r = 15)),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 22),
    legend.key.size = unit(2, "cm"),
    legend.spacing.y = unit(0.4, "cm"),
    plot.title = element_text(size = 22)
  ) +
  labs(fill = "Top 10 Phyla", title = "NM Treatment")

# Combine plots using ggarrange
abs_ab_combined <- ggarrange(abs_ab_am, abs_ab_nm,
                             labels = c("A", "B"),
                             font.label = list(size = 24, face = "bold"),
                             ncol = 1, nrow = 2,
                             common.legend = TRUE, legend = "right")

# Save plot ----
if(save){
  ggsave(filename = "out/Phyla_top10_abs_AM_NM.png", plot = abs_ab_combined, 
         width = 12000, height = 10000, units = "px", dpi = 600, bg = "white")
}