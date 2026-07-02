# PLOT Genus level----
# Use data.table for processing
library(data.table)
library(ggplot2)
library(dplyr)
library(ggpubr)

# Calculate relative abundance for Top 15 genera----
# Use phyloseq object after relative abundance transformation
ps.rel <- transform_sample_counts(pl_un_rarefied, function(x) x/sum(x))

# Get ASV table and metadata
ASV_data_rel <- as(otu_table(ps.rel), "matrix")
metadata_rel <- as(sample_data(ps.rel), "matrix")

# Convert to data.table
ASV_dt_rel <- as.data.table(ASV_data_rel, keep.rownames = "ASV")

# Convert to long format
ASV_long_rel <- melt(ASV_dt_rel, id.vars = "ASV", variable.name = "sample", value.name = "Abundance")

# Get taxonomy information
tax_rel <- as.data.table(as.data.frame(tax_table(ps.rel)), keep.rownames = "ASV")

# Merge ASV and taxonomy
df.m_rel <- merge(ASV_long_rel, tax_rel, by = "ASV") # time-consuming; wait for memory release after all time-consuming steps before running again

# Clean memory
rm(ASV_dt_rel, ASV_long_rel, ASV_data_rel, tax_rel)

# Merge sample metadata
df.m_rel <- merge(df.m_rel, metadata_rel, by = "sample") # time-consuming

# Calculate total relative abundance for each Genus
genus_totals_rel <- df.m_rel[, .(Total = sum(Abundance)), by = Genus]

# Get the top 15 most abundant genera
rel_top15 <- genus_totals_rel[order(-Total)]$Genus[1:15]

# Create Top 15 classification column
df.m_rel[, Genus.top15 := ifelse(Genus %in% rel_top15, Genus, "Others")] # time-consuming

# Ensure Others is last
df.m_rel$Genus.top15 <- factor(df.m_rel$Genus.top15, levels = c(rel_top15, "Others")) # time-consuming

# Verify results
table(df.m_rel$Genus.top15)

# Set factor levels for grouping variables # time-consuming
# Group
df.m_rel$Group <- factor(df.m_rel$Group, 
                         levels = c("S1AM0", "S1NM0", "S1AM2", "S1NM2", "S1AM5", "S1NM5", "S1AM15", "S1NM15",
                                    "S2AM0", "S2NM0", "S2AM2", "S2NM2", "S2AM5", "S2NM5", "S2AM15", "S2NM15",
                                    "S3AM0", "S3NM0", "S3AM2", "S3NM2", "S3AM5", "S3NM5", "S3AM15", "S3NM15",
                                    "S4AM0", "S4NM0", "S4AM2", "S4NM2", "S4AM5", "S4NM5", "S4AM15", "S4NM15",
                                    "S5AM0", "S5NM0", "S5AM2", "S5NM2", "S5AM5", "S5NM5", "S5AM15", "S5NM15"))
# cd
df.m_rel$cd <- factor(df.m_rel$cd, levels = c("Cd 0", "Cd 2", "Cd 5", "Cd 15"))
# treat
df.m_rel$treat <- factor(df.m_rel$treat, levels = c("AM", "NM"))
# stage
df.m_rel$stage <- factor(df.m_rel$stage, levels = c("stage 1", "stage 2", "stage 3", "stage 4", "stage 5"))

levels(df.m_rel$Group)
levels(df.m_rel$cd)
levels(df.m_rel$treat)
levels(df.m_rel$stage)

# Extract AM and NM data
df.am_rel <- df.m_rel[treat == "AM"] # time-consuming
df.nm_rel <- df.m_rel[treat == "NM"] # time-consuming

# Calculate average relative abundance per group (average of 9 replicates)
df.am_avg_rel <- df.am_rel[, .(avg_rel_ab = sum(Abundance)/9), by = .(stage, cd, Genus.top15)]
df.nm_avg_rel <- df.nm_rel[, .(avg_rel_ab = sum(Abundance)/9), by = .(stage, cd, Genus.top15)]

# Set factor levels
df.am_avg_rel$stage <- factor(df.am_avg_rel$stage, levels = levels(df.m_rel$stage))
df.am_avg_rel$cd <- factor(df.am_avg_rel$cd, levels = levels(df.m_rel$cd))
df.nm_avg_rel$stage <- factor(df.nm_avg_rel$stage, levels = levels(df.m_rel$stage))
df.nm_avg_rel$cd <- factor(df.nm_avg_rel$cd, levels = levels(df.m_rel$cd))

# Calculate mean relative abundance for each treatment combination in AM and NM groups
setDT(df.am_avg_rel)

am_relative_abundance <- df.am_avg_rel[, .(
  avg_rel_ab = sum(avg_rel_ab)
), by = .(stage, cd, Genus.top15)]

setDT(df.nm_avg_rel)

nm_relative_abundance <- df.nm_avg_rel[, .(
  avg_rel_ab = sum(avg_rel_ab)
), by = .(stage, cd, Genus.top15)]

# Add treatment type column
am_relative_abundance[, treat := "AM"]
nm_relative_abundance[, treat := "NM"]

# Combine AM and NM data
combined_relative_abundance <- rbindlist(list(am_relative_abundance, nm_relative_abundance))

# Reorder columns
setcolorder(combined_relative_abundance, c("treat", "stage", "cd", "Genus.top15", "avg_rel_ab"))

# Sort by treatment, stage, and Cd concentration
setorder(combined_relative_abundance, treat, stage, cd)

# Also create a wide format table for easier viewing
# Convert to wide format (each Genus.top15 as a column)
wide_format_rel <- dcast(combined_relative_abundance, 
                         treat + stage + cd ~ Genus.top15, 
                         value.var = "avg_rel_ab",
                         fill = 0)

# Export wide format data
fwrite(wide_format_rel, "out/relative_abundance_genus.csv")

# Print result summary
print("Relative abundance data for AM and NM groups calculated and exported")
print(paste("Total rows processed:", nrow(combined_relative_abundance)))
print(paste("Including", length(unique(combined_relative_abundance$Genus.top15)), "genus-level classifications"))
print(paste("AM group:", nrow(am_relative_abundance), "rows"))
print(paste("NM group:", nrow(nm_relative_abundance), "rows"))

# Verify that the sum of relative abundance for each treatment combination is 1 (should be around 1, with small errors possible)
validation <- combined_relative_abundance[, .(total_rel_ab = sum(avg_rel_ab)), by = .(treat, stage, cd)]
print("Validation of relative abundance sums for each treatment combination:")
print(validation)



# Define colors - create color scheme for 15 genera
colors_genus <- c(
  "#8dd3c7", "#ffffb3", "#fccde5", "#80b1d3", "#bebada",
  "#fb8072", "#fdb462", "#b3de69", "#d9d9d9", "#bc80bd",
  "#ccebc5", "#ffed6f", "#e5d8bd", "#fddaec", "#f2f2f2",
  "darkgrey"  # Others
)

# Get actual levels of Genus.top15 present
actual_levels <- unique(c(levels(df.am_avg_rel$Genus.top15), levels(df.nm_avg_rel$Genus.top15)))
names(colors_genus) <- c(actual_levels[1:15], "Others")

# Plot relative abundance for AM data ----
rel_ab_am <- ggplot(df.am_avg_rel, aes(x = stage, y = avg_rel_ab, fill = Genus.top15)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~cd, ncol = 4) +
  scale_fill_manual(values = colors_genus) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1)) +
  xlab("") +
  ylab("Relative abundance") +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
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
    legend.text = element_text(size = 22), 
    legend.key.size = unit(1.7, "cm"),     
    legend.spacing.y = unit(0.4, "cm"),     
    plot.title = element_text(size = 24)
  ) +
  guides(fill = guide_legend(ncol = 1)) +  
  labs(fill = "Top 15 Genera", title = "AM Treatment")

# Plot relative abundance for NM data ----
rel_ab_nm <- ggplot(df.nm_avg_rel, aes(x = stage, y = avg_rel_ab, fill = Genus.top15)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~cd, ncol = 4) +
  scale_fill_manual(values = colors_genus) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1)) +
  xlab("") +
  ylab("Relative abundance") +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
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
    legend.text = element_text(size = 22),  
    legend.key.size = unit(1.7, "cm"),      
    legend.spacing.y = unit(0.4, "cm"),    
    plot.title = element_text(size = 24)
  ) +
  guides(fill = guide_legend(ncol = 1)) +   
  labs(fill = "Top 15 Genera", title = "NM Treatment")

# Combine plots using ggarrange
rel_ab_combined <- ggarrange(rel_ab_am, rel_ab_nm,
                             labels = c("A", "B"),
                             font.label = list(size = 24, face = "bold"),
                             ncol = 1, nrow = 2,
                             common.legend = TRUE, legend = "right")

# Save plot----
if(save){
  ggsave(filename = "out/Genera_top15_rel_AM_NM.png", plot = rel_ab_combined, 
         width = 14000, height = 12000, units = "px", dpi = 600, bg = "white")
}