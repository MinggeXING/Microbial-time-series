# PLOT Phylum level------------------------------------------------------------------------
# Plot relative abundance Phylum top 10----
# Use data.table for processing
library(data.table)
library(ggplot2)
library(dplyr)
library(ggpubr)

# ***Note this step can be time consuming ***
# Transform to relative abundances
ps.rel <- transform_sample_counts(pl_un_rarefied, function(x) x/sum(x))

# Get ASV table and metadata
ASV_data <- as(otu_table(ps.rel), "matrix")
metadata <- as(sample_data(ps.rel), "matrix")

# Convert to data.table
ASV_dt <- as.data.table(ASV_data, keep.rownames = "ASV")

# Convert to long format and merge with metadata
ASV_long <- melt(ASV_dt, id.vars = "ASV", variable.name = "sample", value.name = "Abundance")
merged_ASV <- merge(ASV_long, metadata, by = "sample")  # Merging is time-consuming and requires large memory; wait for memory release from previous task before proceeding.

rm(ASV_dt, ASV_long, ASV_data)

# Get taxonomy information
tax <- as.data.table(as.data.frame(tax_table(ps.rel)), keep.rownames = "ASV")

# Merge ASV and taxonomy
df.m <- merge(merged_ASV, tax, by = "ASV")

rm(ps.rel, tax, merged_ASV)


# Aggregate by sample and Phylum
phylum_abund <- df.m[, .(Abundance = sum(Abundance)), by = .(sample, Phylum)]

# Calculate total abundance for each Phylum
phylum_totals <- phylum_abund[, .(Total = sum(Abundance)), by = Phylum]

# Get the top 10 phyla
top10 <- phylum_totals[order(-Total)]$Phylum[1:10]

# Mark non-top-10 phyla as "Others"
phylum_abund[, Phylum.top10 := ifelse(Phylum %in% top10, Phylum, "Others")]

# Ensure Others is last
phylum_abund$Phylum.top10 <- factor(phylum_abund$Phylum.top10, 
                                    levels = c(top10, "Others"))

# Verify results
unique(phylum_abund$Phylum.top10)


# plot with ggplot2

colors <- c("Pseudomonadota"  = '#8dd3c7',
            "Bacillota"   = '#ffffb3',
            "Bacteroidota"   = '#fccde5',
            "Actinomycetota"   = '#80b1d3',
            "Chloroflexota"  = '#bebada',
            "Myxococcota" = '#fb8072',
            "Planctomycetota"   = '#fdb462',
            "Acidobacteriota"  = '#b3de69',
            "Verrucomicrobiota"  = '#d9d9d9',
            "Gemmatimonadota"  ='#bc80bd',
            "Others" = 'darkgrey'
)

Phylum.top10 <- factor(phylum_abund$Phylum.top10, levels = c("Pseudomonadota" ,
                                                             "Bacillota" ,
                                                             "Bacteroidota" ,
                                                             "Actinomycetota" ,
                                                             "Chloroflexota",
                                                             "Myxococcota" ,
                                                             "Planctomycetota" ,
                                                             "Acidobacteriota" ,
                                                             "Verrucomicrobiota" ,
                                                             "Gemmatimonadota" ,
                                                             "Others" 
))

# abundance of bacterial comm (every rep) ----
rel_ab <- ggplot(phylum_abund, aes(x = sample, y = Abundance, fill = Phylum.top10)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = colors, breaks = c("Pseudomonadota" ,
                                                "Bacillota" ,
                                                "Bacteroidota" ,
                                                "Actinomycetota" ,
                                                "Chloroflexota",
                                                "Myxococcota" ,
                                                "Planctomycetota" ,
                                                "Acidobacteriota" ,
                                                "Verrucomicrobiota" ,
                                                "Gemmatimonadota" ,
                                                "Others"
  )) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

rel_ab # time consuming



# Group
df.m$Group <- factor(df.m$Group, 
                     levels = c("S1AM0", "S1NM0", "S1AM2", "S1NM2", "S1AM5", "S1NM5", "S1AM15", "S1NM15",
                                "S2AM0", "S2NM0", "S2AM2", "S2NM2", "S2AM5", "S2NM5", "S2AM15", "S2NM15",
                                "S3AM0", "S3NM0", "S3AM2", "S3NM2", "S3AM5", "S3NM5", "S3AM15", "S3NM15",
                                "S4AM0", "S4NM0", "S4AM2", "S4NM2", "S4AM5", "S4NM5", "S4AM15", "S4NM15",
                                "S5AM0", "S5NM0", "S5AM2", "S5NM2", "S5AM5", "S5NM5", "S5AM15", "S5NM15"))
# cd
df.m$cd <- factor(df.m$cd, 
                  levels = c("Cd 0", "Cd 2", "Cd 5", "Cd 15"))
# treat
df.m$treat <- factor(df.m$treat, 
                     levels = c("AM", "NM"))
# stage
df.m$stage <- factor(df.m$stage, 
                     levels = c("stage 1", "stage 2", "stage 3", "stage 4", "stage 5"))

levels(df.m$Group)
levels(df.m$cd)
levels(df.m$treat)
levels(df.m$stage)



# create rel ab plot using average value of 9 replicates ----
# create a new column for correct grouping, 5 stage, 4 cd, and 2 Myco, so 5*4*2 = 40 groups for the plot in total
df.m$stage_cd_Myco <- paste(df.m$stage, df.m$cd, df.m$treat,sep = "_")
# check if 40 groups
length(unique(df.m$stage_cd_Myco)) == 40

# Get Phylum.top10 information from phylum_abund and merge into df.m
phylum_top10_info <- unique(phylum_abund[, .(Phylum, Phylum.top10)])
df.m <- merge(df.m, phylum_top10_info, by = "Phylum", all.x = TRUE)
# check if 10 phyla and "others"
length(unique(df.m$Phylum.top10)) == 11

# Extract AM and NM data
df.am <- df.m[treat == "AM"]
df.nm <- df.m[treat == "NM"]

# Calculate average relative abundance per group (average of 9 replicates)
df.am_avg <- df.am %>% 
  group_by(stage, cd, Phylum.top10) %>% 
  summarise(avg_rel_ab = sum(Abundance)/9, .groups = "drop")

df.nm_avg <- df.nm %>% 
  group_by(stage, cd, Phylum.top10) %>% 
  summarise(avg_rel_ab = sum(Abundance)/9, .groups = "drop")

# Set factor levels
df.am_avg$stage <- factor(df.am_avg$stage, levels = levels(df.m$stage))
df.am_avg$cd <- factor(df.am_avg$cd, levels = levels(df.m$cd))
df.nm_avg$stage <- factor(df.nm_avg$stage, levels = levels(df.m$stage))
df.nm_avg$cd <- factor(df.nm_avg$cd, levels = levels(df.m$cd))

# Calculate mean relative abundance for each treatment combination in AM and NM groups
setDT(df.am_avg)

am_relative_abundance <- df.am_avg[, .(
  avg_rel_ab = sum(avg_rel_ab)
), by = .(stage, cd, Phylum.top10)]

setDT(df.nm_avg)

nm_relative_abundance <- df.nm_avg[, .(
  avg_rel_ab = sum(avg_rel_ab)
), by = .(stage, cd, Phylum.top10)]

# Add treatment type column
am_relative_abundance[, treat := "AM"]
nm_relative_abundance[, treat := "NM"]

# Combine AM and NM data
combined_relative_abundance <- rbindlist(list(am_relative_abundance, nm_relative_abundance))

# Reorder columns
setcolorder(combined_relative_abundance, c("treat", "stage", "cd", "Phylum.top10", "avg_rel_ab"))

# Sort by treatment, stage, and Cd concentration
setorder(combined_relative_abundance, treat, stage, cd)

# Also create a wide format table for easier viewing
# Convert to wide format (each Phylum.top10 as a column)
wide_format_rel <- dcast(combined_relative_abundance, 
                         treat + stage + cd ~ Phylum.top10, 
                         value.var = "avg_rel_ab",
                         fill = 0)

# Export wide format data
fwrite(wide_format_rel, "out/relative_abundance_phylum.csv")

# Print result summary
print("Relative abundance data for AM and NM groups calculated and exported")
print(paste("Total rows processed:", nrow(combined_relative_abundance)))
print(paste("Including", length(unique(combined_relative_abundance$Phylum.top10)), "phylum-level classifications"))
print(paste("AM group:", nrow(am_relative_abundance), "rows"))
print(paste("NM group:", nrow(nm_relative_abundance), "rows"))

# Verify that the sum of relative abundance for each treatment combination is 1 (should be around 1, with small errors possible)
validation <- combined_relative_abundance[, .(total_rel_ab = sum(avg_rel_ab)), by = .(treat, stage, cd)]
print("Validation of relative abundance sums for each treatment combination:")
print(validation)


# Plot AM data ----
rel_ab_am <- ggplot(df.am_avg, aes(x = stage, y = avg_rel_ab, fill = Phylum.top10)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~cd, ncol = 4) +
  scale_fill_manual(values = colors) +
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

rel_ab_am

# Plot NM data ----
rel_ab_nm <- ggplot(df.nm_avg, aes(x = stage, y = avg_rel_ab, fill = Phylum.top10)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~cd, ncol = 4) +
  scale_fill_manual(values = colors) +
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

rel_ab_nm

# Combine plots using ggarrange
rel_ab_combined <- ggarrange(rel_ab_am, rel_ab_nm,
                             labels = c("A", "B"),
                             font.label = list(size = 24, face = "bold"),
                             ncol = 1, nrow = 2,
                             common.legend = TRUE, legend = "right")

# Save plot ----
if(save){
  ggsave(filename = "out/Phyla_top10_rel_AM_NM.png", plot = rel_ab_combined, 
         width = 12000, height = 10000, units = "px", dpi = 600, bg = "white")
}
# other class levels, like Class, Order, Phylum, Genus, Species are also plotted as previous
# relative abundance plot end----