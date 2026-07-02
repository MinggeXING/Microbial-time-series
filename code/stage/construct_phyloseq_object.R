# Construct phyloseq object and plot relative abundance bar chart ----

# Analyzing bacterial data using ASV table by unoise3----
# Note: To calculate alpha diversity, a non-normalized asv table should be better. 
# For PCoA, a normalized asv table is preferred. 
# Read: https://www.bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html

# This analysis requires the packages in my_packages. Install them and proceed. 
# Install and load packages, set plot theme, and source own functions ----
# install 'Rtools', search website, download, and install
# install 'phyloseq' and 'DESeq2'

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
BiocManager::install("DESeq2")

# install 'pairwiseAdonis'
install.packages("htmlwidgets") # needed for 'devtools'
install.packages("devtools") # use the package panel of Rstudio to install
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis") # require 'devtools'


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("apeglm")

remotes::install_github("mlverse/chattr")
devtools::install_github("zdk123/SpiecEasi")

my_packages <- c('readxl', 'phyloseq', 'ape', 'Biostrings', 'tidyverse', 'dplyr', 'pairwiseAdonis', 'reshape2', 'ggplot2', 'ggpubr', 'DESeq2', 'GGally', 'igraph', 'vegan', 'ecodist', 'agricolae', 'stats', 'writexl', 'dendextend', 'apeglm', 'pheatmap', 'GOplot','tibble', 'boot', 'chattr', 'SpiecEasi')

install.packages(my_packages) # 'phyloseq', 'DESeq2', and 'pairwiseAdonis' may not be installed using this code. But it is fine as they have been installed above. my_packages is also used for library attaching as follows.

#################################################################################
# rarefy----
# Since the rarefaction value (528969445) and sample size (360) are too large, local computation is not feasible. The asv rarefied data (asv.abs.stage.rarefied.csv) was obtained through rarefaction on the server using Qiime 2 environment; see code in run_rarefy.txt

# START -------------
# If all required packages are installed, start from here

# Clear the environment
rm(list = ls())

# Clear the console
cat("\014")  # This is equivalent to pressing Ctrl+L

# Clear all plots
if (!is.null(dev.list())) dev.off()

# save or not
save <- FALSE

########################################################################################
# asv_abs_stage.csv is filtered from the original data (as corr_networks.R calculated), the original data from major seq company's cloud analysis platform :515F_806R-540samples\workflow_results\3.ASV_Analysis\ASV,which is rep-asv after denoise by DADA2.
# As Phyloseq must be an integer to calculate alpha diversity, the value in the asv_table given by the sequencing company is a decimal place, so it is rounded to the nearest whole number
# un rarefy data rounded----
if (!require("readxl")) install.packages("readxl")
if (!require("dplyr")) install.packages("dplyr")
library(readxl)
library(dplyr)
df <- read.csv("data/stage_data/asv.abs.stage.rarefied.csv", row.names = 1)
df_rounded <- df %>% mutate(across(where(is.numeric), ~ round(., digits = 0)))
head(df_rounded)
#write.csv(df_rounded, "data/stage networks/asv_abs_stage_rounded.csv",row.names = TRUE)
# Open the exported file and change the first column name to "ID"
rm(df,df_rounded)

########################################################################################
# Load packages needed
my_packages <- c( 'phyloseq', 'ape', 'Biostrings', 'tidyverse', 
                  'dplyr', 'pairwiseAdonis', 'reshape2', 'ggplot2', 'ggpubr', 
                  'DESeq2', 'GGally', 'igraph', 'vegan', 'ecodist', 'agricolae', 
                  'dendextend', 'apeglm', 'pheatmap', 
                  'GOplot','tibble', 'boot', 'chattr', 'SpiecEasi')

lapply(my_packages, library, character.only = TRUE) # load my_packages

install.packages()
theme_set(theme_bw()) # Set plot theme of ggplot2. It can be changed to other themes.
source('code/fun.R') # load own function(s)

# check versions of R and other essential packages
R.Version()
packageVersion('phyloseq') # version 1.41.1
packageVersion('DESeq2') # version 1.42.1
packageVersion('vegan') # version 2.7.1

# if use chattr for integrated ChatGPT, run the following
chattr_app()


# ~ Load raw data ----
# prepare 5 files according to required format before proceeding:
# 1. asv_table
# 2. tax
# 3. asv_tree
# 4. metadata
# 5. seq.fasta


# after obtaining the above 5 files in hands, proceed as follows
# ~~ 1. load ASV table ----
asv_table <- read.csv("data/stage_data/asv_abs_stage_rounded.csv") # asv table stage - un-rarefied. 
# The un-rarefied table will be used to calculate alpha-diversity, while the rarefied one will be used to analyze beta-diversity

asv_table_rarefy <- read.csv("data/stage_data/asv.abs.stage.rarefied.csv")  # Server Qiime 2 rarefied, asv table - rarefied # Locally open the file, delete the first line of irrelevant content, rename the first column header to "ID", then read.

# check if the otu table has been normalized
see_if_norm <- asv_table %>% 
  summarise(across(S1AM0_1:S5NM15_9, ~ sum(.x, na.rm = TRUE)))
see_if_norm # not normalized 
min(see_if_norm); max(see_if_norm) # check the smallest read count
rm(see_if_norm) # remove the temp file 'see_if_norm'


see_if_norm <- asv_table_rarefy %>% 
  summarise(across(S1AM0_1:S5NM15_9, ~ sum(.x, na.rm = TRUE)))
see_if_norm # not normalized 
min(see_if_norm); max(see_if_norm) # check the smallest read count
rm(see_if_norm) # remove the temp file 'see_if_norm'

# ~~ 2. load tax info ----
tax <- read.csv('data/stage_data/tax_matched.csv') 
#asv_table should be match with tax,so it is another .csv

# ~~ 3. load tree file ----
tree <- read.tree('data/stage_data/phylo.tre') # ape package is needed

# ~~ 4. load metadata ----
metadata <- read.csv('data/stage_data/stage_metadata.csv')

# ~~ 5. load seqs data ----
# 'ASV_rep-seqs.qza' # obtained using qiime2
# The format is like >ASV1
# ATTGGACAATGGGCGCAAGCCTGATCCAGCCATGCCG....
# There are 148651 ASV in total
seq <- readDNAStringSet('data/stage_data/ASV_reps.fasta', format="fasta",
                        nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)


# ~~ Prepare ASV table matrix ----
# should be matrix, rownames are ASV IDs, colnames are sample IDs
# it is required if need to use phyloseq package.
# convert to matrix
# get the column ASV ID, then assign as row names
asv_table2 <- column_to_rownames(asv_table, var = "ID")
asv_table_rarefy2 <- column_to_rownames(asv_table_rarefy, var = "ID")
# convert data frame to matrix
asv_mat <- data.matrix(asv_table2)
asv_rarefy_mat <- data.matrix(asv_table_rarefy2)
rm(asv_table, asv_table2, asv_table_rarefy,asv_table_rarefy2) # remove temp files. 


# ~~ Prepare tax file matrix ----
# dealing with tax info
# get the column ASV ID, then assign as row names
tax2 <- column_to_rownames(tax, var = "ID")
# convert data frame to matrix
# use as.matrix rather than data.matrix
# otherwise, text will be deleted
tax_mat <- as.matrix(tax2) 
rm(tax2,tax) # remove temp file 'tax2'


# ~~ Prepare metadata dataframe ----
metadata <- as.data.frame(metadata)

metadata_df <- column_to_rownames(metadata, var = "ID")

# convert things in metadata_df to factors
metadata_df$Group <- as.factor(metadata_df$Group)
metadata_df$cd <- as.factor(metadata_df$cd)
metadata_df$treat <- as.factor(metadata_df$treat)
metadata_df$stage <- as.factor(metadata_df$stage)

# set factor class. The later plotting shall follow the set class
# Assuming metadata_df is your dataframe
metadata_df$Group <- factor(metadata_df$Group, levels = c("S1AM0", "S1NM0", "S1AM2", "S1NM2", "S1AM5", "S1NM5", "S1AM15", "S1NM15",
                                                          "S2AM0", "S2NM0", "S2AM2", "S2NM2", "S2AM5", "S2NM5", "S2AM15", "S2NM15",
                                                          "S3AM0", "S3NM0", "S3AM2", "S3NM2", "S3AM5", "S3NM5", "S3AM15", "S3NM15",
                                                          "S4AM0", "S4NM0", "S4AM2", "S4NM2", "S4AM5", "S4NM5", "S4AM15", "S4NM15",
                                                          "S5AM0", "S5NM0", "S5AM2", "S5NM2", "S5AM5", "S5NM5", "S5AM15", "S5NM15"))
metadata_df$cd <- factor(metadata_df$cd, levels = c( "Cd 0", "Cd 2", "Cd 5", "Cd 15"))
metadata_df$treat <- factor(metadata_df$treat, levels = c("AM", "NM"))
metadata_df$stage <- factor(metadata_df$stage, levels = c("stage 1", "stage 2", "stage 3", "stage 4", "stage 5"))

# Check the levels of the different factors
levels(metadata_df$Group)
levels(metadata_df$cd)
levels(metadata_df$treat)
levels(metadata_df$stage)

rm(metadata) # remove temp file


# ~~ Prepare tree and seqs data ----
# no need to prepare(?)


# CHECK POINT 1 ---------------------------------------------------------------
# 1. asv_mat - has asv id and sample id -- OK ******** it should be a matrix
head(asv_mat)
class(asv_mat)

head(asv_rarefy_mat)
class(asv_rarefy_mat)

# 2. tax_mat - has asv id and taxon info   ******** it should be a matrix
head(tax_mat)
class(tax_mat)

# 3. metadata - has sample id and treatments -- OK ******** it should be a data.frame
head(metadata_df)
class(metadata_df)


# 4. tree - has asv id -- OK(?) how to check(?)
taxa_names(tree)
class(tree)

# 5. seq - has DNA sequences and otu id -- OK(?) how to check(?)
taxa_names(seq)
class(seq)


### the above 5 items can be combined as a phylo object for further analysis


# IMPORTANT: telling how to combine the 5 items #### 
ASV <- otu_table(asv_mat, taxa_are_rows = TRUE) # IMPORTANT step ***
rm(asv_mat)

ASV2 <- otu_table(asv_rarefy_mat, taxa_are_rows = TRUE) # IMPORTANT step ***
rm(asv_rarefy_mat)

TAX <- tax_table(tax_mat) # IMPORTANT step ***
rm(tax_mat)

metadata <- sample_data(metadata_df) # IMPORTANT step ***
rm(metadata_df) # remove temp file

# TREE <- phy_tree(tree) # actually no need, seems 'tree' =? 'TREE'

# using the following to check type

class(ASV)
class(TAX)
class(ASV2)
class(metadata)
class(tree)
class(seq)

# Now in the 'Global Environment' panel, there should be 5 objects:ASV, ASV2, TAX, metadata, tree, seq

# Combine to construct phyloseq object ----------------------------------------
pl_un_rarefied <- phyloseq(ASV, TAX, metadata, tree, seq)
pl_rarefied <- phyloseq(ASV2, TAX, metadata, tree, seq)
# if return error saying taxa/OTU names are not match. Use the followings to check
# a <- taxa_names(ASV)
# b <- taxa_names(TAX)
# c <- taxa_names(seq)
# d <- taxa_names(tree)
# setdiff(a, b); setdiff(a, c); setdiff(a, d) etc.

# seems probably the taxa_names of tree are different from other phyloseq objects
# ASV vs asv; but in tree seems OK. For example, check:
# setdiff(taxa_names(ASV),taxa_names(tree))

# see this post: https://github.com/joey711/phyloseq/issues/1044


# CHECK POINT 2 ---------------------------------------------------------------
# inspect different names and variables

pl_un_rarefied
sample_names(pl_un_rarefied)
rank_names(pl_un_rarefied)
sample_variables(pl_un_rarefied)
taxa_names(pl_un_rarefied)
refseq(pl_un_rarefied)

pl_rarefied
sample_names(pl_rarefied)
rank_names(pl_rarefied)
sample_variables(pl_rarefied)
taxa_names(pl_rarefied)
refseq(pl_rarefied)

# check again the min and max read counts
# If the data is not rarefied, min != max.
print(paste('min =', min(sample_sums(pl_un_rarefied)))); print(paste('max =', max(sample_sums(pl_un_rarefied)))) # min != max -- un-rarefied
print(paste('min =', min(sample_sums(pl_rarefied)))); print(paste('max =', max(sample_sums(pl_rarefied)))) # min != max -- rarefied

# remove unnecessary files
rm(metadata, tree, ASV, ASV2, seq)

# now there should be only two phyloseq objects left for further analyses. 
# check
pl_un_rarefied
pl_rarefied
# OK

# Save the created phyloseq object
save.image("phyloseq project.RData")