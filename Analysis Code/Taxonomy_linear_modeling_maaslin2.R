# Analysis of rootstock trial metagenomes
# Code by: Joel F. Swift

#Packages w/ version numbers.
library('tidyverse'); packageVersion('tidyverse')
library('ggpubr'); packageVersion('ggpubr')
library('lmerTest'); packageVersion('lmerTest')
library('phyloseq'); packageVersion('phyloseq')
library('emmeans'); packageVersion('emmeans')
library('lubridate'); packageVersion('lubridate')

#### Color palettes ####
theme_set(theme_bw())
rootstock_palette<- c("101-14 MGT" = "#46436b", "3309 Couderc" = "#7570B3","Schwarzmann" = "#aca9d1",
                      "110 Richter" = "#0a3d2e", "1103 Paulsen" = "#1b9e77",
                      "140Ru" = "#48e0b3",  "775 Paulsen" = "#d3f8ed",
                      "420A" = "#453301", "Kober 5BB" = "#E6AB02", "Teleki 5C" = "#f0cd67")

#### Functions ####
# https://gist.github.com/lwaldron/512d1925a8102e921f05c5b25de7ec94/da56b85f0de9761d4e1a8e28b3319b386132114d
metaphlanToPhyloseq <- function(
    tax,
    metadat=NULL,
    simplenames=TRUE,
    roundtointeger=FALSE,
    split="|"){
  ## tax is a matrix or data.frame with the table of taxonomic abundances, rows are taxa, columns are samples
  ## metadat is an optional data.frame of specimen metadata, rows are samples, columns are variables
  ## if simplenames=TRUE, use only the most detailed level of taxa names in the final object
  ## if roundtointeger=TRUE, values will be rounded to the nearest integer
  xnames = rownames(tax)
  shortnames = gsub(paste0(".+\\", split), "", xnames)
  if(simplenames){
    rownames(tax) = shortnames
  }
  if(roundtointeger){
    tax = round(tax * 1e4)
  }
  x2 = strsplit(xnames, split=split, fixed=TRUE)
  taxmat = matrix(NA, ncol=max(sapply(x2, length)), nrow=length(x2))
  colnames(taxmat) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")[1:ncol(taxmat)]
  rownames(taxmat) = rownames(tax)
  for (i in 1:nrow(taxmat)){
    taxmat[i, 1:length(x2[[i]])] <- x2[[i]]
  }
  taxmat = gsub("[a-z]__", "", taxmat)
  taxmat = phyloseq::tax_table(taxmat)
  otutab = phyloseq::otu_table(tax, taxa_are_rows=TRUE)
  if(is.null(metadat)){
    res = phyloseq::phyloseq(taxmat, otutab)
  }else{
    res = phyloseq::phyloseq(taxmat, otutab, phyloseq::sample_data(metadat))
  }
  return(res)
}


# Given a phyloseq object, return a dataframe of the sample metadata 
# From: https://jacobrprice.github.io/2017/08/26/phyloseq-to-vegan-and-back.html
pssd2veg <- function(physeq) {
  sd <- sample_data(physeq)
  return(as(sd,"data.frame"))
}

# Function to compare the linear modeling framework to an equivalent MaAsLin2 analysis.
compare_lm_maaslin2 <- function(lm_df, maaslin2_df){
  # get significant results from maaslin2
  x <- maaslin2_df$results %>% filter(pval < 0.05 & qval < 0.25) 
  # remove to only columns we need
  y <- data.frame(taxa = x$feature, factor = x$metadata)
  # remove duplicate entries for a given factor
  z <- y[!duplicated(y[c(1,2)]),]
  
  
  # Number unique to Maaslin2
  # Number unique to LM
  # Number in common
  # rootstock
  output_r <- data.frame(MaAsLin2 = length(setdiff(z[z$factor == "rootstock",]$taxa, lm_df[lm_df$factor == "rootstock",]$taxa)),
                         LM_framework = length(setdiff(lm_df[lm_df$factor == "rootstock",]$taxa, z[z$factor == "rootstock",]$taxa)),
                         Conserved = length(intersect(lm_df[lm_df$factor == "rootstock",]$taxa, z[z$factor == "rootstock",]$taxa)))
  # scion
  output_s <- data.frame(MaAsLin2 = length(setdiff(z[z$factor == "scion",]$taxa, lm_df[lm_df$factor == "scion",]$taxa)),
                         LM_framework = length(setdiff(lm_df[lm_df$factor == "scion",]$taxa, z[z$factor == "scion",]$taxa)),
                         Conserved = length(intersect(lm_df[lm_df$factor == "scion",]$taxa, z[z$factor == "scion",]$taxa)))
  # rootstock by scion
  output_rxs <- data.frame(MaAsLin2 = length(setdiff(z[z$factor == "RxS",]$taxa, lm_df[lm_df$factor == "rootstock:scion",]$taxa)),
                           LM_framework = length(setdiff(lm_df[lm_df$factor == "rootstock:scion",]$taxa, z[z$factor == "RxS",]$taxa)),
                           Conserved = length(intersect(lm_df[lm_df$factor == "rootstock:scion",]$taxa, z[z$factor == "RxS",]$taxa)))
  
  output<- rbind(output_r, output_s, output_rxs)
  rownames(output) <- c("Rootstock", "Scion", "RxS")
  
  return(output)
}

# Function to plot comparison table generated from the compare_lm_maaslin2 Function (above).
plot_comparison <- function(compare_lm_maaslin2_df){
  x <- cbind(data.frame(Factor = rownames(compare_lm_maaslin2_df)), compare_lm_maaslin2_df)
  #print(x)
  y <- reshape2::melt(x, id=c("Factor"))
  y$Factor <- factor(y$Factor, levels = c("Rootstock", "Scion", "RxS"))
  
  #print(y)
  ggplot(y, aes(x = variable, y = value, fill = Factor)) +
    geom_bar(stat="identity", color = "black", position=position_dodge()) +
    ylab("Differentailly Abundant Taxa") +
    xlab(NULL) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}





#### 1. Load data ####
metaphlan_file <- read.csv("DATA/metaphlan_taxonomic_profiles.tsv", sep = "\t", header = TRUE, row.names = 1)
# Load and add metadata
metadata_df <- read.csv("DATA/METADATA.tsv", sep = "\t", header = TRUE)
to_factor <- c("rootstock", "scion")
metadata_df[to_factor] <- lapply(metadata_df[to_factor], factor) 
metadata_df$date_collect <- mdy(metadata_df$date_collect)
metadata_df$year <- as.factor(year(metadata_df$date_collect))
rownames(metadata_df) <- metadata_df$Sample_name
summary(metadata_df) #Sanity check
### Reformating and loading into a phyloseq object
# Fix taxonomy strings (inspired by https://github.com/waldronlab/presentations/blob/master/Waldron_2016-06-07_EPIC/metaphlanToPhyloseq.R)
tax_str  <- rownames(metaphlan_file) # vector of taxonomy strings
tax_str2 <- strsplit(tax_str, split='|', fixed=TRUE) # Lists for each taxon
taxmat   <- matrix(NA, ncol=max(sapply(tax_str2, length)), nrow=length(tax_str2)) # Create matrix to hold elements of tax_str2 lists
colnames(taxmat) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")[1:ncol(taxmat)] # col names
rownames(taxmat) <- rownames(metaphlan_file)
for (i in 1:nrow(taxmat)){taxmat[i, 1:length(tax_str2[[i]])] <- tax_str2[[i]]} # Loop to iterate through tax_str2 extracting cells into taxmat
taxmat <- gsub("[a-z]__", "", taxmat)
metaphlan_corrected_df <- cbind(taxmat, metaphlan_file) # Combine the taxonomy matrix to the metaphlan file with relative abundances.

# Because the metaphlan normalizes to relative abundance within a taxonomic level
# I need to extract just the lowest level (aka where we have species assignments). 
# From there I can use phyloseq functions like taxa_glom() to do any subsetting to 
# taxonomic levels needed.

# Take only those assigned to species level
spec_lvl <- metaphlan_corrected_df %>% filter(!is.na(Species))
colSums(as.data.frame(sapply(spec_lvl[,-c(1:7)] , as.numeric))) # Sanity check, Should all be roughly 100% (give or take a decimal)
spec_lvl_tidy <- spec_lvl[,-c(1:7)] # cut off first seven columns (taxa matrix)
colnames(spec_lvl_tidy) <- gsub("_?[N]?_taxonomic_profile", "", colnames(spec_lvl_tidy))
rownames_saved <- rownames(spec_lvl_tidy)
spec_lvl_tidy <- as.data.frame(sapply(spec_lvl_tidy, as.numeric)) #convert all columns to numeric
rownames(spec_lvl_tidy) <- rownames_saved
# Select at the lowest level
phylo_df <- metaphlanToPhyloseq(spec_lvl_tidy, metadat = metadata_df, simplenames = TRUE)
# Relevel the rootstocks to group by pedigree
phylo_df@sam_data$rootstock <- factor(phylo_df@sam_data$rootstock, 
                                      levels = c("101-14 MGT", "3309 Couderc", "Schwarzmann",
                                                 "110 Richter", "1103 Paulsen", "140Ru","775 Paulsen",
                                                 "420A",  "Kober 5BB", "Teleki 5C"))
phylo_df # sanity 563 taxa by 60 samples
rm(metaphlan_corrected_df, metaphlan_file, spec_lvl, spec_lvl_tidy, tax_str2, taxmat)

# Create dataframes for modeling
phylo_df_phy <- psmelt(phyloseq::tax_glom(phylo_df, "Phylum")) # 13 taxa
phylo_df_cla <- psmelt(phyloseq::tax_glom(phylo_df, "Class")) # 29 taxa
phylo_df_ord <- psmelt(phyloseq::tax_glom(phylo_df, "Order")) # 60 taxa
phylo_df_fam <- psmelt(phyloseq::tax_glom(phylo_df, "Family")) # 104 taxa
phylo_df_gen <- psmelt(phyloseq::tax_glom(phylo_df, "Genus")) # 203 taxa
phylo_df_spe <- psmelt(phyloseq::tax_glom(phylo_df, "Species")) # 563 taxa


#### 2. MaAsLin2 #### 
# MaAsLin2 is comprehensive R package for efficiently determining multivariable association
# between phenotypes, environments, exposures, covariates and microbial meta’omic features
library("Maaslin2")
# add interaction term rootstock*scion
metadata_df$RxS <- paste(metadata_df$rootstock, metadata_df$scion, sep = "_")
metadata_df$RxS <- as.factor(metadata_df$RxS)

# Function to generate Maaslin2 input format from phyloseq melt object
Create_Maaslin2_inputdata<- function(phyloseq_melt_obj, taxonomic_level){
  temp <- phyloseq_melt_obj %>%
            select(Sample, Abundance, all_of(taxonomic_level)) %>%
            pivot_wider(names_from = taxonomic_level, values_from = Abundance) %>%
            as.data.frame()
  rownames(temp) <- temp$Sample
  temp$Sample <- NULL
  return(temp)
}

#phylum
phylum_input <- Create_Maaslin2_inputdata(phylo_df_phy, "Phylum")
phyl_maaslin2 <- Maaslin2(input_data = phylum_input, input_metadata =  metadata_df, output = 'Maaslin2/MaasLin2_output_phylum', transform = "NONE",
                          fixed_effects = c('rootstock', 'scion', 'RxS'),
                          random_effects = c('year'),
                          normalization = 'NONE',
                          standardize = FALSE,
                          min_abundance = 0.001, 
                          min_prevalence = 0.25)
#class
class_input <- Create_Maaslin2_inputdata(phylo_df_cla, "Class")
clas_maaslin2 <- Maaslin2(input_data = class_input, input_metadata =  metadata_df, output = 'Maaslin2/MaasLin2_output_class', transform = "NONE",
                          fixed_effects = c('rootstock', 'scion', 'RxS'),
                          random_effects = c('year'),
                          normalization = 'NONE',
                          standardize = FALSE,
                          min_abundance = 0.001, 
                          min_prevalence = 0.25)
#order
order_input <- Create_Maaslin2_inputdata(phylo_df_ord, "Order")
orde_maaslin2 <- Maaslin2(input_data = order_input, input_metadata =  metadata_df, output = 'Maaslin2/MaasLin2_output_order', transform = "NONE",
                          fixed_effects = c('rootstock', 'scion', 'RxS'),
                          random_effects = c('year'),
                          normalization = 'NONE',
                          standardize = FALSE,
                          min_abundance = 0.001, 
                          min_prevalence = 0.25)
#family
family_input <- Create_Maaslin2_inputdata(phylo_df_fam, "Family")
fami_maaslin2 <- Maaslin2(input_data = family_input, input_metadata =  metadata_df, output = 'Maaslin2/MaasLin2_output_family', transform = "NONE",
                          fixed_effects = c('rootstock', 'scion', 'RxS'),
                          random_effects = c('year'),
                          normalization = 'NONE',
                          standardize = FALSE,
                          min_abundance = 0.001, 
                          min_prevalence = 0.25)
#genus
genus_input <- Create_Maaslin2_inputdata(phylo_df_gen, "Genus")
genu_maaslin2 <- Maaslin2(input_data = genus_input, input_metadata =  metadata_df, output = 'Maaslin2/MaasLin2_output_genus', transform = "NONE",
                          fixed_effects = c('rootstock', 'scion', 'RxS'),
                          random_effects = c('year'),
                          normalization = 'NONE',
                          standardize = FALSE,
                          min_abundance = 0.001, 
                          min_prevalence = 0.25)
#species
species_input <- Create_Maaslin2_inputdata(phylo_df_spe, "Species")
spec_maaslin2 <- Maaslin2(input_data = species_input, input_metadata =  metadata_df, output = 'Maaslin2/MaasLin2_output_species', transform = "NONE",
                     fixed_effects = c('rootstock', 'scion', 'RxS'),
                     random_effects = c('year'),
                     normalization = 'NONE',
                     standardize = FALSE,
                     min_abundance = 0.001, 
                     min_prevalence = 0.25)


#### 3. Maaslin2 results ####
# Function to get significant results from maaslin
Sig_table_maaslin2 <- function(maaslin_obj){
  X <- maaslin_obj$results %>%
  filter(pval < 0.05 & qval < 0.05) %>%
  select(taxa = feature, effect = metadata) %>%
  distinct() %>%
  arrange(effect)
  return(X)
}
# Table 2 results
Sig_table_maaslin2(phyl_maaslin2)
Sig_table_maaslin2(clas_maaslin2)
Sig_table_maaslin2(orde_maaslin2)
Sig_table_maaslin2(fami_maaslin2)
Sig_table_maaslin2(genu_maaslin2)
Sig_table_maaslin2(spec_maaslin2)

# Heatmap and selected accession plots
What_to_plot <- spec_maaslin2$results %>%
  filter(pval < 0.05 & qval < 0.05) %>%
  select(Taxa = feature, Factor = metadata, value, coef, pval, qval)

What_to_plot %>%
  select(Taxa) %>%
  distinct()

# Heatmap at species level
a <- phylo_df_spe %>%
  filter(Species %in% c("Devosia_insulae", "Pseudomonas_moraviensis", "Bradyrhizobium_ottawaense", "Bacillus_simplex")) %>%
  mutate(RxS = paste(rootstock, scion, sep = "_")) %>%
  mutate(scion = fct_recode(scion, "Cabernet Sauvignon" = "cabernet sauvignon", Chardonnay = "chardonnay")) %>%
  group_by(Species, rootstock, scion) %>%
  summarize(`Mean Relative Abundance` = mean(Abundance)) 
a <- ggplot(a, aes(y = Species, x = rootstock, fill = `Mean Relative Abundance`)) +
  geom_tile() +
  geom_tile(data = a[a$Species == "Bradyrhizobium_ottawaense",], fill = NA, color = "black", size = 2) +
  scale_fill_viridis_c(limits = c(0, 10), direction = -1) +
  facet_wrap(~scion) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.y = element_text(margin = margin(r = 10)), legend.position = "top", axis.title.x = element_blank())

b <- phylo_df_spe %>%
  filter(Species == "Bradyrhizobium_ottawaense") %>%
  mutate(scion = fct_recode(scion, "Cabernet Sauvignon" = "cabernet sauvignon", Chardonnay = "chardonnay")) %>%{
    ggplot(., aes(x = rootstock, y=Abundance, fill = rootstock)) +
      geom_jitter(width = 0.25, color = "black") +  
      geom_boxplot(outlier.shape = NA, alpha = 0.80, lwd=0.5, color = "black") +
      scale_fill_manual(name="Rootstock", values = rootstock_palette) +
      theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position = "none",
            axis.title.y = element_text(margin = margin(r = 10))) +
      facet_grid(~scion) +
      stat_compare_means(method = "anova", label.x = 7)+        # Add global annova p-value
      labs(x=NULL, y=expression(paste('Relative Abundance', "  ", italic("Bradyrhizobium ottawaense"))))}

Supplement_Tax_maaslin2 <- ggarrange(a,b, labels = "AUTO")

ggsave("figures/Supplemental_Taxonomy_Maaslin2.svg", Supplement_Tax_maaslin2, height = 120, width = 340, units = "mm")

# Some stats for the plots
# Relative abundance
phylo_df_spe %>%
  filter(Species %in% c("Devosia_insulae", "Pseudomonas_moraviensis", "Bradyrhizobium_ottawaense", "Bacillus_simplex")) %>%
  group_by(Species) %>%
  summarize(`Mean Relative Abundance` = mean(Abundance)) 
# pairwise contrasts
spec_maaslin2$results[spec_maaslin2$results$feature == 'Bradyrhizobium_ottawaense',]
spec_maaslin2$results[spec_maaslin2$results$feature == 'Devosia_insulae',]
genu_maaslin2$results[genu_maaslin2$results$feature == 'Nakamurella',]
genu_maaslin2$results[genu_maaslin2$results$feature == 'Paludisphaera',]
# lms and posthocs
lms_spe <- phylo_df_spe %>%
  filter(Species %in% c("Devosia_insulae", "Pseudomonas_moraviensis",
                        "Bradyrhizobium_ottawaense", "Bacillus_simplex")) %>%
  nest_by(Species) %>%
  mutate(mod = list(lmer(Abundance ~ rootstock * scion + (1|year), data=data)))

#Devosia_insulae
anova(lms_spe$mod[[1]])
pairs(emmeans(lms_spe$mod[[1]], ~ rootstock*scion))
#Bradyrhizobium_ottawaense
anova(lms_spe$mod[[2]])
pairs(emmeans(lms_spe$mod[[2]], ~ rootstock*scion))
#Bacillus_simplex
anova(lms_spe$mod[[3]])
pairs(emmeans(lms_spe$mod[[3]], ~ rootstock*scion))
#Pseudomonas_moraviensis
anova(lms_spe$mod[[4]])
pairs(emmeans(lms_spe$mod[[4]], ~ rootstock*scion))


# Relative abundance stats
phylo_df_spe %>%
  filter(Species %in% c("Devosia_insulae", "Pseudomonas_moraviensis",
                        "Bradyrhizobium_ottawaense", "Bacillus_simplex")) %>%
  group_by(Species) %>%
  summarise(mean(Abundance))

phylo_df_spe %>%
  filter(Species == "Bradyrhizobium_ottawaense") %>%
  group_by(rootstock, scion) %>%
  summarise(mean(Abundance))

phylo_df_spe %>%
  filter(Species == "Pseudomonas_moraviensis") %>%
  group_by(rootstock, scion) %>%
  summarise(mean(Abundance))

#### 6. Defining a core microbiome ####
#A taxon's mean relative abundance is log10
#transformed, and then plotted against the proportion of
#discrete samples in which it occurs (with occupancy of 1 to
#be found in all samples). Together, abundance and occu-
#pancy provide rich information for interpreting diversity
#patterns at both population and community levels.
phylo_df_gen

# Determine occupancy and log relative abundance
genus_summary <- phylo_df_gen %>%
  group_by(Genus) %>%
  summarise(OCCUP = mean(Abundance > 0), ABUND = log10(mean(Abundance))) %>%
  arrange(desc(OCCUP))
genus_summary

# Crappy plot
ggplot(genus_summary, aes(x = ABUND, y = OCCUP)) + geom_point() +
  labs(x = "Log10(mean relative abundance)", y = "Occupancy")

# Take 80th precentile abund and genera in all samples
quantile(genus_summary$ABUND, probs = seq(0,1,0.05))
quantile(genus_summary$OCCUP, probs = seq(0,1,0.05))
sum(genus_summary$ABUND >= -0.39131835, na.rm = TRUE) #41 above 80th percentile
sum(genus_summary$OCCUP >= 0.95, na.rm = TRUE) #23 present in 95% samples
CORE <- intersect(which(genus_summary$ABUND >= -0.5594343), which(genus_summary$OCCUP >= 0.95))
genus_summary[CORE,] # list of taxa
length(genus_summary[CORE,]$Genus) # list of taxa

10^-0.39131835

# Plot with core genera colored 
genus_summary$Core <- with(genus_summary, ifelse(ABUND  >= -1.2881450 & OCCUP >= 0.95, "Yes", "No"))
genus_summary$Core <- factor(genus_summary$Core, levels = c("Yes", "No"))

log10_top_95_core <- ggplot(genus_summary, aes(x = ABUND, y = OCCUP, color = Core)) + geom_point() +
  labs(x = "Log10(mean relative abundance)", y = "Occupancy") +
  scale_color_manual(values = c("lightgreen", "firebrick"), labels = c("True", "False")) +
  theme(legend.position = 'right')

ggsave("figures/Core_95_taxa.svg", log10_top_95_core, height = 6, width = 8)

# Make a table for paper
Core_df <- as.data.frame(genus_summary[CORE,])
Core_df$ABUND<- 10^Core_df$ABUND # convert back to relative abundance
Core_df$Core <- NULL
Core_df