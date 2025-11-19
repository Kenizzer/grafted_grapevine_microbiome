# Analysis of rootstock trial metagenomes
# Code by: Joel F. Swift

#Packages w/ version numbers.
library('tidyverse'); packageVersion('tidyverse')
library('ggpubr'); packageVersion('ggpubr')
library('vegan'); packageVersion('vegan')
library('phyloseq'); packageVersion('phyloseq')
library('lubridate'); packageVersion('lubridate')

#### Color pallettes####
theme_set(theme_bw())
rootstock_palette<- c("101-14 MGT" = "#46436b", "3309 Couderc" = "#7570B3","Schwarzmann" = "#aca9d1",
                      "110 Richter" = "#0a3d2e", "1103 Paulsen" = "#1b9e77",
                      "140Ru" = "#48e0b3",  "775 Paulsen" = "#d3f8ed",
                      "420A" = "#453301", "Kober 5BB" = "#E6AB02", "Teleki 5C" = "#f0cd67")
scion_palette <- c("cabernet sauvignon" = '#ed254e', "chardonnay" = '#0e79b2')
Bacterial_palette <- c("Acidobacteriota" = "#88CCEE", "Actinobacteriota" = "#CC6677",
                             "Bacteroidota" = "#DDCC77", "Chloroflexi" = "#AA4499",
                             "Deinococcota" = "#332288","Firmicutes" = "#117733",
                             "Myxococcota" = "#661100", "Planctomycetota" = "#999933",
                             "Proteobacteria" = "#44AA99", "Verrucomicrobiota" = "#882255",
                             "Desulfobacterota" = "#888888","Crenarchaeota" = "#D55E00",
                             "Other" = "black") 
#### Functions ####
# Plot PCoA using ggplot from phyloseq ordination
PLOT_PCoA <- function(plot_data, distance_matrix, axis1, axis2, split_by){
  if(split_by == 'rootstock'){
    temp <- plot_ordination(plot_data, distance_matrix, axes = c(axis1,axis2))
    temp[["layers"]][[1]][["geom"]][["default_aes"]][["colour"]] <- "white" # make the plot_ordination points white to allow me to use alpha without them showing through.
    Plot <- plot_ordination(plot_data, distance_matrix, axes = c(axis1,axis2)) +
      geom_point(aes(fill=rootstock), shape = 21,  size = 3, alpha = 0.95, color = "black") +
      scale_fill_manual(name = "Rootstock", values = rootstock_palette) +
      labs(color= "Rootstock") +
      xlab(paste("PCoA", axis1, sub(".*\\ ", "", temp$labels$x))) +
      ylab(paste("PCoA", axis2, sub(".*\\ ", "", temp$labels$y))) +
      guides(fill = guide_legend(override.aes = c(shape = 21, alpha = 1)), color = guide_legend(override.aes = c(alpha=1))) +
      theme(legend.position = "right")
    return(Plot)
  } else if(split_by == 'scion'){
    temp <- plot_ordination(plot_data, distance_matrix, axes = c(axis1,axis2))
    temp[["layers"]][[1]][["geom"]][["default_aes"]][["colour"]] <- "white" # make the plot_ordination points white to allow me to use alpha without them showing through.
    Plot <- plot_ordination(plot_data, distance_matrix, axes = c(axis1,axis2)) +
      geom_point(aes(fill=scion), shape = 21, size = 3, alpha = 0.95, color = "black") +
      scale_fill_manual(name = "Scion", values = scion_palette, labels = c("Cabernet Sauvignon", "Chardonnay")) +
      labs(color= "Scion") +
      xlab(paste("PCoA", axis1, sub(".*\\ ", "", temp$labels$x))) +
      ylab(paste("PCoA", axis2, sub(".*\\ ", "", temp$labels$y))) +
      guides(fill = guide_legend(override.aes = c(shape = 21, alpha = 1)), color = guide_legend(override.aes = c(alpha=1))) +
      theme(legend.position = "right")
    return(Plot)
  }
}

# Convert a metaphlan file to a phyloseq object
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


# Function from the R cookbook
# From: http://www.cookbook-r.com/Manipulating_data/Summarizing_data/
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  # Rename the "mean" column
  datac <- rename(datac, c("mean" = measurevar))
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  return(datac)
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

#### 2. PCoA plot  ####
# Calculate bray-curtis
out.bray <- ordinate(phylo_df, method = "MDS", distance = "bray")
# Calculate distance matrixes for use in Vegan
out.dist.bray <- phyloseq::distance(phylo_df, method = "bray")
# get variance explained by axes 1-3
sum(out.bray$values$Eigenvalues[1:3])/sum(out.bray$values$Eigenvalues) #52%
# Plot PCoAs 1x2-3
Bray_C_1_2_scion <- PLOT_PCoA(phylo_df, out.bray, 1, 2, split_by = 'scion') + 
  theme(legend.position="none")
Bray_C_1_2_rootstock <- PLOT_PCoA(phylo_df, out.bray, 1, 2, split_by = 'rootstock') + 
  theme(legend.position="none")

Bray12_pcoa <- ggarrange(Bray_C_1_2_scion,
                         Bray_C_1_2_rootstock,
                         align = 'hv', labels = c("A", "B"), nrow = 2)

# combine with another figure below
#ggsave("figures/Figure1_Bray-Curtis_PCOA_scion_rootstock.svg", Bray12_pcoa, height = 120, width = 85, units = "mm")


##### Adonis/PERMANOVA #####
phylo_metadata <- data.frame(phylo_df@sam_data)
perm <- how(nperm = 1000)
setBlocks(perm) <- with(phylo_metadata, year)
adonis2(out.dist.bray ~ rootstock + scion, data = phylo_metadata, by = "margin", permutations = perm)
adonis2(out.dist.bray ~ rootstock*scion, data = phylo_metadata, by = "margin", permutations = perm)
adonis2(out.dist.bray ~ rootstock*scion + scion + rootstock, data = phylo_metadata, by = "margin", permutations = perm)


#### 3. Taxonomic barplots ####
library(fantaxtic)
get_sample_order_by_taxon <- function(ps_obj, tax_level, taxon_name) {
  psmelt(ps_obj) %>%
    group_by(Sample) %>%
    mutate(Abundance = Abundance / sum(Abundance)) %>%           # Convert to relative abundance per sample
    ungroup() %>%
    filter(.data[[tax_level]] == taxon_name) %>%                 # Filter for the specified taxon at the given level
    group_by(Sample) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop") %>%  # Sum abundance across OTUs of that taxon
    arrange(Abundance) %>%                                       # Sort samples by abundance
    pull(Sample) %>%                                             # Extract sample names
    as.character()                                               # Return as character vector
}


# Make Nested dataframe for plotting
top_nested <- nested_top_taxa(phylo_df,
                              top_tax_level = "Phylum",
                              nested_tax_level = "Family",
                              n_top_taxa = 10, 
                              n_nested_taxa = 3,
                              by_proportion = TRUE) 

# Make single row plots by parentage, will combine later.
Rip_Rup_plots <- plot_nested_bar(subset_samples(top_nested$ps_obj, rootstock %in% c("101-14 MGT", "3309 Couderc", "Schwarzmann")),
                top_level = "Phylum",
                nested_level = "Family",
                palette = Bacterial_palette,         
                legend_title = "Phylum and Family", 
                sample_order = rev(get_sample_order_by_taxon(subset_samples(top_nested$ps_obj, rootstock %in% c("101-14 MGT", "3309 Couderc", "Schwarzmann")), tax_level = "Family", "Rhizobiaceae"))) +
  labs(y = "Relative Abuance") +
  facet_wrap(~rootstock, scales = "free_x") +
  # mostly changes made to fit the pubr theme without breaking legend text formatting from fantaxtic
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold"), legend.key.size = unit(10, "points"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'right',
        strip.background = element_rect(fill = "#F2F2F2", colour = "black", size = 0.7),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.text = element_text(color = "black"),
        strip.text = element_text(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line("black"), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.text.x = element_text(size = 8)) +
  guides(fill=guide_legend(ncol=1))


Rip_Ber_plots <- plot_nested_bar(subset_samples(top_nested$ps_obj, rootstock %in% c("420A", "Kober 5BB", "Teleki 5C")),
                                 top_level = "Phylum",
                                 nested_level = "Family",
                                 palette = Bacterial_palette,         
                                 legend_title = "Phylum and Family",
                                 sample_order = rev(get_sample_order_by_taxon(subset_samples(top_nested$ps_obj, rootstock %in% c("420A", "Kober 5BB", "Teleki 5C")), tax_level = "Family", "Rhizobiaceae"))) +
  labs(y = "Relative Abuance") +
  facet_wrap(~rootstock, scales = "free_x") +
  # mostly changes made to fit the pubr theme without breaking legend text formatting from fantaxtic
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold"), legend.key.size = unit(10, "points"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'right',
        strip.background = element_rect(fill = "#F2F2F2", colour = "black", size = 0.7),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.text = element_text(color = "black"),
        strip.text = element_text(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line("black"), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.text.x = element_text(size = 8)) +
  guides(fill=guide_legend(ncol=1))

Ber_Rup_plots <- plot_nested_bar(subset_samples(top_nested$ps_obj, rootstock %in% c("110 Richter", "1103 Paulsen", "140Ru", "775 Paulsen")),
                                 top_level = "Phylum",
                                 nested_level = "Family",
                                 palette = Bacterial_palette,         
                                 legend_title = "Phylum and Family",
                                 sample_order = rev(get_sample_order_by_taxon(subset_samples(top_nested$ps_obj, rootstock %in% c("110 Richter", "1103 Paulsen", "140Ru", "775 Paulsen")), tax_level = "Family", "Rhizobiaceae"))) +
  
  labs(y = "Relative Abuance") +
  facet_wrap(~rootstock, scales = "free_x", nrow = 1) +
  # mostly changes made to fit the pubr theme without breaking legend text formatting from fantaxtic
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold"), legend.key.size = unit(10, "points"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'right',
        strip.background = element_rect(fill = "#F2F2F2", colour = "black", size = 0.7),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.text = element_text(color = "black"),
        strip.text = element_text(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line("black"), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.text.x = element_text(size = 8)) +
  guides(fill=guide_legend(ncol=1))

# All barplots together
Figure2C <- ggarrange(Rip_Rup_plots, Rip_Ber_plots, Ber_Rup_plots,
                     common.legend = TRUE, legend = "right", nrow = 3, labels = c("C"))
# Taxonomic barplot part of figure two
ggsave("figures/Figure1_PART_Taxabarplots.svg", Figure2C, height = 170, width = 130, units = "mm")

ggsave("figures/Figure1_PART_PCoA.svg", Bray12_pcoa, height = 160, width = 80, units = "mm")

# Put together pcoa and barplots
Figure2ABC <- ggarrange(Bray12_pcoa, Figure2C, widths = c(0.5, 1))
ggsave("figures/Figure1_PCOA_with_Taxabarplots_uneditted.svg", Figure2ABC, height = 170, width = 170, units = "mm")




### Taxa percents, paragraph in results ###
#phyla
phylo_df_melted_phy <- psmelt(phyloseq::tax_glom(phylo_df, "Phylum")) # 13 taxa
summarySE(phylo_df_melted_phy[phylo_df_melted_phy$Phylum == "Proteobacteria",], measurevar = 'Abundance') # mean 71.0 +/- 1.61%
summarySE(phylo_df_melted_phy[phylo_df_melted_phy$Phylum == "Actinobacteria",], measurevar = 'Abundance') # mean 16.6 +/- 1.04%
summarySE(phylo_df_melted_phy[phylo_df_melted_phy$Phylum == "Firmicutes",], measurevar = 'Abundance') # mean 6.7 +/- 0.71%
summarySE(phylo_df_melted_phy[phylo_df_melted_phy$Phylum == "Thaumarchaeota",], measurevar = 'Abundance') # mean 2.4 +/- 0.38%

#Family
phylo_df_melted_fam <- psmelt(phyloseq::tax_glom(phylo_df, "Family")) # 104 taxa
summarySE(phylo_df_melted_fam[phylo_df_melted_fam$Family == "Bradyrhizobiaceae",], measurevar = 'Abundance') # mean 6.67 +/- 0.64%
summarySE(phylo_df_melted_fam[phylo_df_melted_fam$Family == "Rhizobiaceae",], measurevar = 'Abundance') # mean 33.9 +/- 1.61%
summarySE(phylo_df_melted_fam[phylo_df_melted_fam$Family == "Sphingomonadaceae",], measurevar = 'Abundance') # mean 8.8 +/- 0.55%

summarySE(phylo_df_melted_fam[phylo_df_melted_fam$Family == "Micrococcaceae",], measurevar = 'Abundance') # mean 6.67 +/- 0.64%
summarySE(phylo_df_melted_fam[phylo_df_melted_fam$Family == "Mycobacteriaceae",], measurevar = 'Abundance') # mean 3.32 +/- 0.45%
summarySE(phylo_df_melted_fam[phylo_df_melted_fam$Family == "Streptomycetaceae",], measurevar = 'Abundance') # mean 2.66 +/- 0.24%

summarySE(phylo_df_melted_fam[phylo_df_melted_fam$Family == "Bacillaceae",], measurevar = 'Abundance') # mean 6.64 +/- 0.71%
