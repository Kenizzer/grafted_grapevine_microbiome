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

## Add parentage
metadata_df <- metadata_df %>%
  data.frame() %>%
  mutate(parentage = case_when(
    rootstock %in% c("Teleki 5C", "Kober 5BB", "420A") ~ "V. berlandieri x V. riparia",
    rootstock %in% c("101-14 MGT", "3309 Couderc", "Schwarzmann") ~ "V. riparia x V. rupestris",
    rootstock %in% c("110 Richter", "1103 Paulsen", "140Ru", "775 Paulsen") ~ "V. berlandieri x V. rupestris"))

metadata_df$parentage <- as.factor(metadata_df$parentage)


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
metadata_df$PxS <- paste(metadata_df$parentage, metadata_df$scion, sep = "_")
metadata_df$PxS <- as.factor(metadata_df$PxS)

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
phyl_maaslin2 <- Maaslin2(input_data = phylum_input, input_metadata =  metadata_df, output = 'Maaslin2/MaasLin2_output_parentage_phylum', transform = "NONE",
                          fixed_effects = c('parentage', 'scion', 'PxS'),
                          random_effects = c('year'),
                          normalization = 'NONE',
                          standardize = FALSE,
                          min_abundance = 0.001, 
                          min_prevalence = 0.25)
#class
class_input <- Create_Maaslin2_inputdata(phylo_df_cla, "Class")
clas_maaslin2 <- Maaslin2(input_data = class_input, input_metadata =  metadata_df, output = 'Maaslin2/MaasLin2_output_parentage_class', transform = "NONE",
                          fixed_effects = c('parentage', 'scion', 'PxS'),
                          random_effects = c('year'),
                          normalization = 'NONE',
                          standardize = FALSE,
                          min_abundance = 0.001, 
                          min_prevalence = 0.25)
#order
order_input <- Create_Maaslin2_inputdata(phylo_df_ord, "Order")
orde_maaslin2 <- Maaslin2(input_data = order_input, input_metadata =  metadata_df, output = 'Maaslin2/MaasLin2_output_parentage_order', transform = "NONE",
                          fixed_effects = c('parentage', 'scion', 'PxS'),
                          random_effects = c('year'),
                          normalization = 'NONE',
                          standardize = FALSE,
                          min_abundance = 0.001, 
                          min_prevalence = 0.25)
#family
family_input <- Create_Maaslin2_inputdata(phylo_df_fam, "Family")
fami_maaslin2 <- Maaslin2(input_data = family_input, input_metadata =  metadata_df, output = 'Maaslin2/MaasLin2_output_parentage_family', transform = "NONE",
                          fixed_effects = c('parentage', 'scion', 'PxS'),
                          random_effects = c('year'),
                          normalization = 'NONE',
                          standardize = FALSE,
                          min_abundance = 0.001, 
                          min_prevalence = 0.25)
#genus
genus_input <- Create_Maaslin2_inputdata(phylo_df_gen, "Genus")
genu_maaslin2 <- Maaslin2(input_data = genus_input, input_metadata =  metadata_df, output = 'Maaslin2/MaasLin2_output_parentage_genus', transform = "NONE",
                          fixed_effects = c('parentage', 'scion', 'PxS'),
                          random_effects = c('year'),
                          normalization = 'NONE',
                          standardize = FALSE,
                          min_abundance = 0.001, 
                          min_prevalence = 0.25)
#species
species_input <- Create_Maaslin2_inputdata(phylo_df_spe, "Species")
spec_maaslin2 <- Maaslin2(input_data = species_input, input_metadata =  metadata_df, output = 'Maaslin2/MaasLin2_output_parentage_species', transform = "NONE",
                     fixed_effects = c('parentage', 'scion', 'PxS'),
                     random_effects = c('year'),
                     normalization = 'NONE',
                     standardize = FALSE,
                     min_abundance = 0.001, 
                     min_prevalence = 0.25)