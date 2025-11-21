# Analysis of rootstock trial metagenomes
# Code by: Joel F. Swift

#Packages w/ version numbers.
library('tidyverse'); packageVersion('tidyverse')
library('ggpubr'); packageVersion('ggpubr')
library('car'); packageVersion('car')
library('lme4'); packageVersion('lme4')
library('lmerTest'); packageVersion('lmerTest')
library('vegan'); packageVersion('vegan')
library('phyloseq'); packageVersion('phyloseq')
library('emmeans'); packageVersion('emmeans')
library('lubridate'); packageVersion('lubridate')
library('RColorBrewer'); packageVersion('RColorBrewer')

#### Color pallettes####
theme_set(theme_pubr())
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
rootstock_palette<- c("101-14 MGT" = "#46436b", "3309 Couderc" = "#7570B3","Schwarzmann" = "#aca9d1",
                      "110 Richter" = "#0a3d2e", "1103 Paulsen" = "#1b9e77",
                      "140Ru" = "#48e0b3",  "775 Paulsen" = "#d3f8ed",
                      "420A" = "#453301", "Kober 5BB" = "#E6AB02", "Teleki 5C" = "#f0cd67")

#### Functions ####
# Function from the R cookbook
# From: http://www.cookbook-r.com/Manipulating_data/Summarizing_data/
## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
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

is_gt <- function(object, dist, threshold){
  samples <- rownames(object)[dist > threshold | dist < -threshold]
  return(samples)
}

# Load data
pathway_file <- read.csv("DATA/pathabundance_relab.tsv", sep = "\t", header = TRUE, row.names = 1)



#### Data Wrangling ####
# Fix pathway strings
# Example BRANCHED-CHAIN-AA-SYN-PWY: superpathway of branched chain amino acid biosynthesis|g__Mesorhizobium.s__Mesorhizobium_sp_LNHC252B00
# What I would like BioCyc ID | Long name | Taxa | measurements ...
path_str  <- rownames(pathway_file) # vector of pathway strings
path_str2 <- strsplit(path_str, split='[|:]+') # Lists for each 

# Checking the number of splits for each string
y <- c()
for (i in 1:length(path_str)){
  x <- lengths(path_str2[i])
  y[[i]] <- x
}
y <- do.call(rbind, y)
summary(as.factor(y))
rm(y, i, x)
# Most of the rownames are split into 3 elements but 5 have only 1 element and 822 have 2 elements
# Looking at these rownames for the cases with 2 it is because they lack an assigned taxon
# The cases of 1 they lack both the long name and assigned taxon
# Will replace with NAs

# Used answer 1 from (https://stackoverflow.com/questions/21103410/irregular-list-of-lists-to-dataframe) to handle this.
#  Find length of each list element
len <- sapply(path_str2, length)
#  Longest gives number of rows
n <- max( len )
#  Number of NAs to fill for column shorter than longest
len <- n - len
path_str3 <- data.frame(t(mapply( function(x,y) c( x , rep( NA , y ) ) , path_str2 , len )))
colnames(path_str3) <- c("BioCyc_ID","Full_name","Taxon")

# Dataframe with first 3 columns representing the pathway string
path_df <- cbind(path_str3, pathway_file)
rm(path_str2, path_str, pathway_file,len, n)

# Fix the 827 without full 3 columns BioCyc ID | Full pathway name | Taxon
# Many of the 827 issues simply lacked the full pathway name (~822) which was readily available via MetaCyc

# Fix 5 with double NAs
path_df[is.na(path_df$Full_name) & is.na(path_df$Taxon),]$BioCyc_ID
# Add these longer pathway names (it seems that most of these are likely deleted from metacyc, atleast from looking at the current release)
# "ARGORNPROST-PWY" L-arginine degradation (Stickland reaction) https://apps.sbri.org/SSGCIDTargetStatus/Pathway/ARGORNPROST-PWY
TEMP <- as.data.frame(path_df[path_df$BioCyc_ID == "ARGORNPROST-PWY",])
TEMP$Taxon <- TEMP$Full_name
TEMP$Full_name <- "L-arginine degradation (Stickland reaction)"
path_df[path_df$BioCyc_ID == "ARGORNPROST-PWY",] <- TEMP
# "PWY-4242"        pantothenate and coenzyme A biosynthesis III http://vm-trypanocyc.toulouse.inra.fr/META/NEW-IMAGE?type=PATHWAY&object=PWY-4242
TEMP <- as.data.frame(path_df[path_df$BioCyc_ID == "PWY-4242",])
TEMP$Taxon <- TEMP$Full_name
TEMP$Full_name <- "pantothenate and coenzyme A biosynthesis III"
path_df[path_df$BioCyc_ID == "PWY-4242",] <- TEMP
# "PWY-5173"        superpathway of acetyl-CoA biosynthesis http://vm-trypanocyc.toulouse.inra.fr/META/NEW-IMAGE?type=PATHWAY&object=PWY-5173
TEMP <- as.data.frame(path_df[path_df$BioCyc_ID == "PWY-5173",])
TEMP$Taxon <- TEMP$Full_name
TEMP$Full_name <- "superpathway of acetyl-CoA biosynthesis"
path_df[path_df$BioCyc_ID == "PWY-5173",] <- TEMP
# "PWY-5791"        1,4-dihydroxy-2-naphthoate biosynthesis II (plants) https://biocyc.org/ARA/NEW-IMAGE?type=PATHWAY&object=PWY-5791
TEMP <- as.data.frame(path_df[path_df$BioCyc_ID == "PWY-5791",])
TEMP$Taxon <- TEMP$Full_name
TEMP$Full_name <- "1,4-dihydroxy-2-naphthoate biosynthesis II (plants)"
path_df[path_df$BioCyc_ID == "PWY-5791",] <- TEMP
# "PWY66-422"       D-galactose degradation V (Leloir pathway) https://biocyc.org/GCF_000294895-HMP/NEW-IMAGE?type=PATHWAY&object=PWY66-422
TEMP <- as.data.frame(path_df[path_df$BioCyc_ID == "PWY66-422",])
TEMP$Taxon <- TEMP$Full_name
TEMP$Full_name <- "D-galactose degradation V (Leloir pathway)"
path_df[path_df$BioCyc_ID == "PWY66-422",] <- TEMP
# Done
path_df[is.na(path_df$Full_name) & is.na(path_df$Taxon),]$BioCyc_ID # Sanity check should return character(0)

# Tidy the dataset 
path_df[is.na(path_df$Taxon),]$Full_name # 493 PATHWAYS (non-species stratified)
Comm_lvl <- path_df[is.na(path_df$Taxon),] # community level (i.e. not stratified to species)
Comm_lvl$Taxon <- NULL # remove unneeded column
colnames(Comm_lvl) <- gsub("_?[N]?_Abundance", "", colnames(Comm_lvl)) # Remove _abundance from the names
Biocyc_IDs_saved <- Comm_lvl$BioCyc_ID
Comm_lvl$BioCyc_ID <- NULL # remove column so I can make the data tidy
TEMP <- Comm_lvl$Full_name 
Comm_lvl <- as.data.frame(t(Comm_lvl)) 
Comm_lvl <- Comm_lvl[-1,] # remove row of pathway names
colnames(Comm_lvl) <- TEMP # add pathway full names to columns

# Super classes of the 493 pathways XXX Not the best, fix this code
Comm_lvl
Biocyc_IDs_saved
Map_file <- read.csv("DATA/map_metacyc-pwy_lineage.tsv", sep = "\t", header = FALSE) # mapping file to superclasses
colnames(Map_file) <- c("Biocyc_ID", "Superclasses") 
superclasses <- Map_file[Map_file$Biocyc_ID %in% Biocyc_IDs_saved, ]

superclasses <- superclasses[!(superclasses$Superclasses=="Super-Pathways"),] # remove duplicates
superclasses <- superclasses[!(superclasses$Superclasses=="Metabolic-Clusters"),] # remove duplicates

# Check pathways with two entries

test <- data.frame(do.call('rbind', strsplit(as.character(superclasses$Superclasses),'|',fixed=TRUE)))
summary(as.factor(test$X1))
summary(as.factor(test[test$X1 == "Biosynthesis", ]$X2)) # break down biosynthesis class to subclasses
summary(as.factor(test[test$X1 == "Degradation", ]$X2))
summary(as.factor(test[test$X1 == "Energy-Metabolism", ]$X2))
test[test$X1 == "Degradation", ]

# Export datasets for WGCNA
saveRDS(Comm_lvl, file = "DATA/Comm_lvl_pathways_WGCNA.rds")

# Adding in the metadata
metadata_df <- read.csv("DATA/METADATA.tsv", sep = "\t", header = TRUE)
to_factor <- c("rootstock", "scion")
metadata_df[to_factor] <- lapply(metadata_df[to_factor], factor) 
metadata_df$date_collect <- mdy(metadata_df$date_collect)
metadata_df$year <- as.factor(year(metadata_df$date_collect))
rownames(metadata_df) <- metadata_df$Sample_name
summary(metadata_df) #Sanity check
# Merge metadata to pathway table
Comm_lvl_w_meta <- merge(Comm_lvl, metadata_df, by = "row.names", all=TRUE)

## Add parentage
metadata_df <- metadata_df %>%
  data.frame() %>%
  mutate(parentage = case_when(
    rootstock %in% c("Teleki 5C", "Kober 5BB", "420A") ~ "V. berlandieri x V. riparia",
    rootstock %in% c("101-14 MGT", "3309 Couderc", "Schwarzmann") ~ "V. riparia x V. rupestris",
    rootstock %in% c("110 Richter", "1103 Paulsen", "140Ru", "775 Paulsen") ~ "V. berlandieri x V. rupestris"))

metadata_df$parentage <- as.factor(metadata_df$parentage)
# add interaction term parentage*scion
metadata_df$PxS <- paste(metadata_df$parentage, metadata_df$scion, sep = "_")
metadata_df$PxS <- as.factor(metadata_df$PxS)
# add interaction term rootstock*scion
metadata_df$RxS <- as.factor(paste(metadata_df$rootstock, metadata_df$scion, sep = "_"))
metadata_df$RxS <- as.factor(metadata_df$RxS)


#### MaAsLin2 #### 
library("Maaslin2")

# Convert to CPM
# From the documentation: We tend to work with CPM units because we find them to
# be more convenient, but they are numerically equivalent to relative abundances
# for modeling purposes (CPM = RA * 1e6).

# Save row names, convert to numeric, replace rownames
saved_rownames <- rownames(Comm_lvl)
Comm_lvl_n <- as.data.frame(sapply(Comm_lvl, as.numeric))
rownames(Comm_lvl_n) <- saved_rownames
Comm_lvl_CPM <- Comm_lvl_n * 1000000


Comm_maaslin2_test <- Maaslin2(input_data = Comm_lvl_CPM, input_metadata =  metadata_df, output = './MaasLin2/MaasLin2_output_pathway', transform = "NONE",
                               fixed_effects = c('rootstock', 'scion', 'RxS'),
                               random_effects = c('year'),
                               normalization = 'NONE',
                               analysis_method = 'LM',
                               cores = 3,
                               standardize = FALSE)

Comm_maaslin2_test <- Maaslin2(input_data = Comm_lvl_CPM, input_metadata =  metadata_df, output = './MaasLin2/MaasLin2_output_pathway_parentage', transform = "NONE",
                               fixed_effects = c('parentage', 'scion', 'PxS'),
                               random_effects = c('year'),
                               normalization = 'NONE',
                               analysis_method = 'LM',
                               cores = 3,
                               standardize = FALSE)


# Replotting some of the significant results a bit better than the MAASLIN2 defaults
# add meta
Comm_lvl_CPM_w_meta <- merge(Comm_lvl_CPM, metadata_df, by = "row.names", all=TRUE)
# Reorder rootstocks to plot with gradient of parentage
Comm_lvl_CPM_w_meta$rootstock <- factor(Comm_lvl_CPM_w_meta$rootstock, levels = c("101-14 MGT", "3309 Couderc", "Schwarzmann",
                                               "110 Richter", "1103 Paulsen", "140Ru",
                                               "420A", "775 Paulsen", "Kober 5BB", "Teleki 5C"))

flavin_plot <- ggplot(Comm_lvl_CPM_w_meta, aes(x=rootstock, y=` flavin biosynthesis I (bacteria and plants)`, fill = rootstock)) +
  geom_jitter(width = 0.25) +
  geom_boxplot(outlier.shape = NA, alpha = 0.80, lwd=1, color = "black") +
  xlab (NULL) + ylab("Flavin biosynthesis I [RIBOSYN2-PWY] (CPM)") +
  scale_fill_manual(name="Rootstock", values = rootstock_palette) +
  theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position = "none",
        axis.title.y = element_text(margin = margin(r = 10))) +
  geom_hline(yintercept = mean(Comm_lvl_CPM_w_meta$` flavin biosynthesis I (bacteria and plants)`), linetype = 3, size = 1, color = "black") + # Add horizontal line at base mean
  annotate("text", x = "Teleki 5C", y = 5200, label = "*", size = 8)

ggsave("figures/flavin_biosyn_I_rootstock.svg", flavin_plot, height = 6, width = 6)
ggsave("figures/flavin_biosyn_I_rootstock.png", flavin_plot, height = 6, width = 6)


# Posthoc testing  
flavin_mod <- lmer(` flavin biosynthesis I (bacteria and plants)` ~ rootstock*scion + (1|year), Comm_lvl_CPM_w_meta)
anova(flavin_mod) # anova, 
# comparisons to basemean
# stats compare to base means
base_mean <- mean(Comm_lvl_CPM_w_meta$` flavin biosynthesis I (bacteria and plants)`)
Comm_lvl_CPM_w_meta %>%
  group_by(rootstock) %>%
  dplyr::summarise(p = t.test(` flavin biosynthesis I (bacteria and plants)`, mu = base_mean)$p.value) %>%
  mutate(p.adj = p.adjust(p, method = "BH"))
# Teleki 5c significant in comparison to basemean

# Stratified flavin plots
stratified_flavin_df <- path_df[path_df$BioCyc_ID == "RIBOSYN2-PWY", ]
stratified_flavin_df <- stratified_flavin_df[!is.na(stratified_flavin_df$Taxon),] # stratified to taxa
colnames(stratified_flavin_df) <- gsub("_?[N]?_Abundance", "", colnames(stratified_flavin_df)) # Remove _abundance from the names
stratified_flavin_df$BioCyc_ID <- NULL # remove column so I can make the data tidy
stratified_flavin_df$Full_name <- NULL # remove column so I can make the data tidy
TEMP <- stratified_flavin_df$Taxon 
stratified_flavin_df <- as.data.frame(t(stratified_flavin_df)) 
stratified_flavin_df <- stratified_flavin_df[-1,] # remove row of pathway names
colnames(stratified_flavin_df) <- TEMP # add pathway full names to columns
saved_rownames <- rownames(stratified_flavin_df)
stratified_flavin_df <- as.data.frame(sapply(stratified_flavin_df, as.numeric))
rownames(stratified_flavin_df) <- saved_rownames
stratified_flavin_CPM <- stratified_flavin_df * 1000000
# add meta
stratified_flavin_CPM_w_meta <- merge(stratified_flavin_CPM, metadata_df, by = "row.names", all=TRUE)

# Reorder rootstocks to plot with gradient of parentage
stratified_flavin_CPM_w_meta$rootstock <- factor(stratified_flavin_CPM_w_meta$rootstock, levels = c("101-14 MGT", "3309 Couderc", "Schwarzmann",
                                                                                  "110 Richter", "1103 Paulsen", "140Ru",
                                                                                  "420A", "775 Paulsen", "Kober 5BB", "Teleki 5C"))
sort(colSums(stratified_flavin_CPM))

ggplot(stratified_flavin_CPM_w_meta, aes(x=rootstock, y=`g__Bradyrhizobium.s__Bradyrhizobium_japonicum`, fill = rootstock)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.90, lwd=0.5, color = "black") +
  geom_jitter(width = 0.25) +
  xlab (NULL) + #ylab("Flavin biosynthesis I [RIBOSYN2-PWY] (CPM)") +
  scale_fill_manual(name="Rootstock", values = rootstock_palette) +
  theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position = "none",
        axis.title.y = element_text(margin = margin(r = 10)))

summarySE(stratified_flavin_CPM_w_meta, measurevar = "unclassified", groupvars = c("rootstock"))
summarySE(stratified_flavin_CPM_w_meta, measurevar = "g__Rhodospirillales_unclassified.s__Rhodospirillales_bacterium_URHD0017", groupvars = c("rootstock"))
summarySE(stratified_flavin_CPM_w_meta, measurevar = "g__Bradyrhizobium.s__Bradyrhizobium_ottawaense", groupvars = c("rootstock"))
summarySE(stratified_flavin_CPM_w_meta, measurevar = "g__Bradyrhizobium.s__Bradyrhizobium_japonicum", groupvars = c("rootstock"))
summarySE(stratified_flavin_CPM_w_meta, measurevar = "g__Variovorax.s__Variovorax_sp_CF079", groupvars = c("rootstock"))
summarySE(stratified_flavin_CPM_w_meta, measurevar = "g__Variovorax.s__Variovorax_sp_WS11", groupvars = c("rootstock"))
summarySE(stratified_flavin_CPM_w_meta, measurevar = "g__Mesorhizobium.s__Mesorhizobium_amorphae", groupvars = c("rootstock"))
summarySE(stratified_flavin_CPM_w_meta, measurevar = "g__Bradyrhizobium.s__Bradyrhizobium_lablabi", groupvars = c("rootstock"))
summarySE(stratified_flavin_CPM_w_meta, measurevar = "g__Sphingomonas.s__Sphingomonas_sp_FARSPH", groupvars = c("rootstock"))
summarySE(stratified_flavin_CPM_w_meta, measurevar = "g__Streptomyces.s__Streptomyces_sp_75", groupvars = c("rootstock"))
summarySE(stratified_flavin_CPM_w_meta, measurevar = "g__Bacillus.s__Bacillus_megaterium", groupvars = c("rootstock"))
summarySE(stratified_flavin_CPM_w_meta, measurevar = "g__Rhodoplanes.s__Rhodoplanes_sp_Z2_YC6860", groupvars = c("rootstock"))
summarySE(stratified_flavin_CPM_w_meta, measurevar = "g__Variovorax.s__Variovorax_paradoxus", groupvars = c("rootstock"))
summarySE(stratified_flavin_CPM_w_meta, measurevar = "g__Mesorhizobium.s__Mesorhizobium_plurifarium", groupvars = c("rootstock"))


# superpathway of N-acetylglucosamine, N-acetylmannosamine and N-acetylneuraminate degradation
Comm_lvl_CPM_w_meta$`superpathway of N-acetylglucosamine, N-acetylmannosamine and N-acetylneuraminate degradation`
# Correct plotting order
levels(Comm_lvl_CPM_w_meta$rootstock) <- names(rootstock_palette)

#plot
N_acetylglucosamine_plot <- ggplot(Comm_lvl_CPM_w_meta, aes(x=rootstock, y=` superpathway of N-acetylglucosamine, N-acetylmannosamine and N-acetylneuraminate degradation`, fill = rootstock)) +
  geom_jitter(width = 0.25) +
  geom_boxplot(outlier.shape = NA, alpha = 0.80, lwd=0.5, color = "black") +
  xlab (NULL) + ylab("GLCMANNANAUT-PWY (CPM)") +
  scale_fill_manual(name="Rootstock", values = rootstock_palette) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position = "none",
        axis.title.y = element_text(margin = margin(r = 10))) +
  facet_wrap(~scion, labeller = labeller(scion = c("cabernet sauvignon" = "Cabernet Sauvignon",
                                                   "chardonnay" = "Chardonnay")))

# Posthoc testing  
N_acetylglucosamine_mod <- lmer(` superpathway of N-acetylglucosamine, N-acetylmannosamine and N-acetylneuraminate degradation` ~ rootstock*scion + (1|year), Comm_lvl_CPM_w_meta)
anova(N_acetylglucosamine_mod) # anova
pairs(emmeans(N_acetylglucosamine_mod, ~ rootstock|scion)) # none signififcant

# stats compare to base means
base_mean <- mean(Comm_lvl_CPM_w_meta$` superpathway of N-acetylglucosamine, N-acetylmannosamine and N-acetylneuraminate degradation`)
Comm_lvl_CPM_w_meta %>%
  group_by(rootstock) %>%
  dplyr::summarise(p = t.test(` superpathway of N-acetylglucosamine, N-acetylmannosamine and N-acetylneuraminate degradation`, mu = base_mean)$p.value) %>%
  mutate(p.adj = p.adjust(p, method = "BH"))


# superpathway.of.pyrimidine.deoxyribonucleotides.de.novo.biosynthesis
Comm_lvl_CPM_w_meta$` superpathway of pyrimidine deoxyribonucleotides de novo biosynthesis`

pyrimidine_plot <- ggplot(Comm_lvl_CPM_w_meta, aes(x=rootstock, y=` superpathway of pyrimidine deoxyribonucleotides de novo biosynthesis`, fill = rootstock)) +
  geom_jitter(width = 0.25) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, lwd=0.5, color = "black") +
  xlab (NULL) + ylab("RIBOSYN2-PWY (CPM)") +
  scale_fill_manual(name="Rootstock", values = rootstock_palette) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position = "none",
        axis.title.y = element_text(margin = margin(r = 10))) +
  facet_wrap(~scion, labeller = labeller(scion = c("cabernet sauvignon" = "Cabernet Sauvignon",
                                                   "chardonnay" = "Chardonnay")))

# Posthoc testing  
pyrimidine_mod <- lmer(` superpathway of pyrimidine deoxyribonucleotides de novo biosynthesis` ~ rootstock*scion + (1|year), Comm_lvl_CPM_w_meta)
anova(pyrimidine_mod) # anova
pairs(emmeans(pyrimidine_mod, ~ rootstock|scion)) # none signififcant
# stats compare to base means
base_mean <- mean(Comm_lvl_CPM_w_meta$` superpathway of pyrimidine deoxyribonucleotides de novo biosynthesis`)
Comm_lvl_CPM_w_meta %>%
  group_by(rootstock) %>%
  dplyr::summarise(p = t.test(` superpathway of pyrimidine deoxyribonucleotides de novo biosynthesis`, mu = base_mean)$p.value) %>%
  mutate(p.adj = p.adjust(p, method = "BH"))


## Combo plot
pathway_RxS_combo_plot <- ggarrange(N_acetylglucosamine_plot, pyrimidine_plot, labels = "AUTO")
ggsave("figures/Pathway_RxS_combo_plot.svg", pathway_RxS_combo_plot, height = 8, width = 18)


#### Defining core pathways ####

# Determine occupancy of pathways
Comm_lvl_CPM 
Occupancy_path <- colSums(Comm_lvl_CPM != 0)/60
# log10 CPM
Log_CPM_path <- log10(colSums(Comm_lvl_CPM)/60)

# What data should look like (i.e. long)
# TAXA Log(CPM) occupancy
length(Occupancy_path) # should be 493
length(Log_CPM_path) # should be 493
Abund_occupancy_df <- data.frame(TAXA = colnames(Comm_lvl_CPM), ABUND = Log_CPM_path, OCCUP = Occupancy_path)
# Crappy plot
ggplot(Abund_occupancy_df, aes(x = ABUND, y = OCCUP)) + geom_point() +
  labs(x = "Log10(mean relative abundance)", y = "Occupancy")

# Take 80th precentile abund and genera in all samples
quantile(Abund_occupancy_df$ABUND, probs = seq(0,1,0.05))
quantile(Abund_occupancy_df$OCCUP, probs = seq(0,1,0.05))
sum(Abund_occupancy_df$ABUND >=  3.6240543, na.rm = TRUE) #25 above 80th percentile 4207 CPM
sum(Abund_occupancy_df$OCCUP >= 0.95, na.rm = TRUE) #280 present in 95% samples
CORE <- intersect(which(Abund_occupancy_df$ABUND >= 3.6240543), which(Abund_occupancy_df$OCCUP >= 0.95))
Abund_occupancy_df[CORE,] # list of taxa
length(Abund_occupancy_df[CORE,]$TAXA) # list of taxa

# Plot with core genera colored 
Abund_occupancy_df$Core <- with(Abund_occupancy_df, ifelse(ABUND  >= 3.6240543 & OCCUP >= 0.95, "Yes", "No"))
Abund_occupancy_df$Core <- factor(Abund_occupancy_df$Core, levels = c("Yes", "No"))

log10_top_95_core <- ggplot(Abund_occupancy_df, aes(x = ABUND, y = OCCUP, color = Core)) + geom_point() +
  labs(x = "Log10(mean CPM)", y = "Occupancy") +
  scale_color_manual(values = c("lightgreen", "firebrick"), labels = c("True", "False")) +
  theme_bw() +
  theme(legend.position = 'right')

ggsave("figures/Core_95_path.svg", log10_top_95_core, height = 6, width = 8)

# Make a table for paper
Core_df <- as.data.frame(Abund_occupancy_df[CORE,])
Core_df$ABUND<- 10^Core_df$ABUND # convert back to CPM
Core_df$Core <- NULL
Core_df

# Relate to superpathway classes
Comm_lvl
Biocyc_IDs_saved
Biocyc_vec <- unique(path_str3[path_str3$Full_name %in% Core_df$TAXA,]$BioCyc_ID) # vector of biocyc IDs
Map_file <- read.csv("DATA/map_metacyc-pwy_lineage.tsv", sep = "\t", header = FALSE) # mapping file to superclasses
colnames(Map_file) <- c("Biocyc_ID", "Superclasses") 
Core_superclasses <- Map_file[Map_file$Biocyc_ID %in% Biocyc_vec, ]
Core_superclasses <- Core_superclasses[!(Core_superclasses$Superclasses=="Super-Pathways"),] # remove duplicates
Core_superclasses <- Core_superclasses[!(Core_superclasses$Superclasses=="Metabolic-Clusters"),] # remove duplicates
length(Core_superclasses$Biocyc_ID)
test <- data.frame(do.call('rbind', strsplit(as.character(Core_superclasses$Superclasses),'|',fixed=TRUE)))
summary(as.factor(test$X1))
summary(as.factor(test[test$X1 == "Biosynthesis", ]$X2)) # break down biosynthesis class to subclasses