# Analysis of rootstock trial metagenomes
# Code by: Joel F. Swift

#Packages w/ version numbers.
library('tidyverse'); packageVersion('tidyverse')
library('rebus'); packageVersion('rebus')
library('viridis'); packageVersion('viridis')
library('car'); packageVersion('car')
library('corrplot'); packageVersion('corrplot')
library('ggpubr'); packageVersion('ggpubr')
library('factoextra'); packageVersion('factoextra')
library('ggthemes'); packageVersion('ggthemes')
library('lmerTest'); packageVersion('lmerTest')
library('lubridate'); packageVersion('lubridate')
library('emmeans'); packageVersion('emmeans')


#### Color pallettes####
theme_set(theme_pubr())
rootstock_palette <- c("101-14 MGT" = "#46436b", "3309 Couderc" = "#7570B3","Schwarzmann" = "#aca9d1",
                      "110 Richter" = "#0a3d2e", "1103 Paulsen" = "#1b9e77", "140Ru" = "#48e0b3",  "775 Paulsen" = "#d3f8ed",
                      "420A" = "#453301", "Kober 5BB" = "#E6AB02", "Teleki 5C" = "#f0cd67")
scion_palette <- c('#ed254e', '#0e79b2')

# Functions
is_gt <- function(object, dist, threshold){
  samples <- rownames(object)[dist > threshold | dist < -threshold]
  return(samples)
}

#### 1. Reading and tidying data ####
# Read data
X <- read.csv("DATA/Ionomics_data_cleaned.csv", header = TRUE)
# Make things factors that should be factors
X$plant_body_site <- as.factor(X$plant_body_site)
X$scion <- as.factor(X$scion)
X$rootstock <- as.factor(X$rootstock)
X$host_virus_stat <- as.factor(X$host_virus_stat)
X$block <- as.factor(X$block)
X$date_collect <- mdy(X$date_collect)
X$year <- as.factor(year(X$date_collect))

# Remove berries and samples from rootstocks not metagenomically sequenced
Root_only_ions <- X[X$plant_body_site == "root",]
Rootstocks_to_keep <- c("101-14 MGT", "3309 Couderc","Schwarzmann", "110 Richter",
                        "1103 Paulsen", "140Ru", "420A", "775 Paulsen","Kober 5BB", "Teleki 5C")
Select_rootstock_ions <- filter(Root_only_ions, Root_only_ions$rootstock %in% Rootstocks_to_keep)
Select_rootstock_ions$rootstock

# Set all negative ion values to zero
Select_rootstock_ions$B[Select_rootstock_ions$B <0] <- 0
Select_rootstock_ions$Na[Select_rootstock_ions$Na <0] <- 0
Select_rootstock_ions$Mg[Select_rootstock_ions$Mg <0] <- 0
Select_rootstock_ions$Al[Select_rootstock_ions$Al <0] <- 0
Select_rootstock_ions$P[Select_rootstock_ions$P <0] <- 0
Select_rootstock_ions$K[Select_rootstock_ions$K <0] <- 0
Select_rootstock_ions$Ca[Select_rootstock_ions$Ca <0] <- 0
Select_rootstock_ions$Fe[Select_rootstock_ions$Fe <0] <- 0
Select_rootstock_ions$Mn[Select_rootstock_ions$Mn <0] <- 0
Select_rootstock_ions$Co[Select_rootstock_ions$Co <0] <- 0
Select_rootstock_ions$Ni[Select_rootstock_ions$Ni <0] <- 0
Select_rootstock_ions$Cu[Select_rootstock_ions$Cu <0] <- 0
Select_rootstock_ions$Zn[Select_rootstock_ions$Zn <0] <- 0
Select_rootstock_ions$As[Select_rootstock_ions$As <0] <- 0
Select_rootstock_ions$Se[Select_rootstock_ions$Se <0] <- 0
Select_rootstock_ions$Rb[Select_rootstock_ions$Rb <0] <- 0
Select_rootstock_ions$Sr[Select_rootstock_ions$Sr <0] <- 0
Select_rootstock_ions$Mo[Select_rootstock_ions$Mo <0] <- 0
Select_rootstock_ions$Cd[Select_rootstock_ions$Cd <0] <- 0

#### 2. Modeling ####
#### 2018 analysis
model_B <- lmer(scale(B) ~ rootstock * scion + (1|block) + (1|year), Select_rootstock_ions)
model_P <- lmer(scale(P) ~ rootstock * scion + (1|block) + (1|year), Select_rootstock_ions)
model_K <- lmer(scale(K) ~ rootstock * scion + (1|block) + (1|year), Select_rootstock_ions)
model_Na <- lmer(scale(Na) ~ rootstock * scion + (1|block) + (1|year), Select_rootstock_ions)
model_Mg <- lmer(scale(Mg) ~ rootstock * scion + (1|block) + (1|year), Select_rootstock_ions)
model_Al <- lmer(scale(Al) ~ rootstock * scion + (1|block) + (1|year), Select_rootstock_ions)
model_Ca <- lmer(scale(Ca) ~ rootstock * scion + (1|block) + (1|year), Select_rootstock_ions)
model_Fe <- lmer(scale(Fe) ~ rootstock * scion + (1|block) + (1|year), Select_rootstock_ions)
model_Mn <- lmer(scale(Mn) ~ rootstock * scion + (1|block) + (1|year), Select_rootstock_ions)
model_Co <- lmer(scale(Co) ~ rootstock * scion + (1|block) + (1|year), Select_rootstock_ions)
model_Ni <- lmer(scale(Ni) ~ rootstock * scion + (1|block) + (1|year), Select_rootstock_ions)
model_Cu <- lmer(scale(Cu) ~ rootstock * scion + (1|block) + (1|year), Select_rootstock_ions)
model_Zn <- lmer(scale(Zn) ~ rootstock * scion + (1|block) + (1|year), Select_rootstock_ions)
model_As <- lmer(scale(As) ~ rootstock * scion + (1|block) + (1|year), Select_rootstock_ions)
model_Se <- lmer(scale(Se) ~ rootstock * scion + (1|block) + (1|year), Select_rootstock_ions)
model_Rb <- lmer(scale(Rb) ~ rootstock * scion + (1|block) + (1|year), Select_rootstock_ions)
model_Sr <- lmer(scale(Sr) ~ rootstock * scion + (1|block) + (1|year), Select_rootstock_ions)
model_Mo <- lmer(scale(Mo) ~ rootstock * scion + (1|block) + (1|year), Select_rootstock_ions)
model_Cd <- lmer(scale(Cd) ~ rootstock * scion + (1|block) + (1|year), Select_rootstock_ions)


#### Heat map of LM results ####
Return_SS_proportion <- function(model){
  SS <- c(anova(model)$Sum, sum(resid(model)^2))
  output <- (SS / sum(SS)) * 100
  output <- as.data.frame(output) 
  rownames(output) <- c(rownames(anova(model)), "Residuals")
  output <- t(output)
  return(as.data.frame(output))
}
# Make list of ion names in order and run the function to generate the % variance
Ion_list <- c("B", "P", "K", "Na", "Mg", "Al", "Ca", "Fe", "Mn", "Co", "Ni", "Cu", "Zn", "As", "Se", "Rb", "Sr", "Mo", "Cd")
Percent_var <- bind_rows(Return_SS_proportion(model_B), Return_SS_proportion(model_P), Return_SS_proportion(model_K), Return_SS_proportion(model_Na), Return_SS_proportion(model_Mg), Return_SS_proportion(model_Al), Return_SS_proportion(model_Ca), Return_SS_proportion(model_Fe), Return_SS_proportion(model_Mn), Return_SS_proportion(model_Co), Return_SS_proportion(model_Ni), Return_SS_proportion(model_Cu), Return_SS_proportion(model_Zn), Return_SS_proportion(model_As), Return_SS_proportion(model_Se), Return_SS_proportion(model_Rb), Return_SS_proportion(model_Sr), Return_SS_proportion(model_Mo), Return_SS_proportion(model_Cd))
rownames(Percent_var) <- Ion_list
# remove residual 
Percent_var_noint_nores <- Percent_var[, !names(Percent_var) %in% "Residuals"]
# Make dataframe of P values to connect to the % var df above
Get_Pval<- function(model){
  Pval<- as.data.frame(anova(model)$`Pr(>F)`)
  rownames(Pval) <- paste(rownames(anova(model)), "p", sep = "_")
  Pval<- t(Pval)
  return(as.data.frame(Pval))
}
# Run function above on all anova and bind output by row
Pval <- bind_rows(Get_Pval(model_B), Get_Pval(model_P), Get_Pval(model_K), Get_Pval(model_Na), Get_Pval(model_Mg), Get_Pval(model_Al), Get_Pval(model_Ca), Get_Pval(model_Fe), Get_Pval(model_Mn), Get_Pval(model_Co), Get_Pval(model_Ni), Get_Pval(model_Cu), Get_Pval(model_Zn), Get_Pval(model_As), Get_Pval(model_Se), Get_Pval(model_Rb), Get_Pval(model_Sr), Get_Pval(model_Mo), Get_Pval(model_Cd))
rownames(Pval) <- Ion_list
# Combine
SS_pval_combo_df <- cbind(Percent_var_noint_nores,Pval)
SS_pval_combo_df$element <- rownames(SS_pval_combo_df)
# Get variance columns
total_var <- SS_pval_combo_df %>% dplyr::select(element, "rootstock":"rootstock:scion")
# Reorganize and rename to format for plotting
total_var <- total_var %>% gather(key=factor, value=var, -element)
total_var <- total_var %>% mutate(factor=str_replace(factor, "_var" %R% END, ""))
# Get p-value columns
total_p_name <- SS_pval_combo_df %>% dplyr::select("rootstock_p":"rootstock:scion_p")
total_p <- data.frame(t(apply(total_p_name, 1, FUN=p.adjust, method='fdr')))
colnames(total_p)<- colnames(total_p_name) #This is hacky but I need to preserve the colnames to correctly join them later on.
total_p$element <- SS_pval_combo_df$element
# Reorganize and rename to format for plotting
total_p <- total_p %>% gather(key=factor, value=p_value, -element)
total_p <- total_p %>% mutate(factor=str_replace(factor,"_p" %R% END, ""))
# Join variance and p-value tables back together 
total_var_p <- full_join(total_var, total_p,by=c("element", "factor"))
# Only  significant p_values
total_var_p_sig <- total_var_p %>% filter(p_value < 0.05) #Only Molybdenum showed a significant effect of rootstock:scion
# postHoc
pairs(emmeans(model_Mo, ~ rootstock|scion))
anova(model_Mo)
# clean env
rm(total_p, total_var_p, total_p_name, total_var, SS_pval_combo_df, Pval, Percent_var_noint_nores, Percent_var, Ion_list)
rm(model_B, model_P, model_K, model_Na, model_Mg, model_Al, model_Ca, 
   model_Fe, model_Mn, model_Co, model_Ni, model_Cu, model_Zn, model_As,
   model_Se, model_Rb, model_Sr, model_Mo, model_Cd)


#### 3. Correlation Matrix ####
library('corrplot'); packageVersion("corrplot")
library('Hmisc')
# Genertate matrix and change P values for 'Same ion x Same ion' comparisons from NA to 1
Correlation_matrix <- rcorr(as.matrix(Select_rootstock_ions %>% dplyr::select(B:Cd)))
Correlation_matrix$P[is.na(Correlation_matrix$P)] <- 1
# Save
svg(file = "figures/Correlation_matrix_ionomics.svg",  height = 6, width = 12)
corrplot(Correlation_matrix$r, type= "upper", method = "color", p.mat = Correlation_matrix$P, sig.level = 0.05, insig = "label_sig", outline = TRUE, mar=c(1,1,1,1))
dev.off()


#### 4. PCA ####
d <- Select_rootstock_ions %>% dplyr::select(B:Cd)
rownames(d) <- Select_rootstock_ions$Sample
d <- as.data.frame(scale(d, scale=T, center=T)) # scale n center
x <- apply(d, 2, is_gt, object=d, threshold=5) # remove samples 5 std from the mean
outliers <- unique(unlist(x))
data_clean <- Select_rootstock_ions[!(Select_rootstock_ions$Sample %in% outliers),]
d_clean <- data_clean %>% dplyr::select(B:Cd)
d_cleanScale <- scale(d_clean, scale=T, center=T)
pca <- prcomp(d_cleanScale)
x <- summary(pca)
x # summary of PC axes
# Add PCs to metadata
d <- as.data.frame(pca$x) %>% dplyr::select(PC1:PC19)
d$Sample <- data_clean$Sample
d$rootstock <- data_clean$rootstock
d$scion <- data_clean$scion
d$year <- as.factor(data_clean$year)

# Regular PCA plots
roots_PC12 <- ggplot(d, aes(x=PC1, y=PC2, fill=rootstock)) +
  geom_point(size=5, shape = 21, alpha = 0.95, color = "black")  +
  scale_fill_manual(name = "Rootstock", values=rootstock_palette) +
  xlab("PC1 (38.8%)") + ylab("PC2 (19.94%)") +
  theme(legend.position = 'right')
scion_PC12 <- ggplot(d, aes(x=PC1, y=PC2, fill=scion)) +
  geom_point(size=5, shape = 21, alpha = 0.95, color = "black")  +
  scale_fill_manual(name = "Scion", values=scion_palette) +
  xlab("PC1 (38.8%)") + ylab("PC2 (19.94%)") +
  theme(legend.position = 'right')
ggarrange(roots_PC12, scion_PC12)

# Split biplots
roots_PC12_biplot <- fviz_pca_biplot(pca, repel = TRUE, geom = "point", fill.ind = d$rootstock, pointshape = 21,
                pointsize = 5, alpha.ind = 0.95, col.var = "black",
                label = "var", palette = rootstock_palette, invisible = "quali", title = NULL) +
                labs(color = "Rootstock", fill = "Rootstock") + xlab("PC1 (38.8%)") +
                ylab("PC2 (19.9%)") + theme_pubr(legend = 'right')
scion_PC12_biplot <- fviz_pca_biplot(pca, repel = TRUE, geom = "point", fill.ind = d$scion, pointshape = 21,
                pointsize = 5, alpha.ind = 0.95, col.var = "black",
                label = "var", palette = scion_palette, invisible = "quali", title = NULL) +
                labs(color = "Scion", fill = "Scion") + xlab("PC1 (38.8%)") +
                ylab("PC2 (19.9%)") + theme_pubr(legend = 'right')
root_scion_biplot <- ggarrange(roots_PC12_biplot, scion_PC12_biplot, labels = "AUTO")
ggsave("figures/Biplot_rootstock_and_scion.svg", root_scion_biplot, height = 8, width = 16)

# Loadings
sort(pca$rotation[,1]) # loadings PC1
sort(pca$rotation[,2]) # loadings PC2

# Linear models for PCs
PC01_MOD<- lmer(PC1 ~ rootstock*scion + (1|year), d)
PC02_MOD<- lmer(PC2 ~ rootstock*scion + (1|year), d)
PC03_MOD<- lmer(PC3 ~ rootstock*scion + (1|year), d)
PC04_MOD<- lmer(PC4 ~ rootstock*scion + (1|year), d)
PC05_MOD<- lmer(PC5 ~ rootstock*scion + (1|year), d)
PC06_MOD<- lmer(PC6 ~ rootstock*scion + (1|year), d)
PC07_MOD<- lmer(PC7 ~ rootstock*scion + (1|year), d)
PC08_MOD<- lmer(PC8 ~ rootstock*scion + (1|year), d)
PC09_MOD<- lmer(PC9 ~ rootstock*scion + (1|year), d)
PC10_MOD<- lmer(PC10 ~ rootstock*scion + (1|year), d)

# Make list of ion names in order and run the function to generate the % variance
PC_list <- c("PC01", "PC02", "PC03", "PC04", "PC05", "PC06", "PC07", "PC08", "PC09", "PC10")
Percent_var <- bind_rows(Return_SS_proportion(PC01_MOD), 
                         Return_SS_proportion(PC02_MOD), 
                         Return_SS_proportion(PC03_MOD), 
                         Return_SS_proportion(PC04_MOD), 
                         Return_SS_proportion(PC05_MOD), 
                         Return_SS_proportion(PC06_MOD), 
                         Return_SS_proportion(PC07_MOD), 
                         Return_SS_proportion(PC08_MOD), 
                         Return_SS_proportion(PC09_MOD), 
                         Return_SS_proportion(PC10_MOD))
rownames(Percent_var) <- PC_list
# remove residual 
Percent_var_noint_nores <- Percent_var[, !names(Percent_var) %in% "Residuals"]
# Make dataframe of P values to connect to the % var df above
Pval <- bind_rows(Get_Pval(PC01_MOD), 
                  Get_Pval(PC02_MOD), 
                  Get_Pval(PC03_MOD), 
                  Get_Pval(PC04_MOD), 
                  Get_Pval(PC05_MOD), 
                  Get_Pval(PC06_MOD), 
                  Get_Pval(PC07_MOD), 
                  Get_Pval(PC08_MOD), 
                  Get_Pval(PC09_MOD), 
                  Get_Pval(PC10_MOD))
rownames(Pval) <- PC_list
# Combine
SS_pval_combo_df <- cbind(Percent_var_noint_nores,Pval)
SS_pval_combo_df$PC <- rownames(SS_pval_combo_df)
# Get variance columns
total_var <- SS_pval_combo_df %>% dplyr::select(PC, "rootstock":"rootstock:scion")
# Reorganize and rename to format for plotting
total_var <- total_var %>% gather(key=factor, value=var, -PC)
total_var <- total_var %>% mutate(factor=str_replace(factor, "_var" %R% END, ""))
# Get p-value columns
total_p_name <- SS_pval_combo_df %>% dplyr::select("rootstock_p":"rootstock:scion_p")
total_p <- data.frame(t(apply(total_p_name, 1, FUN=p.adjust, method='BH')))
colnames(total_p)<- colnames(total_p_name) #This is hacky but I need to preserve the colnames to correctly join them later on.
total_p$PC <- SS_pval_combo_df$PC
# Reorganize and rename to format for plotting
total_p <- total_p %>% gather(key=factor, value=p_value, -PC)
total_p <- total_p %>% mutate(factor=str_replace(factor,"_p" %R% END, ""))
# Join variance and p-value tables back together 
total_var_p <- full_join(total_var, total_p,by=c("PC", "factor"))
# Only  significant p_values
total_var_p_sig <- total_var_p %>% filter(p_value < 0.05) 
total_var_p_sig

anova(PC06_MOD)
anova(PC08_MOD)

# reorder rootstocks for plots 
d$rootstock <- factor(d$rootstock, levels = names(rootstock_palette))

# Split biplots for PCs that siginficant affects of rootstock and scion
roots_PC18_biplot <- fviz_pca_biplot(pca, axes = c(1,8), repel = TRUE, geom = "point", fill.ind = d$rootstock, pointshape = 21,
                                     pointsize = 5, alpha.ind = 0.95, col.var = "black",
                                     label = "var", palette = rootstock_palette, invisible = "quali", title = NULL) +
                                     labs(fill = "Rootstock") + xlab("PC1 (38.8%)") +
                                     ylab("PC8 (2.22%)") + theme_pubr(legend = 'right')
scion_PC16_biplot <- fviz_pca_biplot(pca, axes = c(1,6), repel = TRUE, geom = "point", fill.ind = d$scion, pointshape = 21,
                                     pointsize = 5, alpha.ind = 0.95, col.var = "black",
                                     label = "var", palette = scion_palette, invisible = "quali", title = NULL) +
                                     labs(fill = "Scion") + xlab("PC1 (38.8%)") +
                                     ylab("PC6 (4.68%)") + theme_pubr(legend = 'right')

root_scion_biplot_selected <- ggarrange(roots_PC18_biplot, scion_PC16_biplot, labels = "AUTO")
ggsave("figures/Biplot_rootstock_and_scion_selected_PCs.svg", root_scion_biplot_selected, height = 8, width = 16)


# Loadings
sort(pca$rotation[,8]) # loadings PC8
sort(pca$rotation[,6]) # loadings PC6

# Test if PCs are correlated 
cor(pca$x[,8], pca$x[,6], method = "pearson")
cor.test(pca$x[,8], pca$x[,6], method = "pearson") # They are not correlated


##### 5. Linear modeling of PC6/8 with  MaasLin2 #####
# Load metadata
metadata_df <- read.csv("DATA/METADATA.tsv", sep = "\t", header = TRUE)
to_factor <- c("rootstock", "scion")
metadata_df[to_factor] <- lapply(metadata_df[to_factor], factor) 
metadata_df$date_collect <- mdy(metadata_df$date_collect)
metadata_df$year <- as.factor(year(metadata_df$date_collect))
rownames(metadata_df) <- metadata_df$Sample_name

# Remove berries and samples from rootstocks not metagenomically sequenced
Root_only_ions <- metadata_df[metadata_df$plant_body_site == "root",]
Rootstocks_to_keep <- c("101-14 MGT", "3309 Couderc","Schwarzmann", "110 Richter",
                        "1103 Paulsen", "140Ru", "420A", "775 Paulsen","Kober 5BB", "Teleki 5C")
Select_rootstock_ions <- filter(Root_only_ions, Root_only_ions$rootstock %in% Rootstocks_to_keep)
metadata_df <- Select_rootstock_ions
summary(metadata_df) #Sanity check
metadata_df$rootstock <- factor(metadata_df$rootstock, levels = names(rootstock_palette))


# Genera
Taxa_df <- readRDS("DATA/Genus_lvl_taxa_WGCNA.rds")


# Pathways
Path_df <- readRDS("DATA/Comm_lvl_pathways_WGCNA.rds")
# Convert to CPM
# From the documentation: We tend to work with CPM units because we find them to
# be more convenient, but they are numerically equivalent to relative abundances
# for modeling purposes (CPM = RA * 1e6).
# Save row names, convert to numeric, replace rownames
saved_rownames <- rownames(Path_df)
Comm_lvl_n <- as.data.frame(sapply(Path_df, as.numeric))
rownames(Comm_lvl_n) <- saved_rownames
Comm_lvl_CPM <- Comm_lvl_n * 1000000


library("Maaslin2")
# add interaction term rootstock*scion
metadata_df$RxS <- paste(metadata_df$rootstock, metadata_df$scion, sep = "_")
# Add PC10s to metadata
metadata_df$PC1 <- pca$x[,1]
metadata_df$PC2 <- pca$x[,2]
metadata_df$PC3 <- pca$x[,3]
metadata_df$PC4 <- pca$x[,4]
metadata_df$PC5 <- pca$x[,5]
metadata_df$PC6 <- pca$x[,6]
metadata_df$PC7 <- pca$x[,7]
metadata_df$PC8 <- pca$x[,8]
metadata_df$PC9 <- pca$x[,9]
metadata_df$PC10 <- pca$x[,10]


#genra
taxa_maaslin2 <- Maaslin2(input_data = Taxa_df, input_metadata =  metadata_df, output = './MaasLin2/MaasLin2_output_genera_PCs', transform = "NONE",
                          fixed_effects = c('PC6', 'PC8'),
                          random_effects = c('year'),
                          normalization = 'NONE',
                          min_prevalence = 0.6,
                          min_abundance = 0.5,
                          standardize = FALSE)

#pathways
path_maaslin2 <- Maaslin2(input_data = Comm_lvl_CPM, input_metadata =  metadata_df, output = './MaasLin2/MaasLin2_output_pathways_PCs', transform = "NONE",
                          fixed_effects = c('PC6', 'PC8'),
                          normalization = 'NONE',
                          min_prevalence = 0.5,
                          standardize = FALSE)


# Plots of significant ones
abund <- Taxa_df %>% 
  rownames_to_column(var = "Sample_name") %>% 
  select(Sample_name, Mesorhizobium)

Meso_df <- merge(abund, metadata_df, by = "Sample_name")
summary(lm(Mesorhizobium~PC6, data=Meso_df))
Meso_PC6_plot <- ggplot(Meso_df, aes(x = PC6, y = Mesorhizobium)) +
  geom_smooth(method = 'lm') +
  geom_point(size = 3, shape = 21, alpha = 0.95, aes(fill = scion)) +
  scale_fill_manual(values = scion_palette, name = "Scion", labels = c("Cabernet Sauvignon", "Chardonnay")) +
  annotate("text",x = 2, y = 13, label = "R^2: 0.14", parse = TRUE, size = 5) +
  labs(y = expression(paste("Relative abundance", italic("  Mesorhizobium"))), x = "PC6 (4.68%)") +
  coord_cartesian(ylim=c(0,13)) +
  theme(legend.position = 'right')

Strep_df <- data.frame( Streptomyces = Taxa_df$Streptomyces, PC8 = pca$x[,8], Rootstock = metadata_df$rootstock)
summary(lm(Streptomyces~PC8, data=Strep_df))
Strep_PC8_plot <- ggplot(Strep_df, aes(x = PC8, y = Streptomyces)) +
  geom_smooth(method = 'lm') +
  geom_point(size = 3, shape = 21, alpha = 0.95, aes(fill = Rootstock)) +
  scale_fill_manual(values = rootstock_palette) +
  annotate("text",x = 1, y = 13, label = "R^2: 0.20", parse = TRUE,  size = 5) +
  labs(y = expression(paste("Relative abundance", italic("  Streptomyces"))), x = "PC8 (2.22%)") +
  coord_cartesian(ylim=c(0,13)) +
  theme(legend.position = 'right')

taxa_assocaited_with_elements_plot <- ggarrange(Strep_PC8_plot, Meso_PC6_plot, labels = "AUTO")

ggsave("figures/taxa_assocaited_with_elements.svg", taxa_assocaited_with_elements_plot, height = 8, width = 12)

# pathways
ARGSYN_df <- data.frame( ARGSYN_PWY  = Comm_lvl_CPM$` L-arginine biosynthesis I (via L-ornithine)` , PC6 = pca$x[,6], Scion = metadata_df$scion)
summary(lm(ARGSYN_PWY~PC6, data=ARGSYN_df))
ARGSYN_PC6_plot <- ggplot(ARGSYN_df, aes(x = PC6, y = ARGSYN_PWY)) +
  geom_smooth(method = 'lm') +
  geom_point(size = 3, shape = 21, alpha = 0.95, aes(fill = Scion)) +
  scale_fill_manual(values = scion_palette) +
  annotate("text",x = 2, y = 8500, label = "R^2: 0.15", parse = TRUE, size = 5) +
  labs(y = "ARGSYN-PWY CPM", x = "PC6 (4.68%)") +
  theme(legend.position = 'right')


PWY_7400_df <- data.frame(PWY_7400 = Comm_lvl_CPM$` L-arginine biosynthesis IV (archaebacteria)` , PC6 = pca$x[,6], Scion = metadata_df$scion)
summary(lm(PWY_7400~PC6, data=PWY_7400_df))
PWY_7400_PC6_plot <- ggplot(PWY_7400_df, aes(x = PC6, y = PWY_7400)) +
  geom_smooth(method = 'lm') +
  geom_point(size = 3, shape = 21, alpha = 0.95, aes(fill = Scion)) +
  scale_fill_manual(values = scion_palette) +
  annotate("text",x = 2, y = 8500, label = "R^2: 0.16", parse = TRUE, size = 5) +
  labs(y = "PWY-7400 CPM", x = "PC6 (4.68%)") +
  theme(legend.position = 'right')

Path_assocaited_with_elements_plot <- ggarrange(ARGSYN_PC6_plot, PWY_7400_PC6_plot, labels = "AUTO")

ggsave("figures/Path_assocaited_with_elements_plot.svg", Path_assocaited_with_elements_plot, height = 8, width = 12)

## Testing other PCS
#genra
taxa_maaslin2_allPCs <- Maaslin2(input_data = Taxa_df, input_metadata =  metadata_df, output = './MaasLin2/MaasLin2_output_genera_ALLPCs', transform = "NONE",
                          fixed_effects = c('PC1', 'PC2','PC3', 'PC4','PC5','PC7','PC9', 'PC10'),
                          random_effects = c('year'),
                          normalization = 'NONE',
                          min_prevalence = 0.6,
                          min_abundance = 0.5,
                          standardize = FALSE,
                          plot_scatter = FALSE)

#pathways

Comm_lvl_CPM_clean <- janitor::clean_names(Comm_lvl_CPM)

path_maaslin2_allPCs <- Maaslin2(input_data = Comm_lvl_CPM_clean, input_metadata =  metadata_df, output = './MaasLin2/MaasLin2_output_pathways_ALLPCs', transform = "NONE",
                          fixed_effects = c('PC1', 'PC2','PC3', 'PC4','PC5','PC7','PC9', 'PC10'),
                          normalization = 'NONE',
                          min_prevalence = 0.5,
                          standardize = FALSE,
                          plot_scatter = FALSE)


# Genera
# Make a plot for all the PCs
sig_genera <- taxa_maaslin2_allPCs$results %>% filter(pval < 0.05)

# ensure the columns are treated correctly
df_genra <- sig_genera %>%
  mutate(PC = factor(name, levels = paste0("PC", 1:10)), feature = factor(feature)) %>% 
  mutate(signed_logq = -log10(qval) * sign(coef))

# create the heatmap
genra_plot <- ggplot(df_genra, aes(x = PC, y = feature, fill = signed_logq)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = expression(-log[10](qval) * sign(coef))) +
  theme_classic(base_size = 12) +
  theme( axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) + 
  labs(x = "Principal Component", y = "Genera")
  

# Pathways
# Make a plot for all the PCs
sig_path <- path_maaslin2_allPCs$results %>% filter(qval < 0.05)

# ensure the columns are treated correctly
df_path <- sig_path %>%
  mutate(PC = factor(name, levels = paste0("PC", 1:10)), feature = factor(feature)) %>% 
  mutate(signed_logq = -log10(qval) * sign(coef))

# create the heatmap
pathway_plot <- ggplot(df_path, aes(x = PC, y = feature, fill = signed_logq)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = expression(-log[10](qval) * sign(coef))) +
  theme_classic(base_size = 12) +
  theme( axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) + 
  labs(x = "Principal Component", y = "Pathway")

ggarrange(pathway_plot, align = "h")

ggsave("figures/Supplemental_maaslin2_results.svg", width = 600, height = 500, unit = "mm")
