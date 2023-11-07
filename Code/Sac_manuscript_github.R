#Loading all necessary libraries
library(nortest)
library(tidyverse)
library(rcompanion)
library(phyloseq)
library(corrplot)
library(geosphere)
library(ANCOMBC)
library(microbiome)
library(MicrobiotaProcess)
library(ggplot2)
library(ggrepel)

###Section 1: Assessing statistical assumptions of data, transforming for correlations analyses, and correlational analyses

setwd("~/OneDrive/RStudio") #setting directory to that which contains necessary files
all_taxa <- readRDS(file = "SA_all_taxa") #raw phyloseq object
all_taxa_meta <- sample_data(all_taxa)
all_taxa_meta <- as_tibble(all_taxa_meta)
TRCs_tab <- all_taxa_meta %>% 
  dplyr::select(B1, HMP, cHET, HET, AmMP) %>% 
  mutate_if(is_character, as.numeric)
rownames(TRCs_tab) <- all_taxa_meta$sample

for (column in TRCs_tab) {
  result_1 <- shapiro.test(column)
  result_2 <- ad.test(column)
  print(result_1)
  print(result_2)
} #Testing for normality with two methods (Shapiro and Anderson-Darling tests)

tukey_transform <- function(a){
  rcompanion::transformTukey(a, plotit = F, quiet = F, statistic = 2)
}

x.trans.norm <- function(x.trans) {
  (x.trans-min(x.trans))/(max(x.trans)-min(x.trans))
}

##transforming did not fix normality so let's check for outliers 
ggplot(TRCs_tab, aes(x = B1)) +geom_histogram(binwidth = .1) #remove >20
ggplot(TRCs_tab, aes(x = HMP)) +geom_histogram(binwidth = 1) #remove >200
ggplot(TRCs_tab, aes(x = cHET)) +geom_histogram(binwidth = 50) #remove >1000
ggplot(TRCs_tab, aes(x = HET)) +geom_histogram(binwidth = 10) #drop HET from the analysis because of all the zeros
ggplot(TRCs_tab, aes(x = AmMP)) +geom_histogram(binwidth = 5) #remove > 1000

phylo_plastid <- subset_taxa(phylo_all_taxa_rare, 
                             Kingdom == "Eukaryota" | Phylum == "Cyanobacteria")
phylo_archaea <- subset_taxa(phylo_all_taxa_rare, 
                             Kingdom == "\tArchaea")
phylo_bacteria <- subset_taxa(phylo_all_taxa_rare, 
                              Kingdom == "\tBacteria" & Phylum != "Cyanobacteria")

plastid_div <- microbiome::diversity(phylo_plastid, index = "shannon")
archaea_div <- microbiome::diversity(phylo_archaea, index = "shannon")
bacteria_div <- microbiome::diversity(phylo_bacteria, index = "shannon")

bacteria_richness <- colSums(data.frame(otu_table(phylo_bacteria)))
archaea_richness <- colSums(data.frame(otu_table(phylo_archaea)))
plastid_richness <- colSums(data.frame(otu_table(phylo_plastid)))

all_taxa_meta$"bacteria_div" <- bacteria_div$shannon
all_taxa_meta$"archaea_div" <- archaea_div$shannon
all_taxa_meta$"algae_div" <- plastid_div$shannon
all_taxa_meta$"bacteria_rich" <- richness$bacteria
all_taxa_meta$"archaea_rich" <- richness$archaea
all_taxa_meta$"algae_rich" <- richness$plastid
all_taxa_meta <- as_tibble(all_taxa_meta)

##Diffabund taxa (values taken from section 2)
#sed = hyporheic zone
all_taxa_meta$"Positive_ANCOM_BC_ASVs_b1" <- diffabund_pos_sums_sed_b1
all_taxa_meta$"Negative_ANCOM_BC_ASVs_hmp" <- diffabund_neg_sums_sed_hmp
all_taxa_meta$"Positive_ANCOM_BC_ASVs_chet" <- diffabund_pos_sums_sed_chet
all_taxa_meta$"Negative_ANCOM_BC_ASVs_het" <- diffabund_neg_sums_sed_het
all_taxa_meta$"Positive_ANCOM_BC_ASVs_ammp" <- diffabund_pos_sums_sed_ammp
all_taxa_meta$"Negative_ANCOM_BC_ASVs_b1" <- diffabund_neg_sums_sed_b1
all_taxa_meta$"Positive_ANCOM_BC_ASVs_hmp" <- diffabund_pos_sums_sed_hmp
all_taxa_meta$"Negative_ANCOM_BC_ASVs_chet" <- diffabund_neg_sums_sed_chet
all_taxa_meta$"Positive_ANCOM_BC_ASVs_het" <- diffabund_pos_sums_sed_het
all_taxa_meta$"Negative_ANCOM_BC_ASVs_ammp" <- diffabund_neg_sums_sed_ammp
#

#wat = surface water
all_taxa_meta$"Positive_ANCOM_BC_ASVs_b1_wat" <- diffabund_pos_sums_wat_b1
all_taxa_meta$"Negative_ANCOM_BC_ASVs_hmp_wat" <- diffabund_neg_sums_wat_hmp
all_taxa_meta$"Positive_ANCOM_BC_ASVs_chet_wat" <- diffabund_pos_sums_wat_chet
all_taxa_meta$"Negative_ANCOM_BC_ASVs_het_wat" <- diffabund_neg_sums_wat_het
all_taxa_meta$"Positive_ANCOM_BC_ASVs_ammp_wat" <- diffabund_pos_sums_wat_ammp
all_taxa_meta$"Negative_ANCOM_BC_ASVs_b1_wat" <- diffabund_neg_sums_wat_b1
all_taxa_meta$"Negative_ANCOM_BC_ASVs_chet_wat" <- diffabund_neg_sums_wat_chet
all_taxa_meta$"Positive_ANCOM_BC_ASVs_het_wat" <- diffabund_pos_sums_wat_het
all_taxa_meta$"Negative_ANCOM_BC_ASVs_ammp_wat" <- diffabund_neg_sums_wat_ammp
#

TRCs_tab <- all_taxa_meta %>% 
  dplyr::select(B1, HMP, cHET, HET, AmMP, pH, ,Chl, DO_conc, Temp, Turb,
                Positive_ANCOM_BC_ASVs_b1, Negative_ANCOM_BC_ASVs_b1, 
                Positive_ANCOM_BC_ASVs_chet, Negative_ANCOM_BC_ASVs_chet,
                Positive_ANCOM_BC_ASVs_hmp, Negative_ANCOM_BC_ASVs_hmp,
                Positive_ANCOM_BC_ASVs_het, Negative_ANCOM_BC_ASVs_het,
                Positive_ANCOM_BC_ASVs_ammp, Negative_ANCOM_BC_ASVs_ammp,
                Positive_ANCOM_BC_ASVs_b1_wat, Negative_ANCOM_BC_ASVs_b1_wat, 
                Positive_ANCOM_BC_ASVs_chet_wat, Negative_ANCOM_BC_ASVs_chet_wat,
                Negative_ANCOM_BC_ASVs_hmp_wat, 
                Positive_ANCOM_BC_ASVs_het_wat, Negative_ANCOM_BC_ASVs_het_wat,
                Positive_ANCOM_BC_ASVs_ammp_wat,
                algae_rich, algae_div, archaea_rich, archaea_div, bacteria_rich, bacteria_div) %>% 
  dplyr::mutate_if(is_character, as.numeric) %>% 
  dplyr::mutate(Type = all_taxa_meta$Type) %>% 
  dplyr::mutate(across(.cols = B1:bacteria_div, .fns = tukey_transform), 
                across(.cols = B1:bacteria_div, .fns = x.trans.norm)) %>% 
  filter(B1 < 20, 
         HMP < 200, 
         cHET < 1000)


sed_tab <- TRCs_tab %>% 
  filter(Type == "Sediment") %>% 
  dplyr::select(-Positive_ANCOM_BC_ASVs_b1_wat, -Negative_ANCOM_BC_ASVs_b1_wat, 
                -Positive_ANCOM_BC_ASVs_chet_wat, -Negative_ANCOM_BC_ASVs_chet_wat,
                -Negative_ANCOM_BC_ASVs_hmp_wat, 
                -Positive_ANCOM_BC_ASVs_het_wat, -Negative_ANCOM_BC_ASVs_het_wat,
                -Positive_ANCOM_BC_ASVs_ammp_wat)

wat_tab <- TRCs_tab %>% 
  filter(Type == "Water") %>% 
  dplyr::select(-Positive_ANCOM_BC_ASVs_b1, -Negative_ANCOM_BC_ASVs_b1, 
                -Positive_ANCOM_BC_ASVs_chet, -Negative_ANCOM_BC_ASVs_chet,
                -Negative_ANCOM_BC_ASVs_hmp, -Positive_ANCOM_BC_ASVs_hmp,
                -Negative_ANCOM_BC_ASVs_het, -Positive_ANCOM_BC_ASVs_het,
                -Negative_ANCOM_BC_ASVs_ammp, -Positive_ANCOM_BC_ASVs_ammp)

for (column in TRCs_tab) {
  result_1 <- shapiro.test(column)
  result_2 <- ad.test(column)
  print(result_1)
  print(result_2)
}

sed_tab <- sed_tab %>% 
  dplyr::rename("[DO]" = DO_conc, 
                "+B1 Diffabund" = Positive_ANCOM_BC_ASVs_b1, 
                "-B1 Diffabund" = Negative_ANCOM_BC_ASVs_b1, 
                "+HMP Diffabund" = Positive_ANCOM_BC_ASVs_hmp,
                "-HMP Diffabund" = Negative_ANCOM_BC_ASVs_hmp, 
                "+cHET Diffabund" = Positive_ANCOM_BC_ASVs_chet, 
                "-cHET Diffabund" = Negative_ANCOM_BC_ASVs_chet, 
                "+HET Diffabund" = Positive_ANCOM_BC_ASVs_het, 
                "-HET Diffabund" = Negative_ANCOM_BC_ASVs_het, 
                "+AmMP Diffabund" = Positive_ANCOM_BC_ASVs_ammp, 
                "-AmMP Diffabund" = Negative_ANCOM_BC_ASVs_ammp,
                "Algal Richness" = algae_rich, 
                "Algal Diversity" = algae_div, 
                "Turbidity" = Turb, 
                "Bacterial Richness" = bacteria_rich, 
                "Bacterial Diversity" = bacteria_div, 
                "Archaeal Richness" = archaea_rich, 
                "Archaeal Diversity" = archaea_div)
wat_tab <- wat_tab %>% 
  dplyr::rename("[DO]" = DO_conc, 
                "+B1 Diffabund" = Positive_ANCOM_BC_ASVs_b1_wat, 
                "-B1 Diffabund" = Negative_ANCOM_BC_ASVs_b1_wat, 
                "-HMP Diffabund" = Negative_ANCOM_BC_ASVs_hmp_wat, 
                "+cHET Diffabund" = Positive_ANCOM_BC_ASVs_chet_wat, 
                "-cHET Diffabund" = Negative_ANCOM_BC_ASVs_chet_wat, 
                "+HET Diffabund" = Positive_ANCOM_BC_ASVs_het_wat, 
                "-HET Diffabund" = Negative_ANCOM_BC_ASVs_het_wat, 
                "+AmMP Diffabund" = Positive_ANCOM_BC_ASVs_ammp_wat, 
                "Algal Richness" = algae_rich, 
                "Algal Diversity" = algae_div, 
                "Turbidity" = Turb, 
                "Bacterial Richness" = bacteria_rich, 
                "Bacterial Diversity" = bacteria_div, 
                "Archaeal Richness" = archaea_rich, 
                "Archaeal Diversity" = archaea_div)

sed_tab <- sed_tab %>% 
  dplyr::select(-Type, -Turbidity, -pH, -Chl, -'[DO]', -Temp)

sed_tab <- sed_tab[,var_order]

wat_tab <- wat_tab %>% 
  dplyr::select(-Type)

wat_tab <- wat_tab[,var_order_wat]

set.seed(1252)
cor_trc_sed <- cor(sed_tab, method = "spearman")
cor_trc_wat <- cor(wat_tab, method = "spearman")

testRes_sed <- cor.mtest(cor_trc_sed, conf.level = .95)
testRes_wat <- cor.mtest(cor_trc_wat, conf.level = .95)


var_order <- c("B1", "HMP", "cHET", "HET", "AmMP", "Algal Richness", "Algal Diversity", 
               "Bacterial Richness", "Bacterial Diversity", "Archaeal Richness", "Archaeal Diversity",
               "+B1 Diffabund", "-B1 Diffabund", "+HMP Diffabund",
               "-HMP Diffabund", "+cHET Diffabund", "-cHET Diffabund", 
               "+HET Diffabund", "-HET Diffabund", 
               "+AmMP Diffabund", "-AmMP Diffabund")

var_order_wat <- c("B1", "HMP", "cHET", "HET", "AmMP", "pH", "Chl", "[DO]", "Temp", "Turbidity", "Algal Richness", "Algal Diversity", 
                   "Bacterial Richness", "Bacterial Diversity", "Archaeal Richness", "Archaeal Diversity", "+B1 Diffabund", "-B1 Diffabund", 
                   "-HMP Diffabund", "+cHET Diffabund", "-cHET Diffabund", 
                   "+HET Diffabund", "-HET Diffabund", 
                   "+AmMP Diffabund")

corrplot(cor_trc_sed, method = "color", type = "upper", diag = FALSE, 
         outline = "black", order = "original", 
         rect.col = "black", tl.col = "black", tl.srt = 90, 
         pch.cex = 1, sig.level = c(0.001, 0.01, 0.05),
         insig = 'label_sig', p.mat = testRes_sed$p, tl.cex = 1)

corrplot(cor_trc_wat, method = "color", type = "upper", diag = FALSE, 
         outline = "black", order = "original", 
         rect.col = "black", tl.col = "black", tl.srt = 90, 
         pch.cex = 1, sig.level = c(0.001, 0.01, 0.05),
         insig = 'label_sig', p.mat = testRes_wat$p)

##Mantel tests

mantel_meta$"Lat" <- c(39.80095, 39.80095, 39.80095, 39.80095, 39.80095, 39.80095, 39.80095, 39.80095, 
                       40.49651, 40.49651, 40.49651, 40.49651, 40.49651, 40.49651, 40.49651, 40.49651,
                       40.58355, 40.58355, 40.58355, 40.58355, 40.58355, 40.58355, 40.58355, 40.58355, 
                       39.51092, 39.51092, 39.49980, 39.49980, 
                       39.53529, 39.53529, 39.53529, 39.53529, 39.53529, 39.53529, 39.53529, 39.53529, 
                       39.53529, 39.53529, 39.53529, 39.53529, 39.53529, 39.53529, 39.53529, 39.53529,
                       40.42725, 40.42725, 40.42725, 40.42725, 40.42725, 40.42725, 40.42725, 40.42725, 
                       40.65063, 40.65063, 40.65063, 40.65063, 40.65063, 40.65063, 40.65063, 40.65063, 
                       40.41119, 40.41119, 40.41119, 40.41119, 40.41119, 40.41119, 40.41119, 40.41119,
                       40.75474, 40.75474, 40.75474, 40.75474, 40.75474, 40.75474, 40.75474, 40.75474)
  
mantel_meta$"Long" <- c(-121.65774, -121.65774, -121.65774, -121.65774, -121.65774, -121.65774, -121.65774, -121.65774, 
                        -122.49741, -122.49741, -122.49741, -122.49741, -122.49741, -122.49741, -122.49741, -122.49741,
                        -122.54020, -122.54020, -122.54020, -122.54020, -122.54020, -122.54020, -122.54020, -122.54020, 
                        -121.58048, -121.58048, -121.60829, -121.60829, 
                        -121.55164, -121.55164, -121.55164, -121.55164, -121.55164, -121.55164, -121.55164, -121.55164, 
                        -121.55164, -121.55164, -121.55164, -121.55164, -121.55164, -121.55164, -121.55164, -121.55164, 
                        -121.98701, -121.98701, -121.98701, -121.98701, -121.98701, -121.98701, -121.98701, -121.98701,
                        -122.36914, -122.36914, -122.36914, -122.36914, -122.36914, -122.36914, -122.36914, -122.36914,
                        -121.97952, -121.97952, -121.97952, -121.97952, -121.97952, -121.97952, -121.97952, -121.97952, 
                        -122.35815, -122.35815, -122.35815, -122.35815, -122.35815, -122.35815, -122.35815, -122.35815)


mantel_meta <- TRCs_tab %>% 
  select(B1, AmMP, cHET, HMP, HET, Lat, Long) %>% 
  mutate(across(.cols = B1:HET, .fns = tukey_transform)) %>% 
  mutate(Type = all_taxa_meta$Type)

TRCs_tab_sed <- mantel_meta %>% 
  filter(Type == "Sediment")
TRCs_tab_wat <- mantel_meta %>% 
  filter(Type == "Water")

phylo_tot_sed <- subset_samples(all_taxa_rare, Type == "Sediment") #Rarefied phyloseq object
phylo_tot_wat <- subset_samples(all_taxa_rare, Type == "Water")

mantel_taxa_sed <- data.frame(tax_table(phylo_tot_sed))
mantel_otu_sed <- data.frame(t(otu_table(phylo_tot_sed)))

mantel_taxa_wat <- data.frame(tax_table(phylo_tot_wat))
mantel_otu_wat <- data.frame(t(otu_table(phylo_tot_wat)))

TRCs_tab_sed <- cbind(mantel_otu_sed, TRCs_tab_sed)
TRCs_tab_wat <- cbind(mantel_otu_wat, TRCs_tab_wat)

mantel_otu_sed <- as.matrix(mantel_otu_sed)
mantel_taxa_sed <- as.matrix(mantel_taxa_sed)

mantel_phylo <- phyloseq(otu_table(mantel_otu_sed, taxa_are_rows = FALSE), #HZ samples
                         tax_table(mantel_taxa_sed), 
                         sample_data(TRCs_tab_sed))

mantel_otu_wat <- as.matrix(mantel_otu_wat)
mantel_taxa_wat <- as.matrix(mantel_taxa_wat)

mantel_phylo_wat <- phyloseq(otu_table(mantel_otu_wat, taxa_are_rows = FALSE), 
                             tax_table(mantel_taxa_wat), 
                             sample_data(TRCs_tab_wat))


phylo_dist <- phyloseq::distance(mantel_phylo, method = "bray")
phylo_dist_jac <- phyloseq::distance(mantel_phylo, method = "jaccard")
b1_dist <- dist(TRCs_tab_sed$B1, method = "euclidean")
cHET_dist <- dist(TRCs_tab_sed$cHET, method = "euclidean")
hmp_dist <- dist(TRCs_tab_sed$HMP, method = "euclidean")
AmMP_dist <- dist(TRCs_tab_sed$AmMP, method = "euclidean")
het_dist <- dist(TRCs_tab_sed$HET, method = "euclidean")
geo <- data.frame(TRCs_tab_sed$Long, TRCs_tab_sed$Lat)
geo <- as.matrix(geo)

lat_long_dist <- distm(geo, fun=distHaversine)
lat_long_dist <- as.dist(lat_long_dist)

mat <- data.frame(dist_b1 = as.vector(b1_dist), 
                  dist_cHET = as.vector(cHET_dist), 
                  dist_hmp = as.vector(hmp_dist), 
                  dist_AmMP = as.vector(AmMP_dist), 
                  dist_phylo = as.vector(phylo_dist), 
                  dist_het = as.vector(het_dist), 
                  dist_geo = as.vector(lat_long_dist))

phylo_b1 <- mantel(phylo_dist, b1_dist, permutations = 9999)
phylo_hmp <- mantel(phylo_dist, hmp_dist, permutations = 9999)
phylo_cHET <- mantel(phylo_dist, cHET_dist, permutations = 9999)
phylo_AmMP <- mantel(phylo_dist, AmMP_dist, permutations = 9999)
phylo_het <- mantel(phylo_dist, het_dist, permutations = 9999)
phylo_geo <- mantel(phylo_dist, lat_long_dist, permutations = 9999)

phylo_b1
phylo_hmp
phylo_cHET
phylo_het
phylo_geo

phylo_dist_wat <- phyloseq::distance(mantel_phylo_wat, method = "bray")
phylo_dist_jac_wat <- phyloseq::distance(mantel_phylo_wat, method = "jaccard")
b1_dist_wat <- dist(TRCs_tab_wat$B1, method = "euclidean")
cHET_dist_wat <- dist(TRCs_tab_wat$cHET, method = "euclidean")
hmp_dist_wat <- dist(TRCs_tab_wat$HMP, method = "euclidean")
AmMP_dist_wat <- dist(TRCs_tab_wat$AmMP, method = "euclidean")
het_dist_wat <- dist(TRCs_tab_wat$HET, method = "euclidean")
geo_wat <- data.frame(TRCs_tab_wat$Long, TRCs_tab_sed$Lat)
geo_wat <- as.matrix(geo_wat)

lat_long_dist_wat <- distm(geo_wat, fun=distHaversine)
lat_long_dist_wat <- as.dist(lat_long_dist_wat)

mat_wat <- data.frame(dist_b1 = as.vector(b1_dist_wat), 
                      dist_cHET = as.vector(cHET_dist_wat), 
                      dist_hmp = as.vector(hmp_dist_wat), 
                      dist_AmMP = as.vector(AmMP_dist_wat), 
                      dist_phylo = as.vector(phylo_dist_wat), 
                      dist_het = as.vector(het_dist_wat), 
                      dist_geo = as.vector(lat_long_dist_wat))

phylo_b1_wat <- mantel(phylo_dist_wat, b1_dist_wat, permutations = 9999)
phylo_hmp_wat <- mantel(phylo_dist_wat, hmp_dist_wat, permutations = 9999)
phylo_cHET_wat <- mantel(phylo_dist_wat, cHET_dist_wat, permutations = 9999)
phylo_AmMP_wat <- mantel(phylo_dist_wat, AmMP_dist_wat, permutations = 9999)
phylo_het_wat <- mantel(phylo_dist_wat, het_dist_wat, permutations = 9999)
phylo_geo_wat <- mantel(phylo_dist_wat, lat_long_dist_wat, permutations = 9999)

phylo_b1_wat
phylo_hmp_wat
phylo_cHET_wat
phylo_het_wat
phylo_geo_wat

###Section 2: Differential abundance

#turn "Kingdom" tab into the ASV list because ancbombc only recognizes: Kingdom, Phylum, Class, Order, Family, Genus

setwd("~/OneDrive/RStudio")
all_taxa <- readRDS("SA_all_taxa")
phylo_tot <- all_taxa
asvs <- rownames(otu_table(phylo_tot))
#or
asvs <- taxa_names(phylo_tot) #should give the same results as above as long as rownames in the tax and count tables are ASVs

phylo_tax <- data.frame(tax_table(phylo_tot))
phylo_tax$"Kingdom" <- asvs
phylo_tax <- as.matrix(phylo_tax)
phylo_otu <- data.frame(otu_table(phylo_tot))
phylo_otu <- as.matrix(phylo_otu)
meta <- data.frame(sample_data(phylo_tot))

phylo_tot <- phyloseq(otu_table(phylo_otu, taxa_are_rows = TRUE), 
                      tax_table(phylo_tax), sample_data(meta))

#handy hack to get around ancombc only recognizing the taxonomic levels returned by: <rank_names(phyloseq_object)>

#put in raw data to ancombc, by sample type
phylo_tot_sed <- subset_samples(phylo_tot, Type == "Sediment")
phylo_tot_wat <- subset_samples(phylo_tot, Type == "Water")

set.seed(10023)
diff_abund_b1_sed <- ancombc(phylo_tot_sed, 
                             formula = "B1+HMP+cHET+AmMP+HET", 
                             group = NULL,
                             tax_level = "Kingdom")
diff_abund_b1_wat <- ancombc(phylo_tot_wat, 
                             formula = "B1+HMP+cHET+AmMP+HET", 
                             group = NULL,
                             tax_level = "Kingdom")

saveRDS(diff_abund_b1_sed, file = "ancombc_trc_sed_raw")
saveRDS(diff_abund_b1_wat, file = "ancombc_trc_wat_raw")


diff_abund_b1_sed <- readRDS("ancombc_trc_sed_raw")
diff_abund_b1_wat <- readRDS("ancombc_trc_wat_raw")

#<tax_level = "Kingdom"> is really ASV-by-ASV, not actual Kingdom level. 

res_sed <- diff_abund_b1_sed$res #results
res_wat <- diff_abund_b1_wat$res #results

#extract log-fold changes and standard errors
df_lfc_sed = data.frame(res_sed$lfc * res_sed$diff_abn, check.names = FALSE) %>% 
  rownames_to_column("taxon")

df_se_sed = data.frame(res_sed$se * res_sed$diff_abn, check.names = FALSE) %>% 
  rownames_to_column("taxon")
colnames(df_se_sed)[-1] = paste0(colnames(df_se_sed)[-1], "SE")

df_lfc_wat = data.frame(res_wat$lfc * res_wat$diff_abn, check.names = FALSE) %>% 
  rownames_to_column("taxon")

df_se_wat = data.frame(res_wat$se * res_wat$diff_abn, check.names = FALSE) %>% 
  rownames_to_column("taxon")
colnames(df_se_wat)[-1] = paste0(colnames(df_se_wat)[-1], "SE")
df_se_wat <- data.frame(df_se_wat)
#

#extract results from every individual independent variable, from sediment and water results separately 
lfc_sed_b1 <- df_lfc_sed %>% 
  dplyr::select(B1, taxon)
lfc_wat_b1 <- df_lfc_wat %>% 
  dplyr::select(B1, taxon)

se_sed_b1 <- df_se_sed %>% 
  dplyr::select(B1SE, taxon)
se_wat_b1 <- df_se_wat %>% 
  dplyr::select(B1SE, taxon)


lfc_sed_hmp <- df_lfc_sed %>% 
  dplyr::select(HMP, taxon)
lfc_wat_hmp <- df_lfc_wat %>% 
  dplyr::select(HMP, taxon)

se_sed_hmp <- df_se_sed %>% 
  dplyr::select(HMPSE, taxon)
se_wat_hmp <- df_se_wat %>% 
  dplyr::select(HMPSE, taxon)


lfc_sed_chet <- df_lfc_sed %>% 
  dplyr::select(cHET, taxon)
lfc_wat_chet <- df_lfc_wat %>% 
  dplyr::select(cHET, taxon)

se_sed_chet <- df_se_sed %>% 
  dplyr::select(cHETSE, taxon)
se_wat_chet <- df_se_wat %>% 
  dplyr::select(cHETSE, taxon)



lfc_sed_het <- df_lfc_sed %>% 
  dplyr::select(HET, taxon)
lfc_wat_het <- df_lfc_wat %>% 
  dplyr::select(HET, taxon)

se_sed_het <- df_se_sed %>% 
  dplyr::select(HETSE, taxon)
se_wat_het <- df_se_wat %>% 
  dplyr::select(HETSE, taxon)



lfc_sed_ammp <- df_lfc_sed %>% 
  dplyr::select(AmMP, taxon)
lfc_wat_ammp <- df_lfc_wat %>% 
  dplyr::select(AmMP, taxon)

se_sed_ammp <- df_se_sed %>% 
  dplyr::select(AmMPSE, taxon)
se_wat_ammp <- df_se_wat %>% 
  dplyr::select(AmMPSE, taxon)
#

#make all columns consistent between results so they can eventually all be combined with rbind()
final_b1_sed <- lfc_sed_b1 %>% 
  dplyr::left_join(se_sed_b1, by = "taxon") %>% 
  dplyr::transmute(taxon, B1, B1SE) %>%
  dplyr::filter(B1 != 0) %>% 
  dplyr::arrange(desc(B1)) %>%
  dplyr::mutate(direct = ifelse(B1 > 0, "Positive", "Negative"))

final_b1_wat <- lfc_wat_b1 %>% 
  dplyr::left_join(se_wat_b1, by = "taxon") %>% 
  dplyr::transmute(taxon, B1, B1SE) %>%
  dplyr::filter(B1 != 0) %>% 
  dplyr::arrange(desc(B1)) %>%
  dplyr::mutate(direct = ifelse(B1 > 0, "Positive", "Negative"))


final_hmp_sed <- lfc_sed_hmp %>% 
  dplyr::left_join(se_sed_hmp, by = "taxon") %>% 
  dplyr::transmute(taxon, HMP, HMPSE) %>%
  dplyr::filter(HMP != 0) %>% 
  dplyr::arrange(desc(HMP)) %>%
  dplyr::mutate(direct = ifelse(HMP > 0, "Positive", "Negative"))

final_hmp_wat <- lfc_wat_hmp %>% 
  dplyr::left_join(se_wat_hmp, by = "taxon") %>% 
  dplyr::transmute(taxon, HMP, HMPSE) %>%
  dplyr::filter(HMP != 0) %>% 
  dplyr::arrange(desc(HMP)) %>%
  dplyr::mutate(direct = ifelse(HMP > 0, "Positive", "Negative"))


final_chet_sed <- lfc_sed_chet %>% 
  dplyr::left_join(se_sed_chet, by = "taxon") %>% 
  dplyr::transmute(taxon, cHET, cHETSE) %>%
  dplyr::filter(cHET != 0) %>% 
  dplyr::arrange(desc(cHET)) %>%
  dplyr::mutate(direct = ifelse(cHET > 0, "Positive", "Negative"))

final_chet_wat <- lfc_wat_chet %>% 
  dplyr::left_join(se_wat_chet, by = "taxon") %>% 
  dplyr::transmute(taxon, cHET, cHETSE) %>%
  dplyr::filter(cHET != 0) %>% 
  dplyr::arrange(desc(cHET)) %>%
  dplyr::mutate(direct = ifelse(cHET > 0, "Positive", "Negative"))


final_ammp_sed <- lfc_sed_ammp %>% 
  dplyr::left_join(se_sed_ammp, by = "taxon") %>% 
  dplyr::transmute(taxon, AmMP, AmMPSE) %>%
  dplyr::filter(AmMP != 0) %>% 
  dplyr::arrange(desc(AmMP)) %>%
  dplyr::mutate(direct = ifelse(AmMP > 0, "Positive", "Negative"))

final_ammp_wat <- lfc_wat_ammp %>% 
  dplyr::left_join(se_wat_ammp, by = "taxon") %>% 
  dplyr::transmute(taxon, AmMP, AmMPSE) %>%
  dplyr::filter(AmMP != 0) %>% 
  dplyr::arrange(desc(AmMP)) %>%
  dplyr::mutate(direct = ifelse(AmMP > 0, "Positive", "Negative"))


final_het_sed <- lfc_sed_het %>% 
  dplyr::left_join(se_sed_het, by = "taxon") %>% 
  dplyr::transmute(taxon, HET, HETSE) %>%
  dplyr::filter(HET != 0) %>% 
  dplyr::arrange(desc(HET)) %>%
  dplyr::mutate(direct = ifelse(HET > 0, "Positive", "Negative"))

final_het_wat <- lfc_wat_het %>% 
  dplyr::left_join(se_wat_het, by = "taxon") %>% 
  dplyr::transmute(taxon, HET, HETSE) %>%
  dplyr::filter(HET != 0) %>% 
  dplyr::arrange(desc(HET)) %>%
  dplyr::mutate(direct = ifelse(HET > 0, "Positive", "Negative")) 
#

#change independent variable to the same column name for all, to be combined with rbind()
final_b1_vit <- final_b1_sed %>% 
  dplyr::rename(LFC = B1, 
                SE = B1SE) %>% 
  mutate(vitamer = "B1")

final_hmp_vit <- final_hmp_sed %>% 
  dplyr::rename(LFC = HMP, 
                SE = HMPSE) %>% 
  mutate(vitamer = "HMP")
final_hmp_vit <- final_hmp_vit[c(1:10,119:128),] #top 10 most significant ASVs in either direction

final_chet_vit <- final_chet_sed %>% 
  dplyr::rename(LFC = cHET, 
                SE = cHETSE) %>% 
  mutate(vitamer = "cHET")
final_chet_vit <- final_chet_vit[c(1:10, 35:44),] #top 10 most significant ASVs in either direction

final_ammp_vit <- final_ammp_sed %>% 
  dplyr::rename(LFC = AmMP, 
                SE = AmMPSE) %>% 
  mutate(vitamer = "AmMP")
final_ammp_vit <- final_ammp_vit[c(1:8,14:23),] #top 10 most significant ASVs in either direction

final_het_vit <- final_het_sed %>% 
  dplyr::rename(LFC = HET, 
                SE = HETSE) %>% 
  mutate(vitamer = "HET")
final_het_vit <- final_het_vit[c(1:10,129:138),] #top 10 most significant ASVs in either direction

#rename columns so they all match up, again to be combined with rbind()
final_b1_vit_wat <- final_b1_sed %>% 
  dplyr::rename(LFC = B1, 
                SE = B1SE) %>% 
  mutate(vitamer = "B1")

final_hmp_vit_wat <- final_hmp_wat %>% 
  dplyr::rename(LFC = HMP, 
                SE = HMPSE) %>% 
  mutate(vitamer = "HMP")

final_chet_vit_wat <- final_chet_wat %>% 
  dplyr::rename(LFC = cHET, 
                SE = cHETSE) %>% 
  mutate(vitamer = "cHET")

final_ammp_vit_wat <- final_ammp_wat %>% 
  dplyr::rename(LFC = AmMP, 
                SE = AmMPSE) %>% 
  mutate(vitamer = "AmMP")

final_het_vit_wat <- final_het_wat %>% 
  dplyr::rename(LFC = HET, 
                SE = HETSE) %>% 
  mutate(vitamer = "HET")
#

###the data has now been filtered to only include the most influential ASVs, per independent variable
###all column names correspond so the data can all be combined to plot 

#merge water data
vitamers_wat <- rbind(final_b1_vit_wat, final_hmp_vit_wat)
vitamers_wat <- rbind(vitamers_wat, final_chet_vit_wat)
vitamers_wat <- rbind(vitamers_wat, final_het_vit_wat)
vitamers_wat <- rbind(vitamers_wat, final_ammp_vit_wat) %>% 
  left_join(tot_tax, "taxon") %>% 
  mutate("Sample_Type" = "Water")
saveRDS(vitamers_wat, "vitamers_wat")
vitamers_wat <- readRDS("vitamers_wat")
#

#merge sediment data
vitamers_sed <- rbind(final_b1_vit, final_hmp_vit)
vitamers_sed <- rbind(vitamers_sed, final_chet_vit)
vitamers_sed <- rbind(vitamers_sed, final_het_vit)
vitamers_sed <- rbind(vitamers_sed, final_ammp_vit) %>% 
  left_join(tot_tax, "taxon") %>% 
  mutate("Sample_Type" = "Sediment")
saveRDS(vitamers_sed, "vitamers_sed")
vitamers_sed <- readRDS("vitamers_sed")
#

#get taxonomic and counts info to all be combined 
tot_tax <- data.frame(tax_table(phylo_tot))
tot_tax <- tot_tax %>% 
  dplyr::rename("taxon" = Kingdom) #change to taxon so left_join() can be used to merge differentially abundant taxa with total taxa, based on corresponding ASVs
tot_otu <- data.frame(otu_table(phylo_tot)) 
tot_otu <- tot_otu %>% 
  dplyr::mutate("taxon" = rownames(tot_otu)) #same deal as above ^
#

#final data to feed into ggplot
vitamers_sed$"Sample_Type" <- "HZ"
vitamers_wat$"Sample_Type" <- "SW"
vitamers <- rbind(vitamers_sed, vitamers_wat)

vitamers <- vitamers %>% 
  mutate(across(.cols = "Family", .fns = str_replace, "uncultured", "Vicinamibacterales_uncultured")) %>% 
  mutate(across("Family", str_replace, "uncultured_fa", "Verrucomicrobiae_uncultured"))
vitamers <- vitamers %>% 
  mutate(across("Family", str_replace, "Vicinamibacterales_Verrucomicrobiae_uncultured", "Verrucomicrobiae_uncultured"))

vitamers_pos <- vitamers %>% 
  filter(direct == "Positive")
vitamers_neg <- vitamers %>% 
  filter(direct == "Negative")
#

###this data should now have full taxonomic information, counts information per sample, and sample type information to be plotted to your heart's content
###it should be ONLY differentially abundant ASVs (those matching between raw counts/taxonomy tables and those spit out by ancombc)

#(+) log-fold change plot

families_pos <- vitamers_pos$Family
alph_pos <- rev(sort(families_pos))

(pos_plot <- ggplot(vitamers_pos, aes(x = vitamer, y = Family, color = Sample_Type)) + 
    geom_point(aes(size = LFC)) +
    scale_size_continuous(limits = c(0.0001, 7), breaks = c(.0001, .01, 1, 3, 7), 
                          labels = c("0.0001", "0.01", "1.0", "3.0", "7.0"), 
                          name = "LFC") +
    theme_bw() +
    scale_x_discrete(limits = c("B1", "HMP", "cHET", "HET", "AmMP")) +
    facet_wrap(~Sample_Type, scales = "free_y", nrow = 3) +
    theme(axis.text.x = element_text(angle = 45, hjust=1, face = "bold", color = "black", size = 12), 
          strip.text = element_blank(), 
          axis.text.y = element_text(color = "black", size = 12), 
          axis.title.y = element_text(size = 14, face = "bold", color = "black")) +
    labs(x = "") +
    scale_color_manual(values = c("#d8b365","#5ab4ac")))



ggsave("posplot.png", pos_plot, width = 5.5, height = 10, device = "png")

#(-) log-fold change plot
(neg_plot <- ggplot(vitamers_neg, aes(x = vitamer, y = Family, color = Sample_Type)) + 
    geom_point(aes(size = -1*LFC)) + #I made this positive then manually labeled in <scale_size_continuous()> as negative because ggplot doesn't like negative values
    theme_bw() +  
    scale_size_continuous(limits = c(0.001, 6), breaks = c(.001, .1, 1, 3, 6), 
                          labels = c("-0.001", "-0.1", "-1.0", "-3.0", "-6.0"), name = "LFC") +
    scale_x_discrete(limits = c("B1", "HMP", "cHET", "HET", "AmMP")) +
    facet_wrap(~Sample_Type, scales = "free_y", nrow = 3) +
    theme(axis.text.x = element_text(angle = 45, hjust=1, face = "bold", color = "black", size = 12), 
          strip.text = element_blank(), 
          axis.text.y = element_text(, color = "black", size = 12), 
          axis.title.y = element_text(size = 14, color = "black")) +
    labs(x = "") +
    scale_color_manual(values = c("#d8b365","#5ab4ac"), label = c("HZ", "SW"), 
                       name = "Type"))

neg_plot

ggsave("negplot.png", neg_plot, width = 5.5, height = 10, device = "png")

#I just brought these .png files into powerpoint (or illulstrator) and put the final touches on them there. 

#Richness values for section 1: 

#positive
vitamers_sed_otus_pos_b1 <- vitamers_sed %>% 
  left_join(tot_otu, "taxon") %>% 
  filter(direct == "Positive", 
         vitamer == "B1") %>% 
  dplyr::select(BCCIGS:SSCSW)

vitamers_sed_otus_pos_hmp <- vitamers_sed %>% 
  left_join(tot_otu, "taxon") %>% 
  filter(direct == "Positive", 
         vitamer == "HMP") %>% 
  dplyr::select(BCCIGS:SSCSW)

vitamers_sed_otus_pos_chet <- vitamers_sed %>% 
  left_join(tot_otu, "taxon") %>% 
  filter(direct == "Positive", 
         vitamer == "cHET") %>% 
  dplyr::select(BCCIGS:SSCSW)

vitamers_sed_otus_pos_het <- vitamers_sed %>% 
  left_join(tot_otu, "taxon") %>% 
  filter(direct == "Positive", 
         vitamer == "HET") %>% 
  dplyr::select(BCCIGS:SSCSW)

vitamers_sed_otus_pos_ammp <- vitamers_sed %>% 
  left_join(tot_otu, "taxon") %>% 
  filter(direct == "Positive", 
         vitamer == "AmMP") %>% 
  dplyr::select(BCCIGS:SSCSW)

diffabund_pos_sums_sed_b1 <- colSums(vitamers_sed_otus_pos_b1)
diffabund_pos_sums_sed_hmp <- colSums(vitamers_sed_otus_pos_hmp)
diffabund_pos_sums_sed_chet <- colSums(vitamers_sed_otus_pos_chet)
diffabund_pos_sums_sed_het <- colSums(vitamers_sed_otus_pos_het)
diffabund_pos_sums_sed_ammp <- colSums(vitamers_sed_otus_pos_ammp)

#negative
vitamers_wat_otus_neg_b1 <- vitamers_wat %>% 
  left_join(tot_otu, "taxon") %>% 
  filter(direct == "Negative", 
         vitamer == "B1") %>% 
  dplyr::select(BCCIGS:SSCSW)

vitamers_wat_otus_neg_hmp <- vitamers_wat %>% 
  left_join(tot_otu, "taxon") %>% 
  filter(direct == "Negative", 
         vitamer == "HMP") %>% 
  dplyr::select(BCCIGS:SSCSW)

vitamers_wat_otus_neg_chet <- vitamers_wat %>% 
  left_join(tot_otu, "taxon") %>% 
  filter(direct == "Negative", 
         vitamer == "cHET") %>% 
  dplyr::select(BCCIGS:SSCSW)

vitamers_wat_otus_neg_het <- vitamers_wat %>% 
  left_join(tot_otu, "taxon") %>% 
  filter(direct == "Negative", 
         vitamer == "HET") %>% 
  dplyr::select(BCCIGS:SSCSW)

vitamers_wat_otus_neg_ammp <- vitamers_wat %>% 
  left_join(tot_otu, "taxon") %>% 
  filter(direct == "Negative", 
         vitamer == "AmMP") %>% 
  dplyr::select(BCCIGS:SSCSW)

diffabund_neg_sums_wat_b1 <- colSums(vitamers_wat_otus_neg_b1)
diffabund_neg_sums_wat_hmp <- colSums(vitamers_wat_otus_neg_hmp)
diffabund_neg_sums_wat_chet <- colSums(vitamers_wat_otus_neg_chet)
diffabund_neg_sums_wat_het <- colSums(vitamers_wat_otus_neg_het)
diffabund_neg_sums_wat_ammp <- colSums(vitamers_wat_otus_neg_ammp)

###Section 3: Alpha diversity plots

all_taxa_rare <- phylo_SA
alpha_SA <- data.frame(MicrobiotaProcess::get_alphaindex(all_taxa_rare))
alpha_SA_sed <- alpha_SA %>% filter(Type == "Sediment")
alpha_SA_wat <- alpha_SA %>% filter(Type == "Water")
alpha_SA$Location <- gsub("FRU", "FR", alpha_SA$Location)
alpha_SA$Location <- gsub("FRD", "FR", alpha_SA$Location)
alpha_SA$Location <- gsub("FR2", "FR", alpha_SA$Location)
alpha_SA$Location <- gsub("FR3", "FR", alpha_SA$Location)
alpha_SA <- alpha_SA %>% 
  mutate(across('Region', str_replace, 'Undammed', 'Battle&Butte'))

(alph <- ggplot(data = alpha_SA, mapping = aes(x = Location, y = Shannon)) +
    geom_boxplot(mapping = aes(color = Type)) +
    theme_bw() +
    labs(y = "Shannon Diversity Index", x = "", 
         title = "", color = "Type") +
    scale_color_manual(labels = c("Hyporheic Zone","Surface Water"), values = c("#d8b365","#5ab4ac")) +
    scale_x_discrete(limits = c("BCC", "NFB", "SFB", "SAC", "SSC", "CCM", "CCH", "FR"))) +
  theme(axis.text.x = element_text(face = "bold", size = 14, color = "black"), 
        axis.text.y = element_text(face = "bold", size = 14, colour = "black"), 
        axis.title = element_text(size = 16, colour = "black"))

(alpha <- alph <- ggplot(data = alpha_SA, mapping = aes(x = Region, y = Shannon)) +
    geom_boxplot(mapping = aes(color = Type)) +
    theme_bw() +
    labs(y = "Shannon Diversity Index", x = "", 
         title = "", color = "Type") +
    scale_color_manual(labels = c("Hyporheic Zone","Surface Water"), values = c("#d8b365","#5ab4ac")) +
    theme(axis.text.x = element_text(face = "bold", size = 14, color = "black", angle = 45, vjust = .6), 
          axis.text.y = element_text(face = "bold", size = 14, colour = "black"), 
          axis.title = element_text(size = 16, colour = "black")))

alph.time <- ggplot(data = alpha_SA, 
                    mapping = aes(x = Timepoint, y = Shannon, 
                                  color = Timepoint)) +
  geom_boxplot(position = "dodge") +
  guides(color = "none")


(alpha_time <- alph.time +
    labs(y = "Shannon Diversity Index", x = "", 
         title = "", shape = "Region") +
    scale_color_manual(values = c("#7570b3", "#1b9e77", "#e7298a", "#d95f02")) +
    geom_point(position=position_dodge(width = .5),
               aes(shape = Region), size = 3.2, fill = "darkslategray", color = "black") +
    scale_x_discrete(limits = c("PS", "S", "IG", "P"), 
                     labels = c("Pre-Spawn", "Spawn", "Incubation", "Hatched")) +
    theme_bw() +
    scale_shape_manual(values = c(21,22,23,24), limits = c("Battle&Butte", 
                                                           "Clear Creek", 
                                                           "Feather River", 
                                                           "Sacramento")) +
    theme(axis.text.x = element_text(face = "bold", size = 14, color = "black", angle = 45, vjust = .6), 
          axis.text.y = element_text(face = "bold", size = 14, colour = "black"), 
          axis.title = element_text(size = 16, colour = "black")))

ggsave(alpha_time, height = 5, width = 7, 
       filename = "/Users/kellyshannon/OneDrive/RStudio/SAC_plots/shannondiv_time_.png")
ggsave(alpha, height = 5, width = 7, 
       filename = "/Users/kellyshannon/OneDrive/RStudio/SAC_plots/shannondiv_region_.png")


#Relative Abundance Bubble Plot

#Filter data to get top n Classes

##by region
all_taxa_rare <- phylo_SA
all_taxa_rare_region <- merge_samples(all_taxa_rare, group = 'Region', fun = sum)
SA_fam_reg <- microbiome::aggregate_taxa(all_taxa_rare_region, level = "Family")
SA_fam_reg <- transform_sample_counts(SA_fam_reg, function(OTU) OTU/sum(OTU))
topn_fam_reg <- names(sort(taxa_sums(SA_fam_reg), decreasing = TRUE))[1:30]
SA_topn_fam_reg <- prune_taxa(topn_fam_reg, SA_fam_reg)

top_fam_reg <- data.frame(otu_table(SA_topn_fam_reg))

top_fam_tib_reg <- as_tibble(top_fam_reg) %>% 
  mutate(across(.cols = 'Clear.Creek':Undammed, .fns = decimal_to_percent)) %>% 
  add_column(Family = rownames(top_fam_reg), .before = "Clear.Creek") %>% 
  mutate(across('Family', str_replace, 'Unknown', 'Bacillariophyta_unclassified')) %>% 
  mutate(across('Family', str_replace, '\tBacteria_Verrucomicrobiota_Verrucomicrobiae_uncultured_uncultured_fa', 
                'Verrucomicrobiae_uncultured')) %>% 
  mutate(across('Family', str_replace, 'Clade_III', 'SAR11_Clade_III'))
##

top_fam_tib_reg_seas <- top_fam_tib_reg %>% right_join(top_fam_tib_seas, by = 'Family')

top_fam_tib_reg_seas[30,2] <- 0.37

top_fam_tib_reg_seas[30,3] <- 0.19

top_fam_tib_reg_seas[30,4] <- 0.14

top_fam_tib_reg_seas[30,5] <- 2.39

top_fam_tib_reg <- top_fam_tib_reg_seas %>% 
  select(Family, Clear.Creek, Feather.River, Sacramento, Undammed)

top_fam_tib_seas <- top_fam_tib_reg_seas %>% 
  select(Family, PS, S, IG, P)

families <- top_fam_tib_reg_seas[['Family']]
families_seas <- top_fam_tib_seas[['Family']]

decimal_to_percent <- function(num){
  num*100
}


alph <- sort(families)


top_fam_bub <- top_fam_tib_reg %>% 
  pivot_longer(cols = c('Clear.Creek', 'Undammed', 'Sacramento', 'Feather.River'), names_to = 'Region') %>% 
  mutate(across('Region', str_replace, 'Undammed', 'Battle&Butte')) %>% 
  mutate(across('Region', str_replace, 'Feather.River', 'Feather River')) %>% 
  mutate(across('Region', str_replace, 'Clear.Creek', 'Clear Creek')) %>% 
  mutate(across('value', replace_na, 0))

top_fam_bub_seas <- top_fam_tib_seas %>% 
  pivot_longer(cols = c('PS', 'S', 'IG', 'P'), names_to = 'Season') %>% 
  mutate(across('value', replace_na, 0))

top_fam_bub_reg_seas <- top_fam_bub %>% rename(value_reg = value) %>% 
  left_join(y = top_fam_bub_seas, by = 'Family', multiple = "all")

(region_bub <- top_fam_bub_reg_seas %>% ggplot(aes(x = Region, y = Family)) + 
    geom_point(aes(size = value_reg)) +
    scale_size_continuous(breaks = c(0.1, 1, 10, 20, 30), range = c(0,10),
                          labels = c("0.1%", "1%", "10%", "20%", "30%"), limits = c(.1,100)) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 60, size = 12, color = "black", vjust = .9, hjust = .95, face = "bold"), 
          axis.text.y = element_text(size = 10, color = "black")) +
    ylab("") + xlab("") + labs(size = "Relative Abundance"))

ggsave(region_bub, height = 7, width = 5.5, 
       filename = "/Users/kellyshannon/OneDrive/RStudio/SAC_plots/regions_bub.png")

##by season
all_taxa_rare <- phylo_SA
all_taxa_rare_season <- merge_samples(all_taxa_rare, group = 'Timepoint', fun = sum)
SA_fam_seas <- microbiome::aggregate_taxa(all_taxa_rare_season, level = "Family")
SA_fam_seas <- transform_sample_counts(SA_fam_seas, function(OTU) OTU/sum(OTU))
topn_fam_seas <- names(sort(taxa_sums(SA_fam_seas), decreasing = TRUE))[1:30]
SA_topn_fam_seas <- prune_taxa(topn_fam_seas, SA_fam_seas)

top_fam_seas <- data.frame(otu_table(SA_topn_fam_seas))

top_fam_tib_seas <- as_tibble(top_fam_seas) %>% 
  mutate(across(.cols = IG:S, .fns = decimal_to_percent)) %>% 
  add_column(Family = rownames(top_fam_seas), .before = "IG") %>% 
  mutate(across('Family', str_replace, 'Unknown', 'Bacillariophyta_unclassified')) %>% 
  mutate(across('Family', str_replace, '\tBacteria_Verrucomicrobiota_Verrucomicrobiae_uncultured_uncultured_fa', 
                'Verrucomicrobiae_uncultured')) %>% 
  mutate(across('Family', str_replace, 'Clade_III', 'SAR11_Clade_III'))

top_fam_bub_seas <- top_fam_tib_seas %>% 
  pivot_longer(cols = c('IG', 'P', 'PS', 'S'), names_to = 'Season') %>% 
  mutate(value = as.numeric(value))

(seasons_bub <- top_fam_bub_seas %>% ggplot(aes(x = Season, y = Family)) + 
    geom_point(aes(size = value)) +
    scale_size_continuous(breaks = c(0.1, 1, 10, 20, 30), range = c(0,10),
                          labels = c("0.1%", "1%", "10%", "20%", "30%"), limits = c(.1,100)) + 
    theme_bw() + scale_y_discrete(limits = alph) + 
    scale_x_discrete(limits = c("PS", "S", "IG", "P"), labels = c("Pre-Spawn", "Spawn", "Incubation", "Hatched")) +
    theme(axis.text.x = element_text(angle = 60, size = 12, color = "black", vjust = .9, hjust = .95, face = "bold"), 
          axis.text.y = element_text(size = 10, color = "black")) +
    ylab("") + xlab("") + labs(size = "Relative Abundance"))

ggsave(seasons_bub, height = 7, width = 5.5, 
       filename = "/Users/kellyshannon/OneDrive/RStudio/SAC_plots/seasons_bub.png")

###Section 4: Beta diversity plots

#CAP

phylo_SA.trans <- phyloseq::subset_taxa(phylo_SA, Kingdom == "\tBacteria" | Kingdom == "\tArchaea")

cca.ord_meta <- ordinate(phylo_SA.trans, method = "CAP", 
                         distance = "bray", formula = ~Type+Timepoint+Region)

cap.plot <- plot_ordination(phylo_SA, cca.ord_meta, shape = "Type", color = "Region") + 
  theme_bw() +
  geom_point(size=2.5)

cap.plot + scale_color_manual(values = c("#e41a1c", "#377eb8", "darkgoldenrod1", "#4daf4a"), 
                              labels = c("Clear Creek", "Feather River", 
                                         "Hatched", "Battle&Butte")) +
  scale_shape_discrete(labels = c("HZ","SW")) +
  stat_ellipse(aes(group = Region)) +
  theme_classic()

#NMDS plots
#Prep Sediment and Water Tables 
SA_phylo.sed <- subset_samples(phylo.SA_rarefied, Sample_Type == "S")
SA_phylo.wat <- subset_samples(phylo.SA_rarefied, Sample_Type == "W")

SA_ord.sed <- ordinate(SA_phylo.sed, method = "NMDS", distance = "bray")
SA_ord.wat <- ordinate(SA_phylo.wat, method = "NMDS", distance = "bray")

meta.sed <- meta_SA %>% 
  filter(Sample_Type == "S")

meta.sed_num <- meta.sed %>% 
  dplyr::select(EC, SPC, TDS, Sal, DO_per, DO_conc, pH, Turb, Chl, AmMP, cHET, HMP, B1)
meta.sed_num[] <- lapply(meta.sed_num, function(x) as.numeric(as.character(x)))
meta.sed <- cbind(meta.sed_num,meta.sed[,c(20,21,22,24)])

meta.wat <- meta_SA %>% 
  filter(Sample_Type == "W")
meta.wat_num <- meta.wat %>% 
  dplyr::select(Temp, EC, SPC, TDS, Sal, DO_per, DO_conc, pH, Turb, Chl)
meta.wat_num[] <- lapply(meta.wat_num, function(x) as.numeric(as.character(x)))
meta.wat <- cbind(meta.wat_num,meta.wat[,c(20,21,22,24)])

enfit.sed <- envfit(phylo.sed_ord, meta.sed, perm = 9999, na.rm = TRUE)
enfit.sed

phylo.wat_ord <- ordinate(SA_phylo.wat, method = "NMDS", distance = "bray")
enfit.wat <- envfit(phylo.wat_ord, meta.wat, perm = 9999, na.rm = TRUE)
enfit.wat


#Sediment Samples
SA_data.scores.sed <- scores(SA_ord.sed)
SA_data.scores.sed <- as.data.frame(SA_data.scores.sed$sites)
SA_data.scores.sed$"Season" <- meta.sed$Timepoint
SA_data.scores.sed$"Type" <- meta.sed$Sample_Type
SA_data.scores.sed$"Station" <- meta.sed$Location

set.seed(2001)
SA_enfit.sed <- envfit(SA_ord.sed, meta.sed_num, perm = 9999, na.rm = TRUE)
SA_enfit.sed

A <- as.list(SA_enfit.sed$vectors)
pvals <- as.data.frame(A$pvals)
arrows<-as.data.frame(A$arrows*sqrt(A$r))
C<-cbind(arrows, pvals)
Cred<-subset(C,pvals<0.05)
Cred$"Variable" <- rownames(Cred)

sed_ord <- ggplot(data = SA_data.scores.sed, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = SA_data.scores.sed, aes(shape = Season, color = Station), size = 4) +
  theme_classic() +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = Cred, size =1, colour = "darkgray", 
               arrow = arrow(length = unit(0.2,"cm"))) +
  geom_text_repel(data = Cred, aes(x = NMDS1, y = NMDS2), 
                  colour = "black", size = 4, label = Cred$Variable, 
                  direction = "both") +
  ggtitle("Sacramento River Sediment Community Differences with Metadata Vectors") +
  theme(title = element_text(size = 12))


SA_data.scores.wat <- scores(SA_ord.wat)
SA_data.scores.wat <- as.data.frame(SA_data.scores.wat$sites)
SA_data.scores.wat$"Season" <- meta.wat$Timepoint
SA_data.scores.wat$"Type" <- meta.wat$Sample_Type
SA_data.scores.wat$"Station" <- meta.wat$Location

set.seed(2001)
SA_enfit.wat <- envfit(SA_ord.wat, meta.sed_num, perm = 9999, na.rm = TRUE)
SA_enfit.wat

A <- as.list(SA_enfit.wat$vectors)
pvals <- as.data.frame(A$pvals)
arrows<-as.data.frame(A$arrows*sqrt(A$r))
C<-cbind(arrows, pvals)
Cred<-subset(C,pvals<0.05)
Cred$"Variable" <- rownames(Cred)

water_ord <- ggplot(data = SA_data.scores.wat, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = SA_data.scores.wat, aes(shape = Season, color = Station), size = 4) +
  theme_classic() +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = Cred, size =1, colour = "darkgray", 
               arrow = arrow(length = unit(0.2,"cm"))) +
  geom_text_repel(data = Cred, aes(x = NMDS1, y = NMDS2), 
                  colour = "black", size = 4, label = Cred$Variable, 
                  direction = "both") +
  ggtitle("Sacramento River Sediment Community Differences with Metadata Vectors") +
  theme(title = element_text(size = 12))



















