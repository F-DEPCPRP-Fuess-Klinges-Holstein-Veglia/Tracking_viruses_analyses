library("phyloseq")
library("ggplot2")
library("tidyverse")
library("dplyr")
library("tibble") 
library("speedyseq") 
library("microbiome")
library("ggpubr")
library("vegan")
library("cowplot")
library(DESeq2)
library(scales)
library(grid)
library(reshape2)
library(IHW)
library(phylosmith)
library(ggplotify)

############ Preparing phyloseq objects
#prepare taxonomy and sample info sheets
tax_mat <- read.csv("FL_DT_rnaseq_virus_db_all_seqs_3_15_25_vSUM_classifiedt2_classified_vSUM.csv")
tax_mat$length <- as.numeric(gsub("[^0-9]", "", tax_mat$length))
sample_df <- read.csv("FL_sample_sheet_june25_update.csv")

#define row names for phyloseq
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("seq")
tax_mat <- as.matrix(tax_mat)

sample_df <- sample_df %>% 
  tibble::column_to_rownames("sample")

TAX = tax_table(tax_mat)
samples = sample_data(sample_df)
#prepare TPM normalized counts for phyloseq
##new code
counts <- read.csv("FL_RNA_virus_classseq_counts5_25.isoform.TPM.not_cross_norm.csv", row.names = 1, header = TRUE, check.names=FALSE, sep = ",")
counts_matrix <- as.matrix(counts)
OTU1 = otu_table(counts_matrix, taxa_are_rows = TRUE)

fldnaVirTPM <- phyloseq(OTU1, TAX, samples)

saveRDS(fldnaVirTPM, file="flRNAV_june25_TPM_phylo.rds")
flrnaVirTPM <- readRDS("flRNAV_june25_TPM_phylo.rds")

#prepare TMM normalized counts for phyloseq
##new code
counts1 <- read.csv("FL_RNA_virus_classseq_counts5_25.TMM.EXPR.csv", row.names = 1, header = TRUE, check.names=FALSE, sep = ",")
counts_matrix1 <- as.matrix(counts1)
OTU2 = otu_table(counts_matrix1, taxa_are_rows = TRUE)

fldnaVirTMM <- phyloseq(OTU2, TAX, samples)

saveRDS(fldnaVirTMM, file="flRNAV_june25_TMM_phylo.rds")
flrnaVirTMM <- readRDS("flRNAV_june25_TMM_phylo.rds")

#prepare RAW  counts for phyloseq
counts2 <- read.csv("FL_RNA_virus_classseq_counts5_25.raw.counts.csv", row.names = 1, header = TRUE, check.names=FALSE, sep = ",")
counts_matrix2 <- as.matrix(counts2)
OTU3 = otu_table(counts_matrix2, taxa_are_rows = TRUE)

fldnaVirRAW <- phyloseq(OTU3, TAX, samples)

saveRDS(fldnaVirRAW, file="flRNAV_june25_RAW_phylo.rds")
flrnaVirRAW <- readRDS("flRNAV_june25_RAW_phylo.rds")
################################################### Analyses below

# removing failed sample “DT25_CNAT_1D_June_2022”
grep("DT28_MCAV_3D_Febr_2022", sample_names(flrnaVirRAW), value = TRUE)
to_drop <- "DT28_MCAV_3D_Febr_2022"
stopifnot(to_drop %in% sample_names(flrnaVirRAW))
ps_clean <- prune_samples(sample_names(flrnaVirRAW) != to_drop, flrnaVirRAW)
to_drop %in% sample_names(ps_clean)   # should return FALSE

# Compute α‐diversity (Shannon) for every sample
alpha_df <- estimate_richness(ps_clean, measures = "Shannon") %>%
  rownames_to_column(var = "sample")  # moves sample names into a column

# Merge sample‐level metadata (which must include a “time_point” column)
meta_df <- as(sample_data(flrnaVirRAW), "data.frame") %>%
  rownames_to_column(var = "sample")

alpha_df <- alpha_df %>%
  left_join(meta_df, by = "sample")

# 4. Quick boxplot of Shannon by time_point
ggplot(alpha_df, aes(x = time_point, y = Shannon, fill = time_point)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1) +
  labs(
    x = "Stage (T0 = Pre outbreak, T1 = Early outbreak, T2 = Late outbreak)",
    y = "Shannon diversity (viral transcripts)",
    title = "α‐Diversity across SCTLD stages"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

# Test for overall differences across time_point: Kruskal–Wallis
kw_overall <- kruskal.test(Shannon ~ time_point, data = alpha_df)
print(kw_overall)

## results
#Kruskal-Wallis rank sum test

#data:  Shannon by time_point
#Kruskal-Wallis chi-squared = 20.604, df = 2, p-value = 3.356e-05

pairwise.wilcox.test(
  x = alpha_df$Shannon,
  g = alpha_df$time_point,
  p.adjust.method = "BH"
)

##6. results
#Pairwise comparisons using Wilcoxon rank sum exact test 

#data:  alpha_df$Shannon and alpha_df$time_point 

# T0         T1   
# T1 0.032   -    
# T2 1.5e-05 0.017
#All three pairwise tests are significant after BH correction:

#T0 vs T1: p = 0.032

#T0 vs T2: p = 1.5 × 10⁻⁵

#T1 vs T2: p = 0.017

#In other words, Shannon diversity at T0 (pre‐outbreak) is significantly different 
#from during the outbreak (T1) and from post‐outbreak survivors (T2), and T1 also
#differs significantly from T2.2

#checking within species
#P value adjustment method: BH 
for(sp in unique(alpha_df$host_species)) {
  sub <- filter(alpha_df, host_species == sp)
  cat("\n--- Host species:", sp, "---\n")
  print(kruskal.test(Shannon ~ time_point, data = sub))
}

#--- Host species: natans ---
  
#  Kruskal-Wallis rank sum test

#data:  Shannon by time_point
#Kruskal-Wallis chi-squared = 0.43647, df = 2, p-value = 0.8039


#--- Host species: cavernosa ---
  
#  Kruskal-Wallis rank sum test

#data:  Shannon by time_point
#Kruskal-Wallis chi-squared = 15.625, df = 2, p-value = 0.0008174


#--- Host species: faveolata ---
  
 # Kruskal-Wallis rank sum test

#data:  Shannon by time_point
#Kruskal-Wallis chi-squared = 9.086, df = 2, p-value = 0.01064


#--- Host species: franksi ---
  
#  Kruskal-Wallis rank sum test

#data:  Shannon by time_point
#Kruskal-Wallis chi-squared = 3.5026, df = 2, p-value = 0.1736

#For C. natans and O. franksi, Shannon does not differ across T0/T1/T2 (p = 0.80 and p = 0.17, respectively). In contrast:
  
#  M. cavernosa: χ² = 15.63, p = 0.0004 → significant.

# O. faveolata: χ² = 9.09, p = 0.0106 → significant.

#  M. cavernosa:
mcav <- filter(alpha_df, host_species == "cavernosa")
pairwise.wilcox.test(mcav$Shannon, mcav$time_point, p.adjust.method = "BH")
#T0 vs T1: p = 0.2773 → not significant

#T0 vs T2: p = 0.0032 → significant

#T1 vs T2: p = 0.0009 → significant

#So Shannon diversity does not change from pre‐outbreak to outbreak (T0→T1), 
#but it does change between outbreak and survivor stages (T1→T2), and survivors 
#also differ from pre‐outbreak (T0→T2).

#  O. faveolata:
fave <- filter(alpha_df, host_species == "faveolata")
pairwise.wilcox.test(fave$Shannon, fave$time_point, p.adjust.method = "BH")

#T0 vs T1 (pre‐outbreak vs during outbreak): p = 0.1012 → not significant

#T1 vs T2 (during outbreak vs survivors): p = 0.4342 → not significant

#T0 vs T2 (pre‐outbreak vs survivors): p = 0.0035 → significant

#So Shannon diversity in O. faveolata doesn’t drop immediately when disease appears (T0→T1), 
#nor does it change significantly between the outbreak and survivor stages (T1→T2). However, 
#by the survivor stage (T2), diversity is significantly different from the pre-outbreak 
#baseline (T0). This suggests that for O. faveolata, the community shift occurs sometime 
#between T1 and T2, rather than right when the disease hits.

#I am interesting in the impact of disease on alpha diversity or the interaction between disease and time point

alpha_df$time_point   <- factor(alpha_df$time_point, levels = c("T0","T1","T2"))
alpha_df$health_status <- factor(alpha_df$health_status)

alpha_df$health_status <- factor(alpha_df$health_status)

kw_health <- kruskal.test(Shannon ~ health_status, data = alpha_df)
print(kw_health)
#Kruskal-Wallis rank sum test

#data:  Shannon by health_status
#Kruskal-Wallis chi-squared = 1.0973, df = 2, p-value = 0.5777

alpha_df <- estimate_richness(ps_clean, measures = "Shannon") %>%
  rownames_to_column(var = "sample") %>%
  left_join(
    as(sample_data(flrnaVirRAW), "data.frame") %>% rownames_to_column("sample"),
    by = "sample"
  ) %>%
  mutate(
    time_point    = factor(time_point, levels = c("T0","T1","T2")),
    Site          = factor(Site,       levels = c("SCTLD25","SCTLD26","SCTLD28")),
    host_species  = factor(host_species, levels = c("cavernosa","faveolata","franksi","natans"))
  )


comparisons <- list(
  c("T0", "T1"),
  c("T1", "T2"),
  c("T0", "T2")
)

p <- ggplot(alpha_df, aes(x = time_point, y = Shannon, fill = Site)) +
  geom_boxplot(
    position      = position_dodge(width = 0.75),
    width         = 0.5,
    color         = "black",     # black border
    outlier.shape = NA,          # hide the default outlier points (we’ll overlay jitter)
    alpha         = 0.8
  ) +
  
  geom_jitter(
    aes(color = Site),
    position = position_jitterdodge(
      jitter.width = 0.15,       # horizontal “spread” so points don’t stack
      dodge.width  = 0.75        # match the boxplot dodge width
    ),
    size  = 2,
    alpha = 0.9
  ) +
  
  stat_compare_means(
    aes(label = ..p.signif..),    # use “ns”/“*”/“**” notation automatically
    method          = "wilcox.test",
    comparisons     = comparisons,
    p.adjust.method = "BH",
    hide.ns         = FALSE,      # still plot “ns” if not significant
    size            = 4,          # font size of significance labels
    tip.length      = 0.01,       # short “tips” on the bar ends
    step.increase   = 0.10        # vertical spacing between successive comparisons
  ) +

  facet_wrap(~ host_species, nrow = 2, scales = "free_y") +
  
  scale_fill_manual(
    values = c(
      "SCTLD25" = "#F8766D",     
      "SCTLD26" = "#00BA38",     
      "SCTLD28" = "#619CFF"     
    )
  ) +
  labs(
    x     = "Time Point (T0 = Pre outbreak, T1 = Early outbreak, T2 = Late outbreak)",
    y     = "Shannon Diversity Index",
    fill  = "Site",
    color = "Site",
    title = "Viral Shannon Diversity Across Time Points"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(

    legend.position = "bottom",
    
    strip.background = element_rect(fill = "grey80", color = NA),
    strip.text       = element_text(face = "bold", size = 13),
    
    axis.text.x      = element_text(size = 12),
    axis.text.y      = element_text(size = 12),

    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", size = 0.3),
    
    panel.border = element_rect(color = "black", fill = NA, size = 0.3),

    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  )

print(p)


####################

#BETA DIVErSITY 

# Remove outlier
to_drop <- "DT28_MCAV_3D_Febr_2022"
if (!to_drop %in% sample_names(flrnaVirTMM)) stop("Sample not found")
flrnaVirTMM_clean <- prune_samples(sample_names(flrnaVirTMM) != to_drop, flrnaVirTMM)

# Relative abundance
ps_rel <- transform_sample_counts(flrnaVirTMM_clean, function(x) x / sum(x))

# Bray–Curtis distance
bc_dist <- phyloseq:::distance(ps_rel, method="bray")

meta_df <- as(sample_data(ps_rel), "data.frame") %>%

  tibble::rownames_to_column(var = "sample") %>%

  mutate(
    host_species  = factor(host_species,  levels = c("cavernosa", "faveolata", "franksi", "natans")),
    time_point    = factor(time_point,    levels = c("T0", "T1", "T2")),
    health_status = factor(health_status, levels = c("healthy", "diseased", "apparently_healthy", "quiesced")),
    Site          = factor(Site,          levels = c("SCTLD25", "SCTLD26", "SCTLD28")),
    Month         = factor(Month,         levels = c("February", "June", "August", "September")),
    Year          = factor(as.character(Year)),
    tissue_type   = factor(tissue_type)
  ) %>%

  select(sample, host_species, time_point, health_status, Site, Month, Year, tissue_type)

# 5) Check dispersion
print(permutest(betadisper(bc_dist, meta_df$time_point), permutations=999))
print(permutest(betadisper(bc_dist, meta_df$Site), permutations=999))
print(permutest(betadisper(bc_dist, meta_df$health_status), permutations=999))
##RESULTS
#All three tests return p ≫ 0.05, which means group dispersions do not differ for 
#time_point, Site, or health_status. In other words, there is no evidence of 
#unequal multivariate variance that could bias a PERMANOVA. You can therefore 
#proceed to test for differences in centroid location (i.e. “true” β-diversity 
#differences) without worrying about dispersion artifacts.

# Overall PERMANOVA: species + time + health (to check how species influences)
set.seed(123)
adonis_species_time_health <- adonis2(
  bc_dist ~ host_species + time_point + health_status,
  data        = meta_df,
  permutations = 999,
  by          = "margin"
)
print(adonis_species_time_health)

run_permanova_per_species <- function(spec, bc_dist, meta_df) {
  meta_sub <- meta_df %>% filter(host_species == spec)
  mat      <- as.matrix(bc_dist)
  keep     <- rownames(mat) %in% meta_sub$sample
  mat2     <- mat[keep, keep]
  bc_sub   <- as.dist(mat2)
  
  cat("\n=== Species:", spec, "===\n")
  print(permutest(betadisper(bc_sub, meta_sub$time_point), permutations = 999))
  print(permutest(betadisper(bc_sub, meta_sub$health_status), permutations = 999))
  
  set.seed(123)
  ad <- adonis2(
    bc_sub ~ time_point + health_status,
    data        = meta_sub,
    permutations = 999,
    by          = "margin"
  )
  print(ad)
}

species_list <- levels(meta_df$host_species)
lapply(species_list, function(sp) run_permanova_per_species(sp, bc_dist, meta_df))

# --- Beta diversity PCoA visualizations facet by species

# Overall PCoA colored/shaped by time_point & health_status
p_all <- plot_ordination(ps_rel, ordinate(ps_rel, method="PCoA", distance=bc_dist),
                         color="time_point", shape="health_status") +
  geom_point(size=3, alpha=0.8) +
  scale_color_manual(values=c("T0"="#1B9E77","T1"="#D95F02","T2"="#7570B3")) +
  scale_shape_manual(values=c(16, 17, 15, 1)) +
  labs(title="Overall PCoA (color = Time, shape = Health)") +
  theme_minimal(base_size=14)
print(p_all)

# Faceted PCoA by species
p_facet <- plot_ordination(ps_rel, ordinate(ps_rel, method="PCoA", distance=bc_dist),
                           color="time_point", shape="health_status") +
  geom_point(size=3, alpha=0.8) +
  facet_wrap(~ host_species, nrow=2, scales="free") +
  scale_color_manual(values=c("T0"="#1B9E77","T1"="#D95F02","T2"="#7570B3")) +
  scale_shape_manual(values=c(16, 17, 15, 1)) +
  labs(title="PCoA by Species (color = Time, shape = Health)") +
  theme_minimal(base_size = 14) +
  theme(
    strip.text            = element_text(face = "bold", size = 14),     # bold facet labels
    axis.text             = element_text(size = 12),                    # tick labels
    axis.title            = element_text(face = "bold", size = 14),     # axis titles
    panel.grid.major      = element_line(color = "grey90"),             # light grid lines
    panel.grid.minor      = element_blank(),
    legend.title          = element_text(face = "bold", size = 12),     # legend titles bold
    legend.text           = element_text(size = 11),
    panel.border          = element_rect(color = "black", fill = NA),   # black border around each panel
    axis.line             = element_line(color = "black")               # black axis lines
  )
print(p_facet)

#### CORE VIROME PER SPECIES #####################################################

to_drop <- "DT28_MCAV_3D_Febr_2022"
if (!to_drop %in% sample_names(flrnaVirTMM)) stop("Sample not found")
flrnaVirTMM_clean <- prune_samples(sample_names(flrnaVirTMM) != to_drop, flrnaVirTMM)

sample_data(flrnaVirTMM_clean) %>% colnames()

species_list <- c("cavernosa", "faveolata", "franksi", "natans")

core_by_species <- vector("list", length(species_list))
names(core_by_species) <- species_list

tt_all   <- tax_table(flrnaVirTMM_clean)
classVec <- as.character(tt_all[, "Class"])

keep_idx <- ! (classVec %in% "c__unclassified" | classVec %in% "c_unclassified")

ps_filt  <- prune_taxa(keep_idx, flrnaVirTMM_clean)

rawscore_vec <- as.numeric(as.character(tax_table(ps_filt)[, "rawscore"]))

keep_taxa <- rawscore_vec >= 2


ps_high_rawscore <- prune_taxa(keep_taxa, ps_filt)

ps_noU_class <- prune_taxa(
  !grepl("_U$", tax_table(ps_high_rawscore)[, "Class"]),
  ps_high_rawscore
)

ord_vec <- as.character(tax_table(ps_noU_class)[, "Class"])

keep_viricetes <- grepl("viricetes$", ord_vec, ignore.case = TRUE)

ps_viricetes <- prune_taxa(keep_viricetes, ps_noU_class)

message("Kept ", ntaxa(ps_high_rawscore), " taxa (out of ", ntaxa(ps_viricetes), ") with rawscore ≥ 2")


  for (sp in species_list) {
    message("→ Processing species: ", sp)

    ps_sp <- subset_samples(ps_viricetes, host_species == sp)
    
    # Extract raw OTU counts and Class labels
    otu_mat  <- as(otu_table(ps_sp), "matrix")
    if (!taxa_are_rows(ps_sp)) {
      otu_mat <- t(otu_mat)
    }
    class_vec <- as.character(tax_table(ps_sp)[, "Class"])
    
    # Collapse all OTUs by Class via rowsum
    class_counts <- rowsum(otu_mat, group = class_vec)
    
    # Convert to presence/absence and compute prevalence
    class_pa  <- (class_counts > 0) * 1
    prev_vec  <- rowSums(class_pa) / ncol(class_pa)
    
    # Core classes = those in ≥ 95% of samples
    core_classes <- names(prev_vec)[prev_vec >= 0.95]
    core_by_species[[sp]] <- core_classes
    
    message("   • Found ", length(core_classes), " core classes for ", sp)
  }
  
  core_by_species
  
  df_long <- data.frame(
    species    = rep(names(core_by_species),
                     times = sapply(core_by_species, length)),
    core_class = unlist(core_by_species, use.names = FALSE),
    stringsAsFactors = FALSE
  )
  
  # 2) write to CSV
  write.csv(df_long,
            file      = "core_by_species_long_onlyreal.csv",
            row.names = FALSE)
  
##visualize 
  all_core_classes <- sort(unique(unlist(core_by_species)))

  library(phyloseq)
  tt_all <- tax_table(ps_filt)
  class_vec_all <- as.character(tt_all[, "Class"])
  keep_tx <- rownames(tt_all)[ class_vec_all %in% all_core_classes ]

  ps_core_tx <- prune_taxa(keep_tx, ps_filt)

  ####Cladogram####################

  samdf <- data.frame(sample_data(ps_core_tx))
  species_list <- c("cavernosa","faveolata","franksi","natans")
  
  # Extract relative‐abundance OTU table as a matrix
  otu_all <- as(otu_table(ps_core_tx), "matrix")
  if (!taxa_are_rows(ps_core_tx)) otu_all <- t(otu_all)
  
  # Turn counts into relative abundance per sample
  rel_abund_mat <- sweep(otu_all, 2, colSums(otu_all), FUN = "/")
  
  # Initialize a 4×N matrix
  taxa_names_core <- taxa_names(ps_core_tx)
  sp_abund <- matrix(
    0, 
    nrow = length(species_list),
    ncol = length(taxa_names_core),
    dimnames = list(species_list, taxa_names_core)
  )
  
  for (sp in species_list) {
    # get sample names for that species
    sam_sp <- rownames(samdf)[samdf$host_species == sp]
    # subset the relative-abund matrix to those samples
    sub_mat <- rel_abund_mat[, sam_sp, drop = FALSE]
    # take the mean relative-abundance of each transcript across all those samples
    sp_abund[sp, ] <- rowMeans(sub_mat)
  }
  
  # Compute Bray–Curtis distance on that 4×N abundance matrix
  dist_bc_sp <- vegan::vegdist(sp_abund, method = "bray")
  hc_bc_sp   <- hclust(dist_bc_sp, method = "average")
  
  # Plot the resulting tree
  plot(
    hc_bc_sp,
    main = "Clustering of Coral Species by Core Class Transcript Abundance\n(Bray–Curtis)",
    sub  = "",
    xlab = "Species",
    ylab = "Bray–Curtis distance",
    hang = -1,
    cex  = 1.2
  )
  
  write.csv(df_export, "core_class_by_species_abundance.csv", row.names = FALSE)
  

###### Diff expression 
  
  to_drop <- "DT28_MCAV_3D_Febr_2022"
  if (!to_drop %in% sample_names(flrnaVirRAW)) stop("Sample not found")
  flrnaVirRAW_clean <- prune_samples(sample_names(flrnaVirRAW) != to_drop, flrnaVirRAW)
  
  rawscore_vec <- as.numeric(as.character(tax_table(flrnaVirRAW_clean)[, "rawscore"]))

  keep_taxa <- rawscore_vec >= 2

  flrnaVirRAW_clean_rawscore <- prune_taxa(keep_taxa, flrnaVirRAW_clean)

  message("Kept ", ntaxa(flrnaVirRAW_clean_rawscore), " taxa (out of ", ntaxa(flrnaVirRAW_clean), ") with rawscore ≥ 2")
  
  
  ###health
  
  diagdds = phyloseq_to_deseq2(flrnaVirRAW_clean_rawscore, ~ host_species + tissue_type)
  
  diagdds$host_species <- factor(diagdds$host_species,
                                     levels = c("cavernosa", "faveolata", "franksi", "natans"))
  diagdds$tissue_type  <- factor(diagdds$tissue_type,
                                     levels = c("healthy", "disease_margin", "apparently_healthy"))
  
  dds_local <- DESeq(diagdds, fitType = "local")
  
  resLocal1 <- results(dds_local)
  alpha = 0.1
  resultsNames(dds_local)
  results(dds_local)
  
  alpha_cutoff <- 0.1  # FDR threshold
  
  res_disease <- results(
    dds_local,
    name  = "tissue_type_disease_margin_vs_healthy",
    alpha = alpha_cutoff
  )
  
  summary(res_disease)
  res_dh <- results(dds_local,
                    contrast = c("tissue_type","disease_margin","healthy"),
                    alpha    = 0.1)
  res_da <- results(dds_local,
                    contrast = c("tissue_type","disease_margin","apparently_healthy"),
                    alpha    = 0.1)
  

  # Extract up‐regulated sets (padj<0.1 & log2FC>0)
  up_dh <- rownames(subset(res_dh, padj < 0.1 & log2FoldChange > 0))
  up_da <- rownames(subset(res_da, padj < 0.1 & log2FoldChange > 0))
  
  # Combine into one list
  both_up   <- intersect(up_dh, up_da)  
  either_up <- union(up_dh, up_da)      
  
  length(both_up)  
  length(either_up) 
  
  final_tx <- both_up   
  
  writeLines(both_up, con = "UPregulated_in DM_RNA_seq.txt")
  
  
#### visualize deSEq2 results 

  to_drop <- "DT28_MCAV_3D_Febr_2022"
  if (!to_drop %in% sample_names(flrnaVirTMM)) stop("Sample not found")
  flrnaVirTMM_clean <- prune_samples(sample_names(flrnaVirTMM) != to_drop, flrnaVirTMM)
  ps_final <- prune_taxa(final_tx, flrnaVirTMM_clean)
  
 
  #all by class
  subset_classified_up <- subset_taxa(ps_final, tax_table(ps_final)[, "Class"] != "c__unclassified")
  disease_margin_samples <- subset_samples(subset_classified_up, tissue_type == "disease_margin")
  
  h3 <- abundance_heatmap(disease_margin_samples, classification = 'Class',
                          treatment = "host_species", transformation = 'log10')
  h3 <- as.ggplot(h3)
  h3
  
    ps_noU_order <- prune_taxa(
    !grepl("_U$", tax_table(ps_final)[, "Order"]),
    ps_final
  )
  ord_vec <- as.character(tax_table(ps_noU_order)[, "Order"])
  
  keep_virales <- grepl("virales$", ord_vec, ignore.case = TRUE)

  ps_virales <- prune_taxa(keep_virales, ps_noU_order)
  #all
  subset_classified_up <- subset_taxa(ps_virales, tax_table(ps_virales)[, "Order"] != "o__unclassified")
  disease_margin_samples <- subset_samples(subset_classified_up, tissue_type == "disease_margin")
  
  h3 <- abundance_heatmap(disease_margin_samples, classification = 'Order',
                          treatment = "host_species", transformation = 'log10')
  h3 <- as.ggplot(h3)
  h3
  

###TEM 

library(readxl)
library(lme4)
library(lmerTest)   
library(emmeans)    

# 1) Import data
dat <- read_excel("TEM_DRTO_2.xlsx", sheet = "Sheet1")

tiss  <- dat %>%
  filter(tolower(Type)=="tissue") %>%
  transmute(
    Species,
    Health,
    p_tissue = as.numeric(`Percent_tissue_images_showing_VLPS`),
    asin_t   = asin(sqrt(p_tissue))
  )

zoox <- dat %>%
  filter(startsWith(Type, "Zoox")) %>%
  mutate(
    p_zoox = as.numeric(`Percent_of_zoox_cells_images_showing_VLPS`)
  ) %>%
  filter(!is.na(p_zoox))

# —— POOLED‐SPECIES ANALYSES —— #
# Tissue two‐way ANOVA
mod_tiss <- lm(asin_t ~ Health * Species, data = tiss)
anova(mod_tiss)
#Analysis of Variance Table

#Response: asin_z
#Df  Sum Sq Mean Sq F value Pr(>F)
#Health          2 0.95763 0.47882  1.2985 0.3514
#Species         1 0.28818 0.28818  0.7815 0.4172
#Health:Species  2 0.08969 0.04484  0.1216 0.8880
#Residuals       5 1.84377 0.36875   
# Zoox mixed‐effects model (arcsine‐sqrt transform)
zoox2 <- zoox %>% mutate(asin_z = asin(sqrt(p_zoox)))

mod_pooled <- lmer(asin_z ~ Health + (1|Species/ID), data = zoox2)
anova(mod_pooled)                      # F‐test for Health
#Type III Analysis of Variance Table with Satterthwaite's method
#        Sum Sq  Mean Sq NumDF  DenDF F value Pr(>F)
#Health 0.16182 0.080909     2 6.9466   1.481 0.2913

# —— CNAT‐ONLY ANALYSES —— #

tiss_cnat <- filter(tiss, Species == "Cnat")
zoox_cnat <- filter(zoox, Species == "Cnat")

tiss_cnat %>% kruskal_test(p_tissue ~ Health)
# A tibble: 1 × 6
#.y.          n statistic    df     p method        
#* <chr>    <int>     <dbl> <int> <dbl> <chr>         
#  1 p_tissue     7     0.769     2 0.681 Kruskal-Wallis

zoox_cnat %>% kruskal_test(p_zoox ~ Health)
# A tibble: 1 × 6
#.y.        n statistic    df     p method        
# <chr>  <int>     <dbl> <int> <dbl> <chr>         
#  1 p_zoox    23      4.33     2 0.115 Kruskal-Wallis

# 7) Cnat zoox mixed‐effects model
zoox_cnat2 <- zoox_cnat %>% mutate(asin_z = asin(sqrt(p_zoox)))

mod_cnat <- lmer(asin_z ~ Health + (1|ID), data = zoox_cnat2)
anova(mod_cnat)                        # F‐test for Health in Cnat
#Type III Analysis of Variance Table with Satterthwaite's method
#        Sum Sq  Mean Sq NumDF  DenDF F value Pr(>F)
#Health 0.04769 0.023845     2 3.0433  0.3721 0.7169


#visualize

tem <- read_csv("TEM_DRTO_2.csv", show_col_types = FALSE)
names(tem) <- str_trim(names(tem))

tem <- tem %>% 
  mutate(
    Health  = factor(Health, levels = c("H", "AH", "D"),
                     labels = c("Healthy", "Apparently healthy", "Diseased")),
    Species = factor(Species)
  )

tem <- tem %>%
  mutate(
    `Percent_tissue_images_showing_VLPS` = as.numeric(`Percent_tissue_images_showing_VLPS`),
    `Percent_of_zoox_cells_images_showing_VLPS` = as.numeric(`Percent_of_zoox_cells_images_showing_VLPS`)
  )

tissue_df <- filter(tem, Type == "Tissue")
zoox_df   <- filter(tem, str_starts(Type, "Zoox"))

# Figure Tissue‐level VLP prevalence
p1 <- ggplot(tissue_df, 
             aes(x = Health, 
                 y = `Percent_tissue_images_showing_VLPS`, 
                 fill = Health)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.size = 2) +
  facet_wrap(~ Species) +
  scale_fill_manual(values = c(
    "Healthy" = "#00BA38",
    "Apparently healthy" = "#619CFF",
    "Diseased" = "#F8766D"
  )) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Tissue VLP Prevalence by Health & Species",
    x = "Health State",
    y = "% Tissue Images with VLPs"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

# Figure Zoox‐cell VLP prevalence
p2 <- ggplot(zoox_df, 
             aes(x = Health, 
                 y = `Percent_of_zoox_cells_images_showing_VLPS`, 
                 fill = Health)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.size = 2) +
  facet_wrap(~ Species) +
  scale_fill_manual(values = c(
    "Healthy" = "#00BA38",
    "Apparently healthy" = "#619CFF",
    "Diseased" = "#F8766D"
  )) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Zoox Cell VLP Prevalence by Health & Species",
    x = "Health State",
    y = "% Zoox Cells with VLPs"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

# Print figures
print(p1)
print(p2)

##morphology VLP per health


dat <- read_csv("TEM_DRTO_2.csv")
names(dat)

plot_df <- dat %>%
  filter(Type != "Summary") %>%
  
  mutate(percent = case_when(
    tolower(Type) == "tissue"   ~ as.numeric(`Percent_tissue_images_showing_VLPS`),
    grepl("^Zoox", Type)        ~ as.numeric(`Percent_of_zoox_cells_images_showing_VLPS`),
    TRUE                         ~ NA_real_
  )) %>%
  drop_na(percent) %>%
  
  pivot_longer(
    cols      = c(`Icosahedral (1 or 0)`, `Filamentous (1 or 0)`),
    names_to  = "morphology",
    values_to = "present"
  ) %>%
  filter(present == 1) %>%
  
  mutate(
    morphology = recode(morphology,
                        `Filamentous (1 or 0)` = "filamentous",
                        `Icosahedral (1 or 0)` = "icosahedral"),
    morphology = factor(morphology, levels = c("filamentous","icosahedral")),
    Species    = factor(Species,    levels = c("Cnat","Mcav")),
    Health     = factor(Health,     levels = c("H","AH","D"))
  )

# Plot morphology by species/health
ggplot(plot_df, aes(x = Health, y = percent, fill = Health)) +
  geom_boxplot() +
  facet_grid(
    rows   = vars(morphology),
    cols   = vars(Species),
    switch = "y"               # moves the y‐strip text outside
  ) +
  scale_y_continuous(labels = percent_format(1), limits = c(0,1)) +
  scale_fill_manual(values = c(
    "AH" = "#619CFF",
    "D"  = "#F8766D",
    "H"  = "#00BA38"
  )) +
  labs(x = "health", y = "percent_images_showing_vlps") +
  theme_bw() +
  theme(
    strip.placement  = "outside",
    strip.background = element_rect(fill = "grey80", colour = NA),
    legend.position  = "none",
    panel.grid.major = element_line(colour = "grey90")
  )
