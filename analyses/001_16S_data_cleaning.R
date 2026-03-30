
# Script parameters ----

#Remove saved elements
rm(list = ls())
gc()

# Library
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("phyloseq", quietly = TRUE))
  BiocManager::install("phyloseq")

library("BiocManager")
library("phyloseq")
library("dplyr")
library("ggplot2")

#Function to track number of samples and reads at each 
track_phyloseq <- function(physeq, step_name = NA, step_number = NA) {
  data.frame(
    step_name = step_name,
    step_number = step_number,
    n_taxa = phyloseq::ntaxa(physeq),
    n_samples = phyloseq::nsamples(physeq),
    total_reads = sum(phyloseq::taxa_sums(physeq))
  )
}

## List of rodents identified putative infectious agents ----
#Beware, using this list should not replace prior exploration of the cluster present in the metabarcoding dataset
putative_infectious_genra <- c("Bartonella",
"Neoehrlichia" ,
"Anaplasma" ,
"Francisella", 
"Rickettsia" ,
"Orientia" ,
"Chlamydia" ,
"Mycoplasma",
"Brevinema",
"Leptospira", 
"Borrelia",
"Streptococcus",
"Fusobacterium",
"Neisseria",
"Arcobacter",
"Helicobacter",
"Yersinia",
"Rickettsiella",
"Spiroplasma",
"Coxiella",
"Occidentia",
"Brucella" ,
"Ehrlichia",
"Streptobacillus",
"Pseudostreptobacillus",
"Rodentibacter",
"Treponema",
"Pasteurella",
"Campylobacter",
"Mycobacterium")

putative_infectious_family <- c("Mycoplasmataceae")

putative_infectious_phylum <- c("Apicomplexa")


ranks <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Import data ----

# Import Metabarcoding data
s16_abundance <- data.table::fread(
  file = here::here("data", "raw-data", "diagnostic", "2026.03.24-16S.abundance.final.txt"))
s16_taxo <- data.table::fread(
  file = here::here("data", "raw-data", "diagnostic", "2026.03.24-16S.taxonomy.final.txt"))

#Ensure rowname for metabarcoding data
rownames(s16_abundance) <- s16_abundance$V1
s16_abundance$V1 <- NULL
rownames(s16_taxo) <- s16_taxo$V1
s16_taxo$V1 <- NULL

#Samples data
individuals <- read.delim(file = here::here("data", "raw-data", "20260320_BePrep_samples.csv"), sep = ";")


# Data object creation  ----

## Keep relevant samples  ----

# Keep relevant sample type (Unsure we have only .SP (spleen) samples)
s16_abundance %>% select(!stringr::str_subset(string = colnames(s16_abundance), pattern = ".SP$")) %>% colnames()
s16_abundance <- s16_abundance %>% 
  select(stringr::str_subset(string = colnames(s16_abundance), pattern = ".SP$"))

# Make names concordant with individuals table
s16_abundance <- s16_abundance %>%
  rename_with(~ stringr::str_remove(., "\\.SP$"))

# Keep only BEPREP project samples
s16_abundance <- s16_abundance %>%
  select(dplyr::any_of(individuals %>% dplyr::pull(numero_centre)))



## Generate phyloseq objects  ----

#Convert abundance dataframe to matrix
matrix_ab_s16 <- s16_abundance %>% as.matrix()
rownames(matrix_ab_s16) <- rownames(s16_abundance)

#Convert taxonomy dataframe to matrix
matrix_tax_s16 <- s16_taxo %>%
  tidyr::separate_wider_delim(cols = "blast_taxonomy", delim = ";", names = ranks ) %>%
  as.matrix()
rownames(matrix_tax_s16) <- rownames(s16_taxo)

#Prepare metadata (samples)
#Make numero_centre as rownames for individuals dataframe
individuals <- individuals %>%
  tibble::column_to_rownames(var = "numero_centre")

individuals$numero_centre <- rownames(individuals)


#Convert matrix and metadata to phyloseq objects
OTU_s16 <- phyloseq::otu_table(matrix_ab_s16, taxa_are_rows = TRUE)
TAX_s16 <- phyloseq::tax_table(matrix_tax_s16 )

SAMPL <- phyloseq::sample_data(individuals)


#Generate phyloseq object
physeq_s16 <- phyloseq::phyloseq(OTU_s16, TAX_s16, SAMPL)

#Recount sum of reads for clusters across samples
tax <- as(tax_table(physeq_s16), "matrix") %>% 
  cbind(samples_sum = as.numeric(phyloseq::taxa_sums(physeq_s16)))
tax_table(physeq_s16) <- tax_table(tax)

#Initialize tracking table
tracking_s16 <- track_phyloseq(physeq_s16, step_name = "initial", step_number = 0)



# Filtering ----
#After every filtering steps, samples not containing reads are removed in order
#to follow filtering effect on sample number

#Initialize tracking table
physeq_s16 <- phyloseq::prune_samples(
  x = physeq_s16,
  samples = sample_sums(physeq_s16) != 0)
tracking_s16 <- track_phyloseq(physeq_s16, step_name = "Beprep_samples_with_reads", step_number = 0)


## Remove clusters with 0 reads ----
physeq_s16 <- phyloseq::prune_taxa(
  taxa = phyloseq::taxa_sums(physeq_s16) > 0,
  x =  physeq_s16)

physeq_s16 <- phyloseq::prune_samples(
  x = physeq_s16,
  samples = sample_sums(physeq_s16) != 0)
tracking_s16 <- bind_rows(
  tracking_s16, track_phyloseq(physeq_s16, step_name = "Removed_0_reads_clusters", step_number = 1))


## Remove clusters without affiliation ---- 
physeq_s16 <- phyloseq::prune_taxa( 
  taxa = as.vector(phyloseq::tax_table(physeq_s16)[, "Domain"]) != "no",
  x = physeq_s16)

physeq_s16 <- phyloseq::prune_samples(
  x = physeq_s16,
  samples = sample_sums(physeq_s16) != 0)
tracking_s16 <- bind_rows(
  tracking_s16, track_phyloseq(physeq_s16, step_name = "Removed_NA_clusters", step_number = 2))


## Remove clusters that are not Bacteria or Apicomplexa ---- 
physeq_s16 <- phyloseq::prune_taxa( 
  taxa = as.vector(phyloseq::tax_table(physeq_s16)[, "Domain"]) %in% c("Bacteria", "Apicomplexa"),
  x = physeq_s16)

physeq_s16 <- phyloseq::prune_samples(
  x = physeq_s16,
  samples = sample_sums(physeq_s16) != 0)
tracking_s16 <- bind_rows(
  tracking_s16, track_phyloseq(physeq_s16, step_name = "Removed_not_Bacteria-Apicomplexa_clusters", step_number = 3))


## Identify best affiliation level ---- 

#Extract taxonomy matrix
tax_full <- as(tax_table(physeq_s16), "matrix") 

#Keep only ranks 
tax <- tax_full[, ranks]

#Clean ambiguous annotations (make them NA)
tax[grepl("unknown|Unknown|multi-affiliation|Multi-affiliation", tax)] <- NA

#Reverse matrix to go from Species to Domain
tax_rev <- tax[, rev(seq_len(ncol(tax))), drop = FALSE]
ranks_rev <- rev(ranks)

#Find first non-NA per row
first_valid_index <- max.col(!is.na(tax_rev), ties.method = "first")

#In case some rows contain only NA 
all_na <- rowSums(!is.na(tax_rev)) == 0

#Extract best taxa and level
best_taxa <- tax_rev[cbind(seq_len(nrow(tax_rev)), first_valid_index)]
best_taxa_level <- ranks_rev[first_valid_index]

best_taxa[all_na] <- NA
best_taxa_level[all_na] <- NA

#Replace original rank columns with cleaned ones
tax_full[, ranks] <- tax

#Insert new columns after last rank column
last_rank_pos <- max(match(ranks, colnames(tax_full)))

tax_final <- cbind(
  tax_full[, 1:last_rank_pos, drop = FALSE],
  best_taxa = best_taxa,
  best_taxa_level = best_taxa_level,
  tax_full[, (last_rank_pos + 1):ncol(tax_full), drop = FALSE]
)

#Back to phyloseq
tax_table(physeq_s16) <- tax_table(tax_final)


## Remove clusters that are affiliated above family level ---- 
physeq_s16 <- phyloseq::prune_taxa( 
  taxa = as.vector(phyloseq::tax_table(physeq_s16)[, "best_taxa_level"]) %in% c("Family", "Genus", "Species"),
  x = physeq_s16)

physeq_s16 <- phyloseq::prune_samples(
  x = physeq_s16,
  samples = sample_sums(physeq_s16) != 0)
tracking_s16 <- bind_rows(
  tracking_s16, track_phyloseq(physeq_s16, step_name = "Removed_above_family_affiliation", step_number = 4))


## Remove clusters that are affiliated above genus level ---- 
physeq_s16 <- phyloseq::prune_taxa( 
  taxa = as.vector(phyloseq::tax_table(physeq_s16)[, "best_taxa_level"]) %in% c("Genus", "Species"),
  x = physeq_s16)

physeq_s16 <- phyloseq::prune_samples(
  x = physeq_s16,
  samples = sample_sums(physeq_s16) != 0)
tracking_s16 <- bind_rows(
  tracking_s16, track_phyloseq(physeq_s16, step_name = "Removed_above_genus_affiliation", step_number = 5))


## Remove clusters with reads below 0.001 % of total ---- 
physeq_s16 <- phyloseq::prune_taxa(
  taxa = phyloseq::taxa_sums(physeq_s16) >= (sum(phyloseq::taxa_sums(physeq_s16)) * 0.00001),
  x =  physeq_s16)

physeq_s16 <- phyloseq::prune_samples(
  x = physeq_s16,
  samples = sample_sums(physeq_s16) != 0)
tracking_s16 <- bind_rows(
  tracking_s16, track_phyloseq(physeq_s16, step_name = "Removed_below_0.001%_reads", step_number = 5))


## Remove clusters that are not from putative infectious taxa  ----

#At this step we are using a list of putative genera, family and phylum known to
#contain infectious agents for rodents and/or humans. Groups containing opportunistic 
#agents were not kept as most of them are natural components of environment and microbiota,
#contamination could not be excluded

physeq_s16 <- phyloseq::prune_taxa( 
  taxa = as.vector(phyloseq::tax_table(physeq_s16)[, "Genus"]) %in% putative_infectious_genra | 
    as.vector(phyloseq::tax_table(physeq_s16)[, "Family"]) %in% putative_infectious_family |
    as.vector(phyloseq::tax_table(physeq_s16)[, "Phylum"]) %in% putative_infectious_phylum,
  x = physeq_s16)

physeq_s16 <- phyloseq::prune_samples(
  x = physeq_s16,
  samples = sample_sums(physeq_s16) != 0)
tracking_s16 <- bind_rows(
  tracking_s16, track_phyloseq(physeq_s16, step_name = "Removed_not_infectious_clusters", step_number = 6))



## Remove clusters that may not be infectious ----

#At this step we are taking a look at species and genus kept in order to further
#remove taxa that may be contaminants (from environment or microbiota)

checked_genera_to_remove <- c("Campylobacter", 
                              "Coxiella", 
                              "Fusobacterium", 
                              "Helicobacter",
                              "Mycobacterium",
                              "Neisseria",
                              "Rickettsiella",
                              "Rodentibacter",
                              "Spiroplasma",
                              "Streptococcus",
                              "Treponema")

physeq_s16 <- phyloseq::prune_taxa(
  taxa = !(as.vector(phyloseq::tax_table(physeq_s16)[, "Genus"]) %in% checked_genera_to_remove),
  x = physeq_s16
)

physeq_s16 <- phyloseq::prune_samples(
  x = physeq_s16,
  samples = sample_sums(physeq_s16) != 0)
tracking_s16 <- bind_rows(
  tracking_s16, track_phyloseq(physeq_s16, step_name = "Removed_ambiguous_infectious_clusters", step_number = 7))


## Manual affiliation cleaning ----

#Extract taxonomy matrix
tax_full <- as(tax_table(physeq_s16), "matrix") 


#Rename some affiliations
tax_full[tax_full == "Pseudostreptobacillus"] <- "Streptobacillus"
tax_full[tax_full == "Borreliella"] <- "Borrelia"
tax_full[tax_full == "Borreliella afzelii"] <- "Borrelia afzelii"
tax_full[tax_full == "Mycoplasma coccoides"] <- "Eperythrozoon coccoides"
tax_full[tax_full == "Mycoplasma moatsii"] <- "Mesomycoplasma moatsii"


#Improve  affiliations of clusters
tax_full[tax_full[, "observation_name"] == "Cluster_3", "blast_perc_identity"] <- "99.2"

tax_full[tax_full[, "observation_name"] == "Cluster_5", "Species"] <- "Neoehrlichia mikurensis"

tax_full[tax_full[, "observation_name"] == "Cluster_11", "blast_perc_identity"] <- "91.24"

tax_full[tax_full[, "observation_name"] == "Cluster_23", "Species"] <- "Eperythrozoon coccoides"
tax_full[tax_full[, "observation_name"] == "Cluster_23", "blast_perc_identity"] <- "98"

tax_full[tax_full[, "observation_name"] == "Cluster_109", "Species"] <- "Eperythrozoon coccoides"
tax_full[tax_full[, "observation_name"] == "Cluster_109", "blast_perc_identity"] <- "94.8"

tax_full[tax_full[, "observation_name"] == "Cluster_137", "Species"] <- "Haemobartonella muris"
tax_full[tax_full[, "observation_name"] == "Cluster_137", "blast_perc_identity"] <- "92.4"

tax_full[tax_full[, "observation_name"] == "Cluster_519", "Species"] <- "Haemobartonella muris"
tax_full[tax_full[, "observation_name"] == "Cluster_519", "blast_perc_identity"] <- "94.8"

tax_full[tax_full[, "observation_name"] == "Cluster_334", "Species"] <- "Malacoplasma penetrans"
tax_full[tax_full[, "observation_name"] == "Cluster_334", "blast_perc_identity"] <- "99.2"

tax_full[tax_full[, "observation_name"] == "Cluster_1370", "Species"] <- "Borrelia afzelii"
tax_full[tax_full[, "observation_name"] == "Cluster_1370", "blast_perc_identity"] <- "99.6"

tax_full[tax_full[, "observation_name"] == "Cluster_5922", "Species"] <- "Eperythrozoon coccoides"

tax_full[tax_full[, "observation_name"] == "Cluster_6232", "Species"] <- "Eperythrozoon coccoides"


# Identify best affiliation level again

# Keep only ranks 
tax <- tax_full[, ranks]

# Reverse matrix to go from Species to Domain
tax_rev <- tax[, rev(seq_len(ncol(tax))), drop = FALSE]
ranks_rev <- rev(ranks)

# Find first non-NA per row
first_valid_index <- max.col(!is.na(tax_rev), ties.method = "first")

# In case some rows contain only NA 
all_na <- rowSums(!is.na(tax_rev)) == 0

# Extract best taxa and level
best_taxa <- tax_rev[cbind(seq_len(nrow(tax_rev)), first_valid_index)]
best_taxa_level <- ranks_rev[first_valid_index]

best_taxa[all_na] <- NA
best_taxa_level[all_na] <- NA

# Replace original rank columns with cleaned ones
tax_full[, ranks] <- tax

# REMOVE old best columns BEFORE cbind
tax_full <- tax_full[, !colnames(tax_full) %in% c("best_taxa", "best_taxa_level"), drop = FALSE]

# Insert new columns after last rank column
last_rank_pos <- max(match(ranks, colnames(tax_full)))

tax_final <- cbind(
  tax_full[, 1:last_rank_pos, drop = FALSE],
  best_taxa = best_taxa,
  best_taxa_level = best_taxa_level,
  tax_full[, (last_rank_pos + 1):ncol(tax_full), drop = FALSE]
)

#Drop rdp_taxo
tax_final <- tax_final[, !colnames(tax_final) %in% c("rdp_tax_and_bootstrap"), drop = FALSE]

# Back to phyloseq
tax_table(physeq_s16) <- tax_table(as.matrix(tax_final))

# 
# physeq_s16@tax_table %>% as.data.frame() %>%
# data.table::fwrite(file = "16s_tax_check.csv", sep = ";")
# 
# physeq_s16@otu_table %>% as.data.frame() %>%
#   data.table::fwrite(file = "16s_abund_check.csv", sep = ";")
# 



## Agglomerate by best identified taxa ----
physeq_s16 <- phyloseq::tax_glom(
  physeq = physeq_s16,
  taxrank = "best_taxa",
  bad_empty = c("Multi-affiliation", NA, "NA", "unknown species", "", " ", "\t"),
  NArm = TRUE)


physeq_s16 <- phyloseq::prune_samples(
  x = physeq_s16,
  samples = sample_sums(physeq_s16) != 0)
tracking_s16 <- bind_rows(
  tracking_s16, track_phyloseq(physeq_s16, step_name = "Agglomerate_cluster_best_level", step_number = 8))


#Remove tax columns not useful anymore
tax_full <- as(tax_table(physeq_s16), "matrix") 
tax <- tax_full[, c(ranks, "best_taxa")]
tax_table(physeq_s16) <- tax_table(as.matrix(tax))



## Visualization of filters ----
tracking_all <- tidyr::pivot_longer(data = tracking_s16, 
                                    cols = c("n_taxa", "n_samples", "total_reads"),
                                    names_to = "parameter",
                                    values_to = "parameter_value")
#Plot filtering steps
tracking_all %>% 
  arrange(step_number) %>%
  mutate(step_name = forcats::fct_inorder(step_name)) %>% 
  ggplot(aes(x = step_name, y = parameter_value, group = 1)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 2) +
  facet_grid(parameter ~ ., scales = "free_y") +
  #  expand_limits(y = 0) +
  labs(
    x = "Filtering step",
    y = "Value",
    title = "Tracking of filtering steps"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    strip.placement = "outside"
  )


# Convert to Presence / Absence data ----
physeq_s16 <- phyloseq::transform_sample_counts(physeq_s16, function(x) ifelse(x > 0, 1, 0))



# Save file  ----

saveRDS(object = physeq_s16, file =  here::here("data", "derived-data", "pathos_cleaning", "phyloseq_s16.rds"))

