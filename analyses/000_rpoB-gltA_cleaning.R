
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

ranks <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")



# Import data ----

# Import Metabarcoding data
rpob_abundance <- data.table::fread(
  file = here::here("data", "raw-data", "diagnostic", "2026.03.20-rpoB.abundance.final.txt"))
rpob_taxo <- data.table::fread(
  file = here::here("data", "raw-data", "diagnostic", "2023.03.23-rpoB.taxonomy.final.blastn.clean.txt"))

glta_abundance <-  data.table::fread(
  file = here::here("data", "raw-data", "diagnostic", "2026.03.20-gltA.abundance.final.txt"))
glta_taxo <-  data.table::fread(
  file = here::here("data", "raw-data", "diagnostic", "2026.03.10-gltA.taxonomy.final.blastn.txt"))

#Ensure rowname for metabarcoding data
rownames(rpob_abundance) <- rpob_abundance$V1
rpob_abundance$V1 <- NULL
rownames(rpob_taxo) <- rpob_taxo$V1
rpob_taxo$V1 <- NULL

rownames(glta_abundance) <- glta_abundance$V1
glta_abundance$V1 <- NULL
rownames(glta_taxo) <- glta_taxo$V1
glta_taxo$V1 <- NULL


#Samples data
individuals <- read.delim(file = here::here("data", "raw-data", "20260320_BePrep_samples.csv"), sep = ";")


# Data object creation  ----

## Keep relevant samples  ----

# Check for presence of samples present in multiple runs
rpob_abundance %>% select(stringr::str_subset(string = colnames(rpob_abundance), pattern = "_Run[0-9]+")) %>% colnames()
glta_abundance %>% select(stringr::str_subset(string = colnames(glta_abundance), pattern = "_Run[0-9]+")) %>% colnames()

# Keep duplicated samples from prefered Run
rpob_abundance <- rpob_abundance %>% 
  dplyr::select(-dplyr::matches("_Run01$|_Run05$"))

glta_abundance <- glta_abundance %>% 
  dplyr::select(-dplyr::matches("_Run01$|_Run05$"))


# Rename duplicated samples name (remove Run mention in name)
rpob_abundance <- rpob_abundance %>% 
  dplyr::rename_with(~ stringr::str_replace(.x, "_Run[0-9]+", ""))

glta_abundance <- glta_abundance %>% 
  dplyr::rename_with(~ stringr::str_replace(.x, "_Run[0-9]+", ""))


# Keep relevant sample type (Unsure we have only .SP (spleen) samples)
rpob_abundance %>% select(!stringr::str_subset(string = colnames(rpob_abundance), pattern = ".SP$")) %>% colnames()
rpob_abundance <- rpob_abundance %>% 
  select(stringr::str_subset(string = colnames(rpob_abundance), pattern = ".SP$"))

glta_abundance %>% select(!stringr::str_subset(string = colnames(glta_abundance), pattern = ".SP$")) %>% colnames()
glta_abundance <- glta_abundance %>% 
  select(stringr::str_subset(string = colnames(glta_abundance), pattern = ".SP$"))


# Make names concordant with individuals table
rpob_abundance <- rpob_abundance %>%
  rename_with(~ stringr::str_remove(., "\\.SP$"))

glta_abundance <- glta_abundance %>%
  rename_with(~ stringr::str_remove(., "\\.SP$"))


# Keep only BEPREP project samples
rpob_abundance <- rpob_abundance %>%
  select(dplyr::any_of(individuals %>% dplyr::pull(numero_centre)))

glta_abundance <- glta_abundance %>%
  select(dplyr::any_of(individuals %>% dplyr::pull(numero_centre)))


#Check that both marker contain the same samples
if (!all(colnames(rpob_abundance) %in% colnames(glta_abundance))) {
  warning("Some rpob samples are not found in glta data:\n",
          paste(setdiff(colnames(rpob_abundance), colnames(glta_abundance)), collapse = ", "))
}

if (!all(colnames(glta_abundance) %in% colnames(rpob_abundance))) {
  warning("Some glta samples are not found in rpob data:\n",
          paste(setdiff(colnames(glta_abundance), colnames(rpob_abundance)), collapse = ", "))
}

#Keep only individuals in the rpob and glta abundance tables
individuals <- individuals %>%
  filter(numero_centre %in% colnames(rpob_abundance)
         & numero_centre %in% colnames(glta_abundance))


## Generate phyloseq objects  ----

#Convert abundance dataframe to matrix
matrix_ab_rpob <- rpob_abundance %>% as.matrix()
rownames(matrix_ab_rpob) <- rownames(rpob_abundance)

matrix_ab_glta <- glta_abundance %>% as.matrix()
rownames(matrix_ab_glta) <- rownames(glta_abundance)


#Convert taxonomy dataframe to matrix
matrix_tax_rpob <- rpob_taxo %>%
  tidyr::separate_wider_delim(cols = "blast_taxonomy", delim = ";", names = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
  select(Domain, Phylum, Class, Order, Family, Genus, Species, observation_name,  blast_perc_query_coverage, blast_perc_identity, seed_sequence, observation_sum) %>%
  as.matrix()
rownames(matrix_tax_rpob) <- rownames(rpob_taxo)


matrix_tax_glta <- glta_taxo %>%
  tidyr::separate_wider_delim(cols = "blast_taxonomy", delim = ";", names = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
  select(Domain, Phylum, Class, Order, Family, Genus, Species, observation_name,  blast_perc_query_coverage, blast_perc_identity, seed_sequence, observation_sum) %>%
  as.matrix()
rownames(matrix_tax_glta) <- rownames(glta_taxo)


#Prepare metadata (samples)
#Make numero_centre as rownames for individuals dataframe
individuals <- individuals %>%
  tibble::column_to_rownames(var = "numero_centre")

individuals$numero_centre <- rownames(individuals)


#Convert matrix and metadata to phyloseq objects
OTU_rpob <- phyloseq::otu_table(matrix_ab_rpob, taxa_are_rows = TRUE)
TAX_rpob <- phyloseq::tax_table(matrix_tax_rpob )

OTU_glta <- phyloseq::otu_table(matrix_ab_glta, taxa_are_rows = TRUE)
TAX_glta <- phyloseq::tax_table(matrix_tax_glta )

SAMPL <- phyloseq::sample_data(individuals)


#Generate phyloseq objects
physeq_rpob <- phyloseq::phyloseq(OTU_rpob, TAX_rpob, SAMPL)

physeq_glta <- phyloseq::phyloseq(OTU_glta, TAX_glta, SAMPL)


# Filtering ----
#After every filtering steps, samples not containing reads are removed in order
#to follow filtering effect on sample number

#Initialize tracking table
physeq_rpob <- phyloseq::prune_samples(
  x = physeq_rpob,
  samples = sample_sums(physeq_rpob) != 0)
physeq_glta <- phyloseq::prune_samples(
  x = physeq_glta,
  samples = sample_sums(physeq_glta) != 0)
tracking_rpob <- track_phyloseq(physeq_rpob, step_name = "Beprep_samples_with_reads", step_number = 0)

tracking_glta <- track_phyloseq(physeq_glta, step_name = "Beprep_samples_with_reads", step_number = 0)


## Remove clusters with 0 reads ----
physeq_rpob <- phyloseq::prune_taxa(
  taxa = phyloseq::taxa_sums(physeq_rpob) > 0,
  x =  physeq_rpob)

physeq_glta <- phyloseq::prune_taxa(
  taxa = phyloseq::taxa_sums(physeq_glta) > 0,
  x =  physeq_glta)

physeq_rpob <- phyloseq::prune_samples(
  x = physeq_rpob,
  samples = sample_sums(physeq_rpob) != 0)
physeq_glta <- phyloseq::prune_samples(
  x = physeq_glta,
  samples = sample_sums(physeq_glta) != 0)
tracking_rpob <- bind_rows(
  tracking_rpob, track_phyloseq(physeq_rpob, step_name = "Removed_0_reads_clusters", step_number = 1))
tracking_glta <- bind_rows(
  tracking_glta, track_phyloseq(physeq_glta, step_name = "Removed_0_reads_clusters", step_number = 1))


## Remove clusters with reads below 0.01 % of total ---- 
physeq_rpob <- phyloseq::prune_taxa(
  taxa = phyloseq::taxa_sums(physeq_rpob) >= (sum(phyloseq::taxa_sums(physeq_rpob)) * 0.0001),
  x =  physeq_rpob)

physeq_glta <- phyloseq::prune_taxa(
  taxa = phyloseq::taxa_sums(physeq_glta) >= (sum(phyloseq::taxa_sums(physeq_glta)) * 0.0001),
  x =  physeq_glta)

physeq_rpob <- phyloseq::prune_samples(
  x = physeq_rpob,
  samples = sample_sums(physeq_rpob) != 0)
physeq_glta <- phyloseq::prune_samples(
  x = physeq_glta,
  samples = sample_sums(physeq_glta) != 0)
tracking_rpob <- bind_rows(
  tracking_rpob, track_phyloseq(physeq_rpob, step_name = "Removed_below_0.01%_clusters", step_number = 2))
tracking_glta <- bind_rows(
  tracking_glta, track_phyloseq(physeq_glta, step_name = "Removed_below_0.01%_clusters", step_number = 2))


## Remove clusters without affiliation ---- 
physeq_rpob <- phyloseq::prune_taxa( 
  taxa = as.vector(phyloseq::tax_table(physeq_rpob)[, "Domain"]) != "no",
  x = physeq_rpob)

physeq_glta<- phyloseq::prune_taxa( 
  taxa = as.vector(phyloseq::tax_table(physeq_glta)[, "Domain"]) != "no",
  x = physeq_glta)

physeq_rpob <- phyloseq::prune_samples(
  x = physeq_rpob,
  samples = sample_sums(physeq_rpob) != 0)
physeq_glta <- phyloseq::prune_samples(
  x = physeq_glta,
  samples = sample_sums(physeq_glta) != 0)
tracking_rpob <- bind_rows(
  tracking_rpob, track_phyloseq(physeq_rpob, step_name = "Removed_NA_clusters", step_number = 3))
tracking_glta <- bind_rows(
  tracking_glta, track_phyloseq(physeq_glta, step_name = "Removed_NA_clusters", step_number = 3))


## Remove clusters affiliated to Chordata ---- 
physeq_rpob <- phyloseq::prune_taxa( 
  taxa = as.vector(phyloseq::tax_table(physeq_rpob)[, "Domain"]) != "Chordata",
  x = physeq_rpob)

physeq_glta<- phyloseq::prune_taxa( 
  taxa = as.vector(phyloseq::tax_table(physeq_glta)[, "Domain"]) != "Chordata",
  x = physeq_glta)

physeq_rpob <- phyloseq::prune_samples(
  x = physeq_rpob,
  samples = sample_sums(physeq_rpob) != 0)
physeq_glta <- phyloseq::prune_samples(
  x = physeq_glta,
  samples = sample_sums(physeq_glta) != 0)
tracking_rpob <- bind_rows(
  tracking_rpob, track_phyloseq(physeq_rpob, step_name = "Removed_Chordata_clusters", step_number = 4))
tracking_glta <- bind_rows(
  tracking_glta, track_phyloseq(physeq_glta, step_name = "Removed_Chordata_clusters", step_number = 4))


## Remove clusters that are not Bartonella or Neoehrlichia ---- 
physeq_rpob <- phyloseq::prune_taxa( 
  taxa = as.vector(phyloseq::tax_table(physeq_rpob)[, "Genus"]) %in% c("Bartonella", "Neoehrlichia"),
  x = physeq_rpob)

physeq_glta <- phyloseq::prune_taxa( 
  taxa = as.vector(phyloseq::tax_table(physeq_glta)[, "Genus"]) %in% c("Bartonella", "Neoehrlichia"),
  x = physeq_glta)

physeq_rpob <- phyloseq::prune_samples(
  x = physeq_rpob,
  samples = sample_sums(physeq_rpob) != 0)
physeq_glta <- phyloseq::prune_samples(
  x = physeq_glta,
  samples = sample_sums(physeq_glta) != 0)
tracking_rpob <- bind_rows(
  tracking_rpob, track_phyloseq(physeq_rpob, step_name = "Removed_not_Bartonella_clusters", step_number = 5))
tracking_glta <- bind_rows(
  tracking_glta, track_phyloseq(physeq_glta, step_name = "Removed_not_Bartonella_clusters", step_number = 5))


## Remove clusters badly affiliated ---- 
coverage <- as.numeric(
  phyloseq::tax_table(physeq_rpob)[, "blast_perc_query_coverage"])
identity <- as.numeric(
  phyloseq::tax_table(physeq_rpob)[, "blast_perc_identity"])
physeq_rpob <- phyloseq::prune_taxa(
  coverage >= 98 & identity >= 95,
  physeq_rpob
)

coverage <- as.numeric(
  phyloseq::tax_table(physeq_glta)[, "blast_perc_query_coverage"])
identity <- as.numeric(
  phyloseq::tax_table(physeq_glta)[, "blast_perc_identity"])
physeq_glta <- phyloseq::prune_taxa(
  coverage >= 98 & identity >= 95,
  physeq_glta
)

physeq_rpob <- phyloseq::prune_samples(
  x = physeq_rpob,
  samples = sample_sums(physeq_rpob) != 0)
physeq_glta <- phyloseq::prune_samples(
  x = physeq_glta,
  samples = sample_sums(physeq_glta) != 0)
tracking_rpob <- bind_rows(
  tracking_rpob, track_phyloseq(physeq_rpob, step_name = "Remove_badly_affiliated_clusters", step_number = 6))
tracking_glta <- bind_rows(
  tracking_glta, track_phyloseq(physeq_glta, step_name = "Remove_badly_affiliated_clusters", step_number = 6))


## Merging clusters at species level ----
physeq_rpob <- phyloseq::tax_glom(
  physeq = physeq_rpob, 
  taxrank = "Species",
  bad_empty = c("Multi-affiliation", NA, "NA", "unknown species", "", " ", "\t"), 
  NArm = TRUE)

physeq_glta <- phyloseq::tax_glom(
  physeq = physeq_glta, 
  taxrank = "Species",
  bad_empty = c("Multi-affiliation", NA, "NA", "unknown species", "", " ", "\t"), 
  NArm = TRUE)

physeq_rpob <- phyloseq::prune_samples(
  x = physeq_rpob,
  samples = sample_sums(physeq_rpob) != 0)
physeq_glta <- phyloseq::prune_samples(
  x = physeq_glta,
  samples = sample_sums(physeq_glta) != 0)
tracking_rpob <- bind_rows(
  tracking_rpob, track_phyloseq(physeq_rpob, step_name = "After_merging_at_species_level", step_number = 7))
tracking_glta <- bind_rows(
  tracking_glta, track_phyloseq(physeq_glta, step_name = "After_merging_at_species_level", step_number = 7))


#Remove tax columns not useful anymore
tax_full <- as(tax_table(physeq_rpob), "matrix") 
tax <- tax_full[, ranks]
tax_table(physeq_rpob) <- tax_table(as.matrix(tax))

tax_full <- as(tax_table(physeq_glta), "matrix") 
tax <- tax_full[, ranks]
tax_table(physeq_glta) <- tax_table(as.matrix(tax))



## Visualization of filters ----

#Join tracking dataframes and make it tidy format
tracking_rpob$marker <- "rpoB"
tracking_glta$marker <- "gltA"
tracking_all <- bind_rows(tracking_rpob, tracking_glta)

tracking_all <- tidyr::pivot_longer(data = tracking_all, 
                                    cols = c("n_taxa", "n_samples", "total_reads"),
                                    names_to = "parameter",
                                    values_to = "parameter_value")

#Plot filtering steps
tracking_all %>% 
  arrange(step_number) %>%
  mutate(step_name = factor(step_name, levels = unique(step_name))) %>% 
  ggplot(aes(x = step_name, y = parameter_value, color = marker, group = marker)) +
  geom_line() +
  geom_point() +
  facet_grid(parameter ~ marker, scales = "free_y") +
  expand_limits(y = 0) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



# Merging markers ----

## Differences between markers ----

### Of remaining samples ---- 
samples_rpob <- phyloseq::sample_names(physeq_rpob)
samples_glta <- phyloseq::sample_names(physeq_glta)
samples_beprep <- individuals %>% dplyr::pull(numero_centre)

length(intersect(samples_rpob, samples_glta))
length(setdiff(samples_rpob, samples_glta))
length(setdiff(samples_glta, samples_rpob))
length(setdiff(samples_beprep, samples_rpob))
length(setdiff(samples_beprep, samples_glta))

euler_plot_samples <- eulerr::euler(list(
  BEPREP = samples_beprep,
  rpoB = samples_rpob,
  gltA = samples_glta
))
plot(euler_plot_samples,
     alpha = 1,
     labels = TRUE,
     quantities = TRUE,
     main = "Sample retention and overlap across BEPREP, rpoB and gltA")


### Of remaining taxa ----
taxa_names(physeq_rpob) <-  as.vector(phyloseq::tax_table(physeq_rpob)[, "Species"])
taxa_names(physeq_glta) <-  as.vector(phyloseq::tax_table(physeq_glta)[, "Species"])

taxa_rpob <- phyloseq::taxa_names(physeq_rpob)
taxa_glta <- phyloseq::taxa_names(physeq_glta)

data.frame(
  category = c("shared", "rpob_only", "glta_only"),
  n = c(
    length(intersect(taxa_rpob, taxa_glta)),
    length(setdiff(taxa_rpob, taxa_glta)),
    length(setdiff(taxa_glta, taxa_rpob))
  )
) %>%
  ggplot(aes(x = category, y = n)) +
  geom_col()


### Of occurences across shared samples and taxa ----

#Convert to 0/1
pa_rpob <- phyloseq::transform_sample_counts(physeq_rpob, function(x) ifelse(x > 0, 1, 0))
pa_glta <- phyloseq::transform_sample_counts(physeq_glta, function(x) ifelse(x > 0, 1, 0))

# keep only shared taxa and samples shared among markers
shared_samples <-intersect(phyloseq::sample_names(pa_rpob), phyloseq::sample_names(pa_glta))
pa_rpob <- phyloseq::prune_samples(shared_samples, pa_rpob)
pa_glta <- phyloseq::prune_samples(shared_samples, pa_glta)

shared_taxa <- intersect(phyloseq::taxa_names(pa_rpob), phyloseq::taxa_names(pa_glta))
pa_rpob <- phyloseq::prune_taxa(shared_taxa, pa_rpob)
pa_glta <- phyloseq::prune_taxa(shared_taxa, pa_glta)

#Extract abundance matrix and align rows
mat_rpob <- as(otu_table(pa_rpob), "matrix")
mat_glta <- as(otu_table(pa_glta), "matrix")

mat_rpob <- mat_rpob[shared_taxa, , drop = FALSE]
mat_glta <- mat_glta[shared_taxa, , drop = FALSE]

#Comparison matrix into dataframe(1 = rpob only, 2 = glta only, 3 = shared occurrences)
comparison <- mat_rpob + 2 * mat_glta

summary_indiv <- apply(comparison, 2, function(x) {
  c(
    rpob_only = sum(x == 1),
    glta_only = sum(x == 2),
    shared = sum(x == 3)
  )
})

summary_indiv <- as.data.frame(t(summary_indiv))
summary_indiv$sample <- rownames(summary_indiv)
summary_indiv$total <- summary_indiv$rpob_only +
  summary_indiv$glta_only +
  summary_indiv$shared

summary_indiv$sample <- factor(
  summary_indiv$sample,
  levels = summary_indiv$sample[order(summary_indiv$shared, decreasing = TRUE)])

summary_long <- summary_indiv %>%
  tidyr::pivot_longer(cols = c(rpob_only, glta_only, shared),
               names_to = "category",
               values_to = "n")

#Visualize comparison
ggplot(summary_long, aes(x = sample, y = n, fill = category)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


## Merging ----
#Rename taxa in accordance to species (for merging)
phyloseq::taxa_names(physeq_rpob) <- as.data.frame(tax_table(physeq_rpob)) %>% pull(Species)
phyloseq::taxa_names(physeq_glta) <- as.data.frame(tax_table(physeq_glta)) %>% pull(Species)

#Union merging
physeq_merged <- phyloseq::merge_phyloseq(physeq_rpob, physeq_glta)


## Convert to Presence / Absence data ----
physeq_merged <- phyloseq::transform_sample_counts(physeq_merged, function(x) ifelse(x > 0, 1, 0))


# Save file  ----
dir.create(
  here::here("data", "derived-data", "infectious_cleaning"),
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(object = physeq_merged, file =  here::here("data", "derived-data", "infectious_cleaning", "phyloseq_rpob_glta.rds"))
