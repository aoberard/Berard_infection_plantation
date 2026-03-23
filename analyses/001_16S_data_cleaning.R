
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

## List of rodents identified putative pathogens ----
#Beware, using this list should replace prior exploration of the cluster present in the metabarcoding dataset
putative_pathos_genra <- c("Bartonella", "Neoehrlichia", "Mycoplasma", "Eperythrozoon",
                           "Borrelia", "Borreliella", "Streptobacillus", "Ehrlichia", "Rodentibacter", 
                           "Spiroplasma", "Brevinema", "Orientia", "Rickettsia", "Leptospira",
                           "Yersinia", "Actinobacillus", "Treponema", "Chlamydia", "Neisseria",
                           "Pasteurella", "Francisella", "Brucella", "Coxiella", "Anaplasma")
putative_pathos_family <- c("Sarcocystidae")


# Import data ----

# Import mMetabarcoding data
s16_abundance <- data.table::fread(
  file = here::here("data", "raw-data", "diagnostic", "2026.03.04-16S.abundance.final.txt"))
s16_taxo <- data.table::fread(
  file = here::here("data", "raw-data", "diagnostic", "2026.03.04-16S.taxonomy.final.txt"))

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
  tidyr::separate_wider_delim(cols = "blast_taxonomy", delim = ";", names = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
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



