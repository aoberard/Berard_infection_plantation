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



# Import data ----

#Cleaned Metabarcoding phyloseq objects
s16_phylo <- readRDS(file =  here::here("data", "derived-data", "pathos_cleaning", "phyloseq_s16.rds"))
rpob_glta_phylo <- readRDS(file =  here::here("data", "derived-data", "pathos_cleaning", "phyloseq_rpob_glta.rds"))

#Leptospira pcr data
lepto <- data.table::fread(file = here::here("data", "derived-data", "pathos_cleaning", "leptospira_pcr.csv"))

#Samples data
individuals <- read.delim(file = here::here("data", "raw-data", "20260320_BePrep_samples.csv"), sep = ";")


# Metabarcoding objects join ----

## Generate markers tables ----

# #Rename taxa in accordance to best_taxa
phyloseq::taxa_names(s16_phylo) <- as.data.frame(tax_table(s16_phylo)) %>% pull(best_taxa)

#Pivot otu_tables
s16 <- phyloseq::otu_table(s16_phylo) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(numero_centre = rownames(.))

barto <- phyloseq::otu_table(rpob_glta_phylo) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(numero_centre = rownames(.))


## Join to samples data ----
#Starting with 16S data
infectious <- individuals %>% 
  select(numero_centre) %>% 
  left_join(s16, by = c("numero_centre"))

#Adding rpoB and gltA data
infectious <- infectious %>%
  left_join(barto, by = c("numero_centre"))

#Replace NA by 0 in new data
infectious <- infectious %>%
  mutate(across(-numero_centre, ~ if_else(is.na(.), 0, .)))

#Transform names for easier manipulation
infectious <- infectious %>%
  rename_with(~ gsub(" ", "_", .))

#See benefit of Neoehrlichia mikurensis columns combining
infectious %>%
  filter(Neoehrlichia_mikurensis.y > 0) %>%
  filter(Neoehrlichia_mikurensis.x == 0) %>%
  select(numero_centre) %>%
  pull()

infectious <- infectious %>%
  mutate(
    Neoehrlichia_mikurensis = case_when(
      Neoehrlichia_mikurensis.x == 0 & Neoehrlichia_mikurensis.y > 0 ~ Neoehrlichia_mikurensis.y,
      TRUE ~ Neoehrlichia_mikurensis.x
    )
  ) %>%
  select(-Neoehrlichia_mikurensis.x, -Neoehrlichia_mikurensis.y)


#Identify if some samples were positives for Bartonella in 16S but are not with rpoB-gltA
infectious %>%
  filter(
    rowSums(
      across(
        phyloseq::taxa_names(rpob_glta_phylo) %>%
          gsub(" ", "_", .) %>%
          setdiff("Neoehrlichia_mikurensis")
      )
    ) == 0
  ) %>%
  filter(Bartonella > 0) %>%
  select(numero_centre) %>%
  pull()

#And the opposite
infectious %>%
  filter(
    rowSums(
      across(
        phyloseq::taxa_names(rpob_glta_phylo) %>%
          gsub(" ", "_", .) %>%
          setdiff("Neoehrlichia_mikurensis")
      )
    ) > 0
  ) %>%
  filter(Bartonella == 0) %>%
  select(numero_centre) %>%
  pull()

#Drop genus level Bartonella taxa
infectious <- infectious %>%
  select(-Bartonella)


#Specific diagnostic detection  join ----

## Leptospira PCR ----
infectious <- infectious %>%
  left_join(lepto)

#Are there NA for this specific diagnostic ?
infectious %>%
  filter(
    if_any(
      setdiff(names(lepto), "numero_centre"),
      is.na
    )
  ) %>%
  pull(numero_centre)

#Make choice regarding NA samples for the specific taxa (keep as NA or else ?)
#Here we make them 0 (because only two samples with low global prevalence)
infectious <- infectious %>%
  mutate(
    across(
      all_of(setdiff(names(lepto), "numero_centre")),
      ~ tidyr::replace_na(.x, 0)
      )
    )


# EXPLORATION TEMPO MYCOPLASMA

# infectious <- infectious %>% 
#   mutate(across(-numero_centre, ~ if_else(is.na(.), 0, .)))
# 
# infectious <- infectious %>% 
#   left_join(individuals)



#Goal a la fin tableau clean qui peut par la suite etre join a dataframe line var (pas travailler avec le foutoir ici - dans analyses ok)

# Attention à garder NA si pas regarder pour X échantillons ou pas dispo ou je ne sais quoi
