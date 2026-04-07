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


#Graphical parameters 
category_palette <- c("pine_edge"          = "#3f826d",
                      "hedgerows"          = "#f2d0a4",
                      "broadleaved_forest" = "#a4a8d1")

common_rodent_palette <- c(
  "Clethrionomys glareolus" = "#748CAB",  
  "Apodemus sylvaticus"      = "#D0A98F"
)

category_order <- c("pine_edge", "hedgerows", "broadleaved_forest")
mission_order <- c("Juin_2023", "Octobre_2023", "Juin_2024", "Septembre_2024", "Juin_2025")
year_order <- c("2023", "2024", "2025")


# Import data ----

#Samples data
dissected <- read.delim(file = here::here("data", "raw-data", "dissected_indiv.csv"), sep = ",")

#Site environmental data
envi <- data.table::fread(file = here::here("data", "raw-data", "lines_variability.csv"))

#Cleaned infectious agents data
infecti <- data.table::fread(file = here::here("data", "derived-data", "infectious_cleaning", "infectious_data.csv"))


# Generate data ----

## Join dataframes

#Infectious agents and individuals
infected <- infecti %>% 
  left_join(dissected)

#Adding environment
infected <- infected %>% 
  left_join(envi, by = c("numero_ligne" = "num",
                         "code_mission" = "code_mission"))


# Exploration ----

#Infectious agents name
names_infectious <- setdiff(colnames(infecti), y = "numero_centre")


## Across rodent species ----

#Number of individuals per rodent species
infected %>%
  group_by(taxon_mamm) %>%
  count()

#Number infected per taxa
infected %>%
  mutate(across(all_of(names_infectious), ~ replace(., . > 0, 1))) %>%
  group_by(taxon_mamm) %>%
  summarise(across(all_of(names_infectious), sum))

#List of pathogen per species
list_pathos_per_species <- list()
for (i in unique(infected$taxon_mamm)) {
  list_pathos_per_species[[i]] <- infected %>%
    filter(taxon_mamm == i) %>%
    select(all_of(names_infectious)) %>%
    select(where(~ any(. > 0))) %>%
    colnames()
}
list_pathos_per_species


# Compute prevalence for every species

#observed number of cases
prev_df <- infected %>%
  tidyr::pivot_longer(
    cols = all_of(names_infectious),
    names_to = "pathogen",
    values_to = "infected"
  ) %>%
  group_by(taxon_mamm, pathogen) %>%
  summarise(
    n = n(),
    cases = sum(infected),
    .groups = "drop"
  )

#prevalence IC
prev_df <- prev_df %>%
  rowwise() %>%
  mutate(
    ci = list(prevalence::propCI(cases, n, method = "exact")),
    prevalence = ci[, "p"],
    lower = ci[, "lower"],
    upper = ci[, "upper"]
  ) %>%
  ungroup() %>%
  select(-ci)

#remove useless rodent species
prev_df <- prev_df %>% 
  filter(!taxon_mamm %in% c("Mus musculus", "Micromys minutus"))

#Plot graph
prev_df %>%
  mutate(
    pathogen = as.factor(reorder(pathogen, cases)),
    taxon_mamm = as.factor(taxon_mamm)
  ) %>%
  ggplot(aes(x = prevalence, y = pathogen, color = taxon_mamm)) +
  geom_point(
    position = position_dodge(width = 0.5),
    size = 2
  ) +
  geom_errorbarh(
    aes(xmin = lower, xmax = upper),
    position = position_dodge(width = 0.5),
    width = 0.3
  ) +
  scale_color_manual(values = common_rodent_palette) +
  theme_bw()


## For A. sylvaticus ----

#Vector of pathogens present in A. sylvaticus
names_infectious_apo


#Vector of pathogens present in A. sylvaticus with at least 10% global prevalence 
names_infectious_apo






#Prevalence variable
#selon ?

#nbr d'individus infectés variable (on regarde que la prevalence a chaque fois c'est bon wsh, c'est pas ça le "danger")
# mais attention implique de standardiser data pour meme effort capture : + probleme c'est que c'est pas des vrais
#abondances qu'on regarde mais du relatif - pas forcement grave là aussi on peut dire c'est relatif non ?
#pour faire ça regarder les individus supprimés dans autre script quand retrait 4 eme nuit de piegeage
#peut-être autre truc a check


#réseaux co-occurences 

#coinfection representation ?



# specificité hote 
#table, reseau, heatmap


#pour l'instant y'a pas : les variables des individus potentiellment a integrer : maturité sex, classe d'age, condition physique

#pour les variables environnements : jamais regardé hobo


