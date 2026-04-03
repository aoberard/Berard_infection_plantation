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

#Cleaned infectious agents data
infecti <- data.table::fread(file = here::here("data", "derived-data", "infectious_cleaning", "infectious_data.csv"))

#Samples data
individuals <- read.delim(file = here::here("data", "raw-data", "20260320_BePrep_samples.csv"), sep = ";")


# Generate data ----

## Join dataframe
infected <- infecti %>% 
  left_join(individuals)


# Exploration 



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




