# Script parameters ----

# Remove saved elements (clean workspace)
rm(list = ls())
gc()

# Load libraries
library("dplyr")


# Import data ----

#Apodemus sylvaticus infectious data
inf <- data.table::fread(file = here::here("data", "derived-data", "asylvaticus_data", "asylvaticus_infectious.csv"))

names_inf <- readRDS(file = here::here("data", "derived-data", "asylvaticus_data", "names_infectious_apo.rds"))

#Infectious agents traits
traits <- readxl::read_excel(path = here::here("data", "raw-data", "20260520_agents_traits.xlsx"))



# Parameter choice  ----

# Choose which explanatory/random variables to consider ▲
namesvar <- c("poids", "sexe", "maturity", "season", "category", "nb_apo", "herb_cover", "shrub_cover", "tree_cover", "gauss500brd")
namesrandom <- c( "year", "field_session", "numero_ligne", "numero_centre")


# Data modification ----

## Ensure relevant columns are factors or numeric----

##Convert all character columns to factors
data_hmsc <- inf %>%
  mutate(across(where(is.character), as.factor))


##Convert numeric columns that should be treated as factors 
#numero_ligne
data_hmsc <- data_hmsc %>%
  mutate(numero_ligne = factor(sprintf("%02d", numero_ligne)))

#Every random effects)
data_hmsc <- data_hmsc %>%
  mutate(across(all_of(namesrandom), as.factor))


##Ensure response variables are numeric
data_hmsc <- data_hmsc %>%
  mutate(across(all_of(names_inf), as.numeric))


## Ensure wide format ----
#Data should be in wide format for hmsc
head(data_hmsc %>%  select(all_of(names_inf))) #fine here


## Remove rows with NA ----
#We remove rows with NA for either response or explanatory/random variables

#response
data_hmsc %>%
  filter(if_any(any_of(names_inf), is.na)) %>% 
  select(all_of( c("numero_centre", names_inf)))

data_hmsc <- data_hmsc %>%
  filter(! if_any(any_of(names_inf), is.na))


#explanatory
data_hmsc %>%
  filter(if_any(any_of(namesvar), is.na)) %>% 
  select(all_of( c("numero_centre", namesvar)))

data_hmsc <- data_hmsc %>%
  filter(! if_any(any_of(namesvar), is.na))


#random
data_hmsc %>%
  filter(if_any(any_of(namesrandom), is.na)) %>% 
  select(all_of( c("numero_centre", namesrandom)))

data_hmsc <- data_hmsc %>%
  filter(! if_any(any_of(namesrandom), is.na))


##Ensure no residual factor levels remain (as this would generate error during model runs)
data_hmsc <- data_hmsc %>%
  mutate(across(where(is.factor), droplevels))


## Verify the structure of data
data_hmsc %>%
  select(any_of(c(names_inf, namesvar, namesrandom))) %>%
  str()


## Traits dataframe modification ----
traits <- traits %>% 
  mutate(across(all_of(colnames(traits)), as.factor))



# Matrix setup ----

## Y table ----
Y <- data_hmsc %>% 
  tibble::column_to_rownames("numero_centre") %>% 
  select(all_of(names_inf))
  

head(Y)

prev.out <- (colSums(Y)/nrow(Y))*100
prev.out %>% sort()


## X table ----
X <- data_hmsc %>% 
  tibble::column_to_rownames("numero_centre") %>% 
  select(all_of(namesvar))


## T table ----
Tm <- traits %>% 
  as.data.frame()

rownames(Tm) <- Tm$names


## Random effect (study design) data frame ----
studyDesign <- data_hmsc %>% 
  select(all_of(namesrandom)) %>% 
  rename( individual = numero_centre,
          location = numero_ligne,
          mission = field_session,
          year = year)

rownames(studyDesign) <- studyDesign$individual 



## Random effect (spatial) matrix ----

#Defined at the studyDesign chosen "location" scale
locations.xy <- data_hmsc %>%
  select(numero_ligne, longitude_lambert93, latitude_lambert93) %>%
  distinct() %>% 
  tibble::column_to_rownames(var = "numero_ligne") %>% 
  as.matrix()


# Save objects ----
dirpath <- here::here("data", "derived-data", "hmsc_out")

#Generate repository for future saved objects
if (!dir.exists(dirpath)) {
  dir.create(dirpath, recursive = TRUE)
}

save(Y,
     X,
     Tm,
     studyDesign,
     locations.xy,
     data_hmsc,
     file = here::here(dirpath, "allData.RData"))
