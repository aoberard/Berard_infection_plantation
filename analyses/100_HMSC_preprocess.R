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
names_inf10 <- readRDS(file = here::here("data", "derived-data", "asylvaticus_data", "names_infectious_apo10.rds"))


# Parameter choice  ----

# Choose which infectious agent to consider ▲
namesagent <- names_inf10

# Choose which explanatory/random variables to consider ▲
namesvar <- c("poids", "sexe", "maturity", "season", "category", "nb_apo")
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
  mutate(across(all_of(namesagent), as.numeric))


## Ensure wide format ----
#Data should be in wide format for hmsc
head(data_hmsc %>%  select(all_of(namesagent))) #fine here


## Remove rows with NA ----
#We remove rows with NA for either response or explanatory/random variables

#response
data_hmsc %>%
  filter(if_any(any_of(namesagent), is.na)) %>% 
  select(all_of( c("numero_centre", namesagent)))

data_hmsc <- data_hmsc %>%
  filter(! if_any(any_of(namesagent), is.na))


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
  select(any_of(c(namesagent, namesvar, namesrandom))) %>%
  str()


# Matrix setup ----

## Y matrix ----
Y <- data_hmsc %>% 
  tibble::column_to_rownames("numero_centre") %>% 
  select(all_of(namesagent))
  

head(Y)

prev.out <- (colSums(Y)/nrow(Y))*100
prev.out %>% sort()


## X matrix ----
X <- data_hmsc %>% 
  tibble::column_to_rownames("numero_centre") %>% 
  select(all_of(namesvar))

## T matrix ----


## Random effect (study design) data frame ----
studyDesign <- data_hmsc %>% 
  select(all_of(namesrandom)) %>% 
  rename( individual = numero_centre,
          site = numero_ligne,
          mission = field_session)

rownames(studyDesign) <- studyDesign$individual 

#of possible use if want to add spatial random effect
# data_hmsc %>% 
#   select(numero_centre, longitude_piege, latitude_piege) %>%  
#   sf::st_as_sf(
#     coords = c("longitude_piege", "latitude_piege"),
#     crs = "EPSG:2154",
#     remove = FALSE
#     )


# Save objects ----
dirpath <- here::here("data", "derived-data", "hmsc_out")

#Generate repository for future saved objects
if (!dir.exists(dirpath)) {
  dir.create(dirpath, recursive = TRUE)
}

save(Y,
     X,
     studyDesign,
     data_hmsc,
     file = here::here(dirpath, "allData.RData"))
