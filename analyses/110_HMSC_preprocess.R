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
traits <- readxl::read_excel(path = here::here("data", "raw-data", "agents_traits.xlsx"))

#Apodemus sylvaticus number per site-session (over the three first trapping nights)
effectif <- data.table::fread(file = here::here("data", "raw-data", "apodemus_number.csv"))

#Site environmental data
envi <- data.table::fread(file = here::here("data", "raw-data", "lines_variability.csv"))


# Parameter choice  ----

# Choose which explanatory/random variables to consider ▲
namesvar <- c("season", "category")
namesrandom <- c( "year", "field_session", "numero_ligne", "site_session")


# Data modification ----

## Keep individuals trapped within 3 nights ----
inf_3night <- inf %>% 
  filter(after_3days == FALSE)


## Summarise data per site-session ----

#Regroup infected individuals data
infect_counts <- inf_3night %>%
  group_by(numero_ligne, field_session) %>% 
  summarise(across(all_of(names_inf),
                   ~sum(.x, na.rm = TRUE)),
           .groups = "drop")

#Join to effectif data (otherwise 0 abundance sites would be excluded)
data_hmsc <- effectif %>%
  left_join(infect_counts,
            by = c("numero_ligne", "field_session"))

#Join back environmental and sampling session information
data_hmsc <- data_hmsc %>% 
  left_join(envi %>%
              select(
                num,
                category,
                longitude_lambert93,
                latitude_lambert93
              ) %>%
              distinct(),
            by = c("numero_ligne" = "num"))

data_hmsc <- data_hmsc %>%
  left_join(inf_3night %>%
              select(
                field_session,
                season,
                year,
              ) %>%
              distinct(),
            by = join_by("field_session")) %>%
  mutate(site_session = interaction(numero_ligne, field_session))

#Ensure 0 abundance site are not NA but 0 infected individuals
data_hmsc <- data_hmsc %>%
  mutate(
    across(
      all_of(names_inf),
      ~tidyr::replace_na(.x, 0)
    )
  )


## Remove site-session with too many individuals not dissected

#Site-session with individuals not dissected
data_hmsc %>% 
  filter(total > dissected) %>% 
  nrow()
  
#Remove site-session with more than 1 individual not dissected
data_hmsc <- data_hmsc %>%
  filter(dissected >= total - 1)


## Ensure relevant columns are factors or numeric----

##Convert all character columns to factors
data_hmsc <- data_hmsc %>%
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
  select(all_of( c("site_session", names_inf)))

data_hmsc <- data_hmsc %>%
  filter(! if_any(any_of(names_inf), is.na))


#explanatory
data_hmsc %>%
  filter(if_any(any_of(namesvar), is.na)) %>% 
  select(all_of( c("site_session", namesvar)))

data_hmsc <- data_hmsc %>%
  filter(! if_any(any_of(namesvar), is.na))


#random
data_hmsc %>%
  filter(if_any(any_of(namesrandom), is.na)) %>% 
  select(all_of( c("site_session", namesrandom)))

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
  tibble::column_to_rownames("site_session") %>% 
  select(all_of(names_inf))


head(Y)

mean.out <- (colSums(Y)/nrow(Y))
mean.out %>% sort()


## X table ----
X <- data_hmsc %>% 
  tibble::column_to_rownames("site_session") %>% 
  select(all_of(namesvar))


## T table ----
Tm <- traits %>% 
  as.data.frame()

rownames(Tm) <- Tm$names


## Random effect (study design) data frame ----
studyDesign <- data_hmsc %>% 
  select(all_of(namesrandom)) %>% 
  rename( site_session = site_session,
          site = numero_ligne,
          mission = field_session,
          year = year) %>% 
  mutate(location = site)%>% 
  as.data.frame()

rownames(studyDesign) <- studyDesign$site_session 



## Random effect (spatial) matrix ----

#Defined at the studyDesign chosen "site" scale
locations.xy <- data_hmsc %>%
  select(numero_ligne, longitude_lambert93, latitude_lambert93) %>%
  distinct() %>% 
  tibble::column_to_rownames(var = "numero_ligne") %>% 
  as.matrix()


# Save objects ----
dirpath <- here::here("data", "derived-data", "hmsc_ab_out")

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

