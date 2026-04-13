
# Script parameters ----

#Remove saved elements
rm(list = ls())
gc()

# Library
library("dplyr")
library("ggplot2")


# Import data ----
lepto_pcr <- readxl::read_excel(path = here::here("data", "raw-data", "diagnostic", "20260330_data_qPCR_lepto_BePrep.xlsx"))


# Generate data ----
lepto_pcr <- lepto_pcr %>%
  mutate(Leptospira = 
           case_when(
             !is.na(Lfb1_PCR) & Lfb1_PCR == 1 ~ 1,
             !is.na(Lfb1_PCR) & Lfb1_PCR == 0 ~ 0,
             is.na(Lfb1_PCR) & lipL32_result == 1 ~ 1,
             is.na(Lfb1_PCR) & lipL32_result == 0 ~ 0,
             .default  = NA
           ))


# Save data ----
dir.create(
  here::here("data", "derived-data", "infectious_cleaning"),
  recursive = TRUE,
  showWarnings = FALSE
)

lepto_pcr %>% 
  rename(numero_centre = ID_rodent) %>% 
  select(numero_centre, Leptospira) %>% 
  data.table::fwrite(file = here::here("data", "derived-data", "infectious_cleaning", "leptospira_pcr.csv"))
