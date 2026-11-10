#' Berard_infection_script
#' 
#' @author Aloïs Berard \email{aloisberard@outlook.fr}
#' 
#' @date 2026/11/09



## Install Dependencies ----
if (!("renv" %in% installed.packages())) {
  install.packages("renv")
}

renv::restore()




## Run Project ----

### Generate data ----

#Bartonella data
source(here::here("analyses", "000_rpoB-gltA_cleaning.R"))
#16S data
source(here::here("analyses", "001_16S_spleen_cleaning.R"))
#Leptospira data
source(here::here("analyses", "002_lepto_data_cleaning.R"))
#Combine datasets
source(here::here("analyses", "003_combine_infectious_agents.R"))
#Filter dataset (keeping only A. sylvaticus)
source(here::here("analyses", "010_data_join_and_exploration.R"))


### Analysis ----

#### Diversity ----

#Richness linear modelling
source(here::here("analyses", "011_alpha_div.R"))


#### HMSC infection probability ----
#HMSC was done on genotoul cluster as it may take a long time with current parameters (>20 days)

#Prepare data
source(here::here("analyses", "100_HMSC_preprocess.R"))
#Define HMSC models
source(here::here("analyses", "101_HMSC_define_models.R"))
#Fitting models (time consuming)
source(here::here("analyses", "102_HMSC_fit_models.R"))
#Evaluating convergence
source(here::here("analyses", "103_HMSC_evaluate_convergence.R"))
#Computing fit (time consuming)
source(here::here("analyses", "104_HMSC_compute_model_fit.R"))
#Checking fit outputs for model selection
source(here::here("analyses", "105_HMSC_show_model_fit.R"))
#Checking parameter outputs
source(here::here("analyses", "106_HMSC_show_parameter_estimates.R"))
#Making and checking predictions
source(here::here("analyses", "107_HMSC_predictions.R"))


#### HMSC infected host abundance ----
#HMSC was done on genotoul cluster as it may take a long time with current parameters (>20 days)

#Prepare data
source(here::here("analyses", "110_HMSC_preprocess.R"))
#Define HMSC models
source(here::here("analyses", "111_HMSC_define_models.R"))
#Fitting models (time consuming)
source(here::here("analyses", "112_HMSC_fit_models.R"))
#Evaluating convergence
source(here::here("analyses", "113_HMSC_evaluate_convergence.R"))
#Computing fit (time consuming)
source(here::here("analyses", "114_HMSC_compute_model_fit.R"))
#Checking fit outputs for model selection
source(here::here("analyses", "115_HMSC_show_model_fit.R"))
#Checking parameter outputs
source(here::here("analyses", "116_HMSC_show_parameter_estimates.R"))
#Making and checking predictions
source(here::here("analyses", "117_HMSC_predictions.R"))


#### Combine and compare HMSC results  ----
source(here::here("analyses", "200_HMSC_compare_approaches.R"))

