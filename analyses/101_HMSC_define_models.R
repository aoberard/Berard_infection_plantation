# Script parameters ----

#INPUT : Original datafiles of the case study, placed in the data folder.
#OUTPUT : Unfitted models, i.e., the list of Hmsc model(s) that have been defined
#but not fitted yet, stored in the file "models/unfitted_models.RData".

# Remove saved elements (clean workspace)
rm(list = ls())
gc()

# Load libraries
library("dplyr")
library("Hmsc")

# Make script reproducible
set.seed(1)

# Set directories 
dataDir <- here::here("data", "derived-data", "hmsc_out")
modelDir <- file.path(dataDir, "models")
if(!dir.exists(modelDir)) {
  dir.create(modelDir, recursive = TRUE)
}

# Import data ----
load(here::here(dataDir, "allData.RData"))


# Explore data ----

#Here done previously but always make sure no strong collinearity exists among explanatory variables used in models ! 
#Data exploration is highly advised also

# Just make sure matrix X, Y, and the studyDesign data frame have aligned row names
stopifnot(
  identical(rownames(X), rownames(Y)),
  identical(rownames(X), rownames(studyDesign))
)
message("Row names of X, Y and studyDesign are aligned.")


# Set up the models ----

## Models' fixed effects formula ----
XFormula.null = ~ 1
XFormula.host = ~ sexe + poids + maturity


## Random effects ----

#Define unstructured random effect at spatiotemporal unit
rL.year <- Hmsc::HmscRandomLevel(units = levels(studyDesign$year))
rL.mission <- Hmsc::HmscRandomLevel(units = levels(studyDesign$mission))
rL.site <- Hmsc::HmscRandomLevel(units = levels(studyDesign$site))


#Define unstructured random effect at spatiotemporal unit
rL.ind <- Hmsc::HmscRandomLevel(units = levels(studyDesign$individual))


# For another example, you may define location as a spatial random effect
# rL.location = HmscRandomLevel(sData = locations.xy)
# Here locations.xy would be a matrix (one row per unique location)
# where row names are the levels of studyDesign$location,
# and the columns are the xy-coordinates
# For another example, you may define the sample = sampling unit = row of matrix Y
# as a random effect, in case you are interested in co-occurrences at that level
# rL.sample = HmscRandomLevel(units = levels(studyDesign$sample))


## Candidate models specification ----
#using hmsc constructor

m.null <- Hmsc::Hmsc(Y = Y,
                     distr = "probit", 
                     XData = X,
                     XFormula = XFormula.null,
                     studyDesign = studyDesign, 
                     # TrData = TP,
                     # TrFormula = TrFormula,
                     ranLevels = list(site = rL.site, mission = rL.mission, individual = rL.ind))

m.host <- Hmsc::Hmsc(Y = Y,
                      distr = "probit", 
                      XData = X,
                      XFormula = XFormula.host,
                      studyDesign = studyDesign, 
                      # TrData = TP,
                      # TrFormula = TrFormula,
                      ranLevels = list(site = rL.site, mission = rL.mission, individual = rL.ind))


# Combining and saving models 
models = list(m.null,
              m.host)

names(models) = c("nullmodel",
                  "hostmodel")

save(models, file = here::here(modelDir, "unfitted_models.RData"))



# Test run ----
#Just to see if would fit without errors (before doing real fitting procedure)
for(i in 1:length(models)){
 print(i)
 sampleMcmc(models[[i]],samples=2)
}
