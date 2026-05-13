# Script parameters ----

#INPUT : Original datafiles of the case study, placed in the data folder

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


# Data check ----

## Exploration ----
#Here done previously but always make sure no strong collinearity exists among explanatory variables used in models ! 
#Data exploration is highly advised also, notably in case some variables have weird distributions (too influential outliers p.ex)
#Or if levels in qualitative variables from X does not contain enough rows (and maybe they should regrouped)
#Same for traits, maybe some traits are not represented by enough species 
plot(X)
hist(X$nb_apo)
hist(X$poids)
plot(X$category)
plot(X$season)
count(X, season, category)

colSums(Y) %>% sort(decreasing = T)
hist(colSums(Y))

table(Tm$direct_transmission)
table(Tm$tick_transmission)


## Safeguards ----

#Just make sure matrix X, Y and the studyDesign dataframe have aligned row names
stopifnot(
  identical(rownames(X), rownames(Y)),
  identical(rownames(X), rownames(studyDesign))
)
message("Row names of X, Y and studyDesign are aligned.")


#Also make sure all species in Y are present in Tm
stopifnot(
  all(colnames(Y) %in% rownames(Tm))
)
#And that they appear in same order
Tm <- Tm[colnames(Y), , drop = FALSE]
stopifnot(
  identical(rownames(Tm), colnames(Y))
  )
message("All species in Y are present in T, and T has been reordered to match Y.")


#And make sure locations.xy match the studyDesign location

#Check uniqueness of row names
stopifnot(!anyDuplicated(rownames(locations.xy)))

#Check that all studyDesign locations are represented
stopifnot(all(as.character(studyDesign$location) %in% rownames(locations.xy)))

#Reorder to match the factor levels 
locations.xy <- locations.xy[levels(studyDesign$location), , drop = FALSE]

#Final check
stopifnot(identical(rownames(locations.xy), levels(studyDesign$location)))
message("Matrix location.xy has been reordered to match studyDesign 'location' levels")



# Set up the models ----

## Models' fixed effects formula ----
XFormula.null = ~ 1
XFormula.host = ~ sexe + poids + maturity

## Models' traits formula ----
TrFormula = ~ direct_transmission + tick_transmission + flea_transmission 



## Random effects ----

#Define unstructured random effect at spatiotemporal unit
rL.year <- Hmsc::HmscRandomLevel(units = levels(studyDesign$year))
rL.mission <- Hmsc::HmscRandomLevel(units = levels(studyDesign$mission))
rL.location <- Hmsc::HmscRandomLevel(units = levels(studyDesign$location))


#Define unstructured random effect at spatiotemporal unit
rL.ind <- Hmsc::HmscRandomLevel(units = levels(studyDesign$individual))

#Define studyDesign "location" spatial random effect
rL.location.xy = HmscRandomLevel(sData = locations.xy, sMethod = "Full")


## Phylogenetic tree ----
#I think only ultrametric tree are pertinent for HMSC (all leaves at same distances from root)

#Construction of a simple taxonomical tree as proxy for phylogenetic tree
taxonomicTree <- ape::as.phylo(~ domain / phylum / class / order / family / genus / names, 
                               data = Tm)
taxonomicTree$edge.length <- rep(1, length(taxonomicTree$edge))
plot(taxonomicTree, cex = 0.5)


## Candidate models specification ----
#using hmsc constructor

m.null <- Hmsc::Hmsc(Y = Y,
                     distr = "probit", 
                     XData = X,
                     XFormula = XFormula.null,
                     studyDesign = studyDesign, 
                     # TrData = Tm,
                     # TrFormula = TrFormula,
                     # phyloTree = taxonomicTree,
                     ranLevels = list(
                       location = rL.location, 
                       mission = rL.mission, 
                       individual = rL.ind
                       ))

m.host <- Hmsc::Hmsc(Y = Y,
                      distr = "probit", 
                      XData = X,
                      XFormula = XFormula.host,
                      studyDesign = studyDesign, 
                      # TrData = Tm,
                      # TrFormula = TrFormula,
                      # phyloTree = taxonomicTree,
                      ranLevels = list(
                        location = rL.location, 
                        mission = rL.mission,
                        individual = rL.ind
                        ))


#Check models structure
m.null
m.host
head(m.host$X) #factor should be subdivided into levels Boolean variables

# Combining and saving models ----
models = list(m.null,
              m.host)

names(models) = c("nullmodel",
                  "hostmodel")

save(models, file = here::here(modelDir, "unfitted_models.RData"))



# Test run ----
#Just to see if would fit without errors (before doing real fitting procedure)
for(i in 1:length(models)){
 print(i)
 sampleMcmc(models[[i]], samples = 2)
}
