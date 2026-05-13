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


# Explore data ----

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



#Just make sure matrix X, Y, and the studyDesign and locations.xy dataframes/matrix have aligned row names
stopifnot(
  identical(rownames(X), rownames(Y)),
  identical(rownames(X), rownames(studyDesign)),
  identical(rownames(X), rownames(locations.xy))
)
message("Row names of X, Y, studyDesign and locations.xy are aligned.")


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
rL.site <- Hmsc::HmscRandomLevel(units = levels(studyDesign$site))


#Define unstructured random effect at spatiotemporal unit
rL.ind <- Hmsc::HmscRandomLevel(units = levels(studyDesign$individual))


#Define location spatial random effect
#Add a tiny jitter to duplicated coordinates only (as sData locations must be unique)
dup <- duplicated(locations.xy) | duplicated(locations.xy, fromLast = TRUE)

locations.xy[dup, ] <- locations.xy[dup, ] +
  matrix(
    rnorm(sum(dup) * 2, mean = 0, sd = 0.01),  # 1 cm SD
    ncol = 2
  )

# Check that all coordinates are now unique
stopifnot(!anyDuplicated(locations.xy))

rL.location = HmscRandomLevel(sData = locations.xy)

# For another example, you may define location as a spatial random effect
# rL.location = HmscRandomLevel(sData = locations.xy)
# Here locations.xy would be a matrix (one row per unique location)
# where row names are the levels of studyDesign$location,
# and the columns are the xy-coordinates
# For another example, you may define the sample = sampling unit = row of matrix Y
# as a random effect, in case you are interested in co-occurrences at that level
# rL.sample = HmscRandomLevel(units = levels(studyDesign$sample))


## Phylogenetic tree ----

#Construction of a simple taxonomical tree as proxy for phylogenetic tree

taxonomicTree <- ape::as.phylo(~ domain / phylum / class / order / family / genus / names, 
                               data = Tm)




## Candidate models specification ----
#using hmsc constructor

m.null <- Hmsc::Hmsc(Y = Y,
                     distr = "probit", 
                     XData = X,
                     XFormula = XFormula.null,
                     studyDesign = studyDesign, 
                     # TrData = Tm,
                     # TrFormula = TrFormula,
                     ranLevels = list(site = rL.site, mission = rL.mission, individual = rL.ind))

m.host <- Hmsc::Hmsc(Y = Y,
                      distr = "probit", 
                      XData = X,
                      XFormula = XFormula.host,
                      studyDesign = studyDesign, 
                      # TrData = Tm,
                      # TrFormula = TrFormula,
                      ranLevels = list(site = rL.site, mission = rL.mission, individual = rL.ind))


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
