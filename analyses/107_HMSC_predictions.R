# Script parameters ----


#	INPUT: the Fitted models.

#	OUTPUT: Predictions over environmental gradients (for highest RUN of S2) in the file
# "results/prediction/predictions.pdf".


# Remove saved elements (clean workspace)
rm(list = ls())
gc()

# Load libraries
library("dplyr")
library("Hmsc")
library("ggplot2")

# Make script reproducible
set.seed(1)




## Prediction parameters ----

#Setting Default parameter (Null implies default)

species.list = NULL #one example species shown for each model,
#selected as prevalence closest to 0.5 (probit models) or most abundant species (other models)
trait.list = NULL #community weighted mean shown for all traits
env.list = NULL #predictions constructed over all environmental gradients

#Change default parameters by removing comments and setting value
#Typically only section of the script that nee to be modified     ▲

#use species.list to select which species are used as examples for which predictions are shown
#species.list should be a list of length the number of models. 
#for each element provide either 0 (use default) or a vector of species indices
#species.list = list()
#species.list[[1]] = 0
#species.list[[2]] = c(1,2)

#use trait.list to select for which traits predictions for community weighted mean traits are shown
#trait.list should be a list of length the number of models. 
#for each element provide either 0 (use default) or a vector of trait indices
#see models[[j]]$trNames to see which trait each index corresponds to
#trait.list = list()
#trait.list[[1]] = c(2,10)
#trait.list[[2]] = 0

#use env.list to select over which environmental gradients predictions are generated
#env.list should be a list of length the number of models. 
#for each element provide either 0 (use default) or a vector of environmental variables
#env.list = list()
#env.list[[1]] = 0
#env.list[[2]] = c("tree","decay")


# Load saved elements ----

#Setting directories
dataDir <- here::here("data", "derived-data", "hmsc_out")
modelDir <- file.path(dataDir, "models")
resultDir <- file.path(dataDir, "results", "prediction")

if(!dir.exists(resultDir)) {
  dir.create(resultDir, recursive = TRUE)
}

#Loading fitting parameters
load(file = here::here(modelDir, "fitting_parameters.RData"))



# Execute prediction ----

nst = length(thin_list)

for (Lst in nst:1) {
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  filename = file.path(modelDir,paste("models_thin_", as.character(thin),
                                      "_samples_", as.character(samples),
                                      "_chains_",as.character(nChains),
                                      ".Rdata",sep = ""))
  if(file.exists(filename)){break}
}
if(file.exists(filename)){
  load(filename)
  nm = length(models)
  modelnames = names(models)
  if(is.null(species.list)){
    species.list = list()
    for(j in 1:nm) species.list[[j]] = 0
  }
  if(is.null(trait.list)){
    trait.list = list()
    for(j in 1:nm) trait.list[[j]] = 0
  }
  if(is.null(env.list)){
    env.list = list()
    for(j in 1:nm) env.list[[j]] = 0
  }
  
  pdf(file= file.path(resultDir,"predictions.pdf"))
  for(j in 1:nm){
    m = models[[j]]
    if(all(env.list[[j]]==0)){
      if(m$XFormula=="~."){
        covariates = colnames(m$XData)
      } else {
        covariates = all.vars(m$XFormula)
      }
    } else {
      covariates = env.list[[j]]
    }
    ex.sp = which.max(colMeans(m$Y,na.rm = TRUE)) #most common species as example species
    if(m$distr[1,1]==2){
      ex.sp = which.min(abs(colMeans(m$Y,na.rm = TRUE)-0.5))
    }
    if(!all(species.list[[j]])==0){
      ex.sp = species.list[[j]]
    }
    if(length(covariates)>0){
      for(k in 1:(length(covariates))){
        covariate = covariates[[k]]
        Gradient = constructGradient(m,focalVariable = covariate)
        Gradient2 = constructGradient(m,focalVariable = covariate,non.focalVariables = 1)
        predY = predict(m, Gradient=Gradient, expected = TRUE)  
        predY2 = predict(m, Gradient=Gradient2, expected = TRUE)  
        par(mfrow=c(2,1))
        pl = plotGradient(m, Gradient, pred=predY, yshow = 0, measure="S", showData = TRUE, 
                          main = paste0(modelnames[j],": summed response (total effect)"))
        if(inherits(pl, "ggplot")){
          print(pl + labs(title=paste0(modelnames[j],": summed response (total effect)")))
        }
        pl = plotGradient(m, Gradient2, pred=predY2, yshow = 0, measure="S", showData = TRUE, 
                          main = paste0(modelnames[j],": summed response (marginal effect)"))
        if(inherits(pl, "ggplot")){
          print(pl + labs(title=paste0(modelnames[j],": summed response (marginal effect)")))
        }
        for(l in 1:length(ex.sp)){
          par(mfrow=c(2,1))
          pl = plotGradient(m, Gradient, pred=predY, yshow = if(m$distr[1,1]==2){c(-0.1,1.1)}else{0}, measure="Y",index=ex.sp[l], showData = TRUE, 
                            main = paste0(modelnames[j],": example species (total effect)"))
          if(inherits(pl, "ggplot")){
            print(pl + labs(title=paste0(modelnames[j],": example species (total effect)")))
          }
          pl = plotGradient(m, Gradient2, pred=predY2, yshow = if(m$distr[1,1]==2){c(-0.1,1.1)}else{0}, measure="Y",index=ex.sp[l], showData = TRUE, 
                            main = paste0(modelnames[j],": example species (marginal effect)"))
          if(inherits(pl, "ggplot")){
            print(pl + labs(title=paste0(modelnames[j],": example species (marginal effect)")))
          }
        }
        if(m$nt>1){
          traitSelection = 2:m$nt
          if(!all(trait.list[[j]]==0)) traitSelection = trait.list[[j]]
          for(l in traitSelection){
            par(mfrow=c(2,1))
            pl = plotGradient(m, Gradient, pred=predY, measure="T",index=l, showData = TRUE,yshow = 0,
                              main = paste0(modelnames[j],": community weighted mean trait (total effect)"))
            if(inherits(pl, "ggplot")){
              print(pl + labs(title=paste0(modelnames[j],": community weighted mean trait (total effect)")))
            }
            pl = plotGradient(m, Gradient2, pred=predY2, measure="T",index=l, showData = TRUE, yshow = 0,
                              main = paste0(modelnames[j],": community weighted mean trait (marginal effect)"))
            if(inherits(pl, "ggplot")){
              print(pl + labs(title=paste0(modelnames[j],": community weighted mean trait (marginal effect)")))
            }
          }
        }
      }
    }
  }
  dev.off()
}