# Global parameters ----


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



#Loading "best" MCMC model list
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
}


## Selecting model of choice ----
#Keeping only selected models from the list of models (models object)
if(file.exists(filename)){
  
  models <- models[c("habitat&hostmodel")]
  
  nm = length(models)
}



# Prediction parameters ----

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




# Execute prediction ----

if(file.exists(filename)){
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




# Attempt to predict with interaction ----

# Grille season x category, autres variables fixées
XNew <- expand.grid(
  season   = levels(m$XData$season),
  category = levels(m$XData$category),
  sexe     = levels(m$XData$sexe)[1],
  maturity = levels(m$XData$maturity)[2],
  poids    = mean(m$XData$poids),
  nb_apo   = mean(m$XData$nb_apo)
)

# studyDesign avec nouvelles unités pour chaque random level
studyDesignNew <- data.frame(
  year       = factor("new_unit"),
  mission    = factor("new_unit"),
  location   = factor("new_unit"),
  individual = factor("new_unit")
)[rep(1, nrow(XNew)), ]
rownames(studyDesignNew) <- NULL

# Gradient minimal (comme constructGradient le fait)
rLNew <- list()
for(r in names(m$rL)){
  rL1 <- m$rL[[r]]
  rL1$pi <- c(rL1$pi, rep("new_unit", nrow(XNew)))
  rL1$N  <- rL1$N + nrow(XNew)
  rLNew[[r]] <- rL1
}

# Prédictions (effets aléatoires fixés à leur moyenne = 0)
preds <- predict(m,
                 XData          = XNew,
                 studyDesign    = studyDesignNew,
                 ranLevels      = rLNew,
                 expected       = TRUE,
                 predictEtaMean = TRUE)

# preds = liste de n_samples matrices [nrow(XNew) x 12 espèces]
# Résumer les posterieurs
pred_array <- abind::abind(preds, along = 0)  # [n_samples, nrow(XNew), 12]

pred_mean <- apply(pred_array, c(2,3), mean)
pred_q025 <- apply(pred_array, c(2,3), quantile, 0.025)
pred_q975 <- apply(pred_array, c(2,3), quantile, 0.975)

# Nommer
colnames(pred_mean) <- m$spNames
colnames(pred_q025) <- m$spNames
colnames(pred_q975) <- m$spNames







# Mettre pred_mean en format lisible d'abord
rownames(pred_mean) <- interaction(XNew$category, XNew$season)
rownames(pred_q025) <- interaction(XNew$category, XNew$season)
rownames(pred_q975) <- interaction(XNew$category, XNew$season)



# Comparaisons Spring vs Autumn pour chaque habitat
habitats <- levels(m$XData$category)

results <- list()
for(hab in habitats){
  idx_spring <- which(XNew$season == "Spring" & XNew$category == hab)
  idx_autumn <- which(XNew$season == "Autumn" & XNew$category == hab)
  
  # Différence pour chaque sample MCMC
  diff_samples <- pred_array[, idx_spring, ] - pred_array[, idx_autumn, ]
  # [1000 x 12 espèces]
  
  # Support bayésien : P(Spring > Autumn)
  support <- colMeans(diff_samples > 0)
  support_neg <- colMeans(diff_samples < 0)
  diff_mean <- colMeans(diff_samples)
  
  results[[hab]] <- data.frame(
    habitat    = hab,
    species    = m$spNames,
    diff_mean  = diff_mean,
    support    = support,
    supportNeg = support_neg,
    significant = support > 0.95 | support_neg > 0.95
  )
}

df_results <- do.call(rbind, results)
df_results %>% filter(significant)












# 
# 
# 
# # Effet de season pour chaque habitat fixé
# results_gradient <- list()
# 
# for(hab in levels(m$XData$category)){
#   
#   Gradient <- constructGradient(
#     m,
#     focalVariable = "season",
#     non.focalVariables = list(
#       category = list(3, hab),      # type 3 = valeur fixée
#       poids    = list(1, NULL),     # type 1 = fixé à la moyenne
#       sexe     = list(1, NULL),
#       maturity = list(1, NULL),
#       nb_apo   = list(1, NULL)
#     )
#   )
#   
#   preds <- predict(m, Gradient = Gradient, 
#                    expected = TRUE, 
#                    predictEtaMean = TRUE)
#   
#   pred_array <- abind::abind(preds, along = 0)
#   
#   results_gradient[[hab]] <- list(
#     XDataNew   = Gradient$XDataNew,
#     pred_mean  = apply(pred_array, c(2,3), mean),
#     pred_q025  = apply(pred_array, c(2,3), quantile, 0.025),
#     pred_q975  = apply(pred_array, c(2,3), quantile, 0.975),
#     pred_array = pred_array
#   )
# }
# 
# df_comparisons <- list()
# 
# for(hab in levels(m$XData$category)){
#   
#   pred_arr <- results_gradient[[hab]]$pred_array
#   XNew_hab <- results_gradient[[hab]]$XDataNew
#   
#   idx_spring <- which(XNew_hab$season == "Spring")
#   idx_autumn <- which(XNew_hab$season == "Autumn")
#   
#   diff <- pred_arr[, idx_spring, ] - pred_arr[, idx_autumn, ]
#   
#   df_comparisons[[hab]] <- data.frame(
#     habitat    = hab,
#     species    = m$spNames,
#     diff_mean  = colMeans(diff),
#     support    = colMeans(diff > 0),
#     supportNeg = colMeans(diff < 0)
#   ) %>%
#     mutate(significant = support > 0.95 | supportNeg > 0.95)
# }
# 
# df_all <- do.call(rbind, df_comparisons)
# 
# # Effets significatifs
# df_all %>% filter(significant) %>% arrange(habitat, species)
