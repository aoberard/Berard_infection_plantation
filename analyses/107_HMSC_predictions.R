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

#Graphical parameters
category_palette <- c("pine_edge"          = "#3f826d",
                      "hedgerows"          = "#f2d0a4",
                      "broadleaved_forest" = "#a4a8d1")

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
        Gradient = constructGradient(m, focalVariable = covariate)
        Gradient2 = constructGradient(m, focalVariable = covariate, non.focalVariables = 1)
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


## Generate manual gradient ----

#New table of covariates used for prediction (X matrix)
#Here Season x category table is the effect of interest

#other variables are fixed (we are considering marginal effects)
XNew <- expand.grid(
  season   = levels(m$XData$season),
  category = levels(m$XData$category),
  sexe     = levels(m$XData$sexe)[1],
  maturity = levels(m$XData$maturity)[2],
  poids    = mean(m$XData$poids),
  nb_apo   = mean(m$XData$nb_apo)
)

#New table of random effects used for prediction
#We generate an artificial unit for every random effect
#in order to predict on a unit that does not exist in observed data
names(m$rL)

studyDesignNew <- data.frame(
  year       = factor("new_unit"),
  mission    = factor("new_unit"),
  location   = factor("new_unit"),
  individual = factor("new_unit")
)[rep(1, nrow(XNew)), ]
rownames(studyDesignNew) <- NULL

#The artificial unit is added in existing levels of hmsc random levels 
rLNew <- list()
for(r in names(m$rL)){
  rL1 <- m$rL[[r]]
  rL1$pi <- c(rL1$pi, "new_unit")
  rL1$N  <- rL1$N + 1
  rLNew[[r]] <- rL1
}


## Generate predictions and distributions ----

#For each MCMC samples (samples * nbChains)
#Predict occurence probability for every combination of XNew and every species
#So gives n_samplesMCMC matrices of probability (samples*nbChains matrices of [nrow(XNew) x n_species])
preds <- predict(m,
                 XData = XNew,
                 studyDesign = studyDesignNew,
                 ranLevels = rLNew,
                 expected = TRUE,
                 predictEtaMean = FALSE) #random effect fixed at their means (TRUE) or not (FALSE) 

#Bind all values from each element within the list of matrices into an array
pred_array <- abind::abind(preds, along = 0)  #[n_samples, nrow(XNew), n_species]

#Quick renaming
dimnames(pred_array) <- list(
  NULL,
  interaction(XNew$category, XNew$season),
  m$spNames)

#Compute means and distribution
pred_mean <- apply(pred_array, c(2,3), mean)
pred_q025 <- apply(pred_array, c(2,3), quantile, 0.025)
pred_q975 <- apply(pred_array, c(2,3), quantile, 0.975)
pred_q25 <- apply(pred_array, c(2,3), quantile, 0.25)
pred_q75 <- apply(pred_array, c(2,3), quantile, 0.75)

#Convert each matrix into long format
df_mean <- as.data.frame(pred_mean) %>%
  tibble::rownames_to_column("covariate") %>%
  tidyr::pivot_longer(-covariate, names_to = "species", values_to = "mean")

df_q025 <- as.data.frame(pred_q025) %>%
  tibble::rownames_to_column("covariate") %>%
  tidyr::pivot_longer(-covariate, names_to = "species", values_to = "lower")

df_q975 <- as.data.frame(pred_q975) %>%
  tibble::rownames_to_column("covariate") %>%
  tidyr::pivot_longer(-covariate, names_to = "species", values_to = "upper")

df_q25 <- as.data.frame(pred_q25) %>%
  tibble::rownames_to_column("covariate") %>%
  tidyr::pivot_longer(-covariate, names_to = "species", values_to = "q25")

df_q75 <- as.data.frame(pred_q75) %>%
  tibble::rownames_to_column("covariate") %>%
  tidyr::pivot_longer(-covariate, names_to = "species", values_to = "q75")

#Generate dataframe based on predicted probabilities
pred_distributions <- df_mean %>%
  left_join(df_q025, by = c("covariate", "species")) %>%
  left_join(df_q975, by = c("covariate", "species")) %>%
  left_join(df_q25, by = c("covariate", "species")) %>%
  left_join(df_q75, by = c("covariate", "species"))

#Join back some metadata for plotting
pred_distributions <- pred_distributions %>% 
  left_join(XNew %>% 
              mutate(category_season = interaction(category, season) )%>% 
              select(category_season, category, season),
            by = c("covariate" = "category_season"))

## Visualize predictions ----

pred_distributions %>% 
  ggplot(aes(x = mean, y = category, color = species, group = species)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = lower, xmax = upper),
                 width = 0,
                 linewidth = 0.9,
                 position = position_dodge(width = 0.85)) +
  geom_errorbarh(aes(xmin = q25, xmax = q75),
                 width = 0, linewidth = 1.2,
                 position = position_dodge(width = 0.85)) +
  geom_point(aes(color = species),
             size = 2.5,
             position = position_dodge(width = 0.85)) +
  scale_color_viridis_d(option = "magma", end = 0.9) +
  theme_bw() +
  theme(
    strip.text = element_text(color = "white", face = "bold"),
    strip.background = element_rect(fill = "grey40"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.5),
    legend.position = "right" ) +
  facet_wrap( .~ season, 
             #scales = "free_x"
  )

pred_distributions %>% 
  ggplot(aes(x = mean, y = species, color = category, group = category)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = lower, xmax = upper),
                 width = 0,
                 linewidth = 0.9,
                 position = position_dodge(width = 0.85)) +
  geom_errorbarh(aes(xmin = q25, xmax = q75),
                 width = 0, linewidth = 1.2,
                 position = position_dodge(width = 0.85)) +
  geom_point(aes(color = category),
             size = 2.5,
             position = position_dodge(width = 0.85)) +
  scale_color_manual(values = category_palette) +
  theme_bw() +
  theme(
    strip.text = element_text(color = "white", face = "bold"),
    strip.background = element_rect(fill = "grey40"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.5),
    legend.position = "right" ) +
  facet_wrap( .~ season, 
              #scales = "free_x"
  )


## Statistical investigation ----

### Season effect comparison across habitats ----
#For each habitat - are there differences in infection probability between seasons (for each species)

results_season <- list()

#For each habitat
for(hab in levels(XNew$category)){
  
  #Identify which row of XNew correspond to the right season and habitat
  idx_spring <- which(
    XNew$season == "Spring" & XNew$category == hab)
  
  idx_autumn <- which(
    XNew$season == "Autumn" & XNew$category == hab)
  
  #For each habitat
  #Substract autumn prediction from spring prediction (for every species we model)
  diff_samples <- pred_array[, idx_spring, ] - pred_array[, idx_autumn, ]
  
  #Save it into a dataframe saved within a list (one element per habitat)
  results_season[[hab]] <- data.frame(
    habitat = hab,
    species = colnames(diff_samples),
    diff_mean = colMeans(diff_samples),
    #P(Spring > Autumn)
    support_pos = colMeans(diff_samples > 0),
    #P(Spring < Autumn)
    support_neg = colMeans(diff_samples < 0),
    #Significance
    significant = if_else(colMeans(diff_samples > 0) >= 0.95 | colMeans(diff_samples < 0) >= 0.95, true = TRUE, false = FALSE)
    
  )
}

df_season <- do.call(rbind, results_season)


### Habitat effect comparison across season ----
#For each season - are there differences in infection probability between habitats (for each species)

hab_pairs <- combn(levels(XNew$category), 2, simplify = FALSE)
results_habitat <- list()
k <- 1

#For each season
for(season_i in levels(XNew$season)){
  
  #For each combination of habitats
  for(pair in hab_pairs){
    
    hab1 <- pair[1]
    hab2 <- pair[2]
    
    #Identify which row of XNew correspond to the right season and habitat OF THE PAIR of interest
    idx1 <- which(
      XNew$season == season_i & XNew$category == hab1)
    
    idx2 <- which(
      XNew$season == season_i & XNew$category == hab2)
    
    #For each season, and each pair of habitat
    #Substract prediction between habitat from the pair 
    diff_samples <- pred_array[ , idx1, ] - pred_array[ , idx2, ]
    
    #Save it into a dataframe saved within a list (one element per pair of habitat)
    results_habitat[[k]] <- data.frame(
      season = season_i,
      habitat1 = hab1,
      habitat2 = hab2,
      species = m$spNames,
      
      diff_mean = colMeans(diff_samples),
      
      support_pos = colMeans(diff_samples > 0),
      support_neg = colMeans(diff_samples < 0),
      significant = if_else(colMeans(diff_samples > 0) >= 0.95 | colMeans(diff_samples < 0) >= 0.95, true = TRUE, false = FALSE))
    
    k <- k + 1
  }
}
df_habitat <- do.call(rbind, results_habitat)

df_habitat
df_habitat %>% filter(species == "Neoehrlichia_mikurensis")



### Interaction effect ----
#Testing if the interaction between season and habitat is significant

results_interaction <- list()
k <- 1

#For each combination of habitats
for(pair in hab_pairs){
  
  hab1 <- pair[1]
  hab2 <- pair[2]
  
  #Identify which row of XNew correspond to the right season and habitat OF THE PAIR of interest
  #For each habitat both season will be extracted as we are interested about interaction
  idx_h1_spring <- which(
    XNew$category == hab1 & XNew$season == "Spring")
  
  idx_h1_autumn <- which(
    XNew$category == hab1 & XNew$season == "Autumn")
  
  idx_h2_spring <- which(
    XNew$category == hab2 & XNew$season == "Spring")
  
  idx_h2_autumn <- which(
    XNew$category == hab2 & XNew$season == "Autumn")
  
#For each habitat of the pair, compute difference of probability between seasons
  #For first habitat of the pair
  delta_h1 <-
    pred_array[, idx_h1_spring, ] - pred_array[, idx_h1_autumn, ]
  #For second habitat of the pair
  delta_h2 <-
    pred_array[, idx_h2_spring, ] - pred_array[, idx_h2_autumn, ]
  
  #Then compute the differences of differences :
  #Difference between season among habitat1 vs difference between season among habitat2
  interaction_samples <-
    delta_h1 - delta_h2
  
  #Save this into a dataframe within a list
  results_interaction[[k]] <- data.frame(
    habitat1 = hab1,
    habitat2 = hab2,
    species = m$spNames,
    
    interaction_mean = colMeans(interaction_samples),
    
    support_pos = colMeans(interaction_samples > 0),
    support_neg = colMeans(interaction_samples < 0),
    significant = if_else(colMeans(interaction_samples > 0) >= 0.95 | colMeans(interaction_samples < 0) >= 0.95,
                          true = TRUE, false = FALSE))
  
  k <- k + 1
}

df_interaction <- do.call(rbind, results_interaction)
rownames(df_interaction) <- NULL
df_interaction
 

#Plot interaction support
df_interaction %>%
  mutate(
    pair_label = paste(habitat1, "vs", habitat2),
    support_display = pmax(support_pos, support_neg)) %>%
  ggplot(aes(x = species, y = pair_label, fill = support_display)) +
  geom_tile(color = "white") +
  geom_text(aes(label = if_else(significant, "*", "")), 
            size = 7, vjust = 0.8) +
  scale_fill_viridis_c(option = "cividis", end = 0.9,
    name = "Support") +
  labs(x = NULL, y = NULL, 
       title = "Maximum interaction season x category support ") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    panel.grid = element_blank())
