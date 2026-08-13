# Global parameters ----

#	INPUT. the Fitted models.

# OUTPUT: Parameter estimates illustrated (for highest RUN of among fitted models) in the file
# "results/parameter_estimates.pdf", the text file "results/parameters/parameter_estimates.txt", 
# as well as given numerically in multiple csv files (one per parameter type) named 
# "results/parameter_estimates_[parameter_name].csv".


# Remove saved elements (clean workspace)
rm(list = ls())
gc()

# Load libraries
library("dplyr")
library("ggplot2")
library("Hmsc")
library("colorspace")
library("corrplot")
library("writexl")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("ggtree", quietly = TRUE))
BiocManager::install("ggtree")

# Make script reproducible
set.seed(1)


# Graphical parameters
effect_direction_palette <- c("positive"          = "#2C7796",
                              "negative"          = "#D55763",
                              "non_significant"   = "gray60")


# Load saved elements ----

#Setting directories
dataDir <- here::here("data", "derived-data", "hmsc_ab_out")
modelDir <- file.path(dataDir, "models")
resultDir <- file.path(dataDir, "results", "parameters")

if(!dir.exists(resultDir)) {
  dir.create(resultDir, recursive = TRUE)
}

#Loading fitting parameters
load(file = here::here(modelDir, "fitting_parameters.RData"))


#Loading "best" MCMC model list
nst = length(thin_list)

text.file <- file.path(resultDir,"parameter_estimates.txt")
cat(c("This file contains additional information regarding parameter estimates.","\n","\n", sep=""), file=text.file)

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
  cat(c("\n",filename,"\n","\n"), file=text.file, sep="", append = TRUE)
}


## Selecting model of choice ▲ ----
#Keeping only selected models from the list of models (models object)
#Make as a comment in order to explore all models, but last part of the script will not work
# if(file.exists(filename)){
# 
#   models <- models[c("habitat_season_int")]
# }

nm = length(models)





# Parameters ----

#Setting Default parameter (Null implies default)

support.level.beta = NULL #Default: 0.95
support.level.gamma = NULL #Default: 0.95
support.level.omega = NULL #Default: 0.9
var.part.order.explained = NULL #Default: in variance partitioning of explained variance, species are shown in the order they are in the model
var.part.order.raw = NULL #Default: in variance partitioning of raw variance, species are shown in the order they are in the model
show.sp.names.beta = NULL #Default: species names shown in beta plot if there are at most 30 species and no phylogeny
plotTree = NULL #Default: tree is plotted in Beta plot if the model includes it
omega.order = NULL #Default: species shown in the order they are in the model
show.sp.names.omega = NULL #Default: species names shown in beta plot if there are at most 30 species


#Change default parameters by removing comments and setting value
#Typically only section of the script that nee to be modified     ▲

#use support levels to selected the level of statistical support shown
#support.level.beta = 0.8
#support.level.gamma = 0.8
#support.level.omega = 0.8

#use var.part.order.explained to select which order species are shown in the raw variance partitioning
#var.part.order.raw should be a list of length the number of models. 
#for each element provide either 0 (use default);
#or a vector of species indices;
#or "decreasing" if you wish to order according to explanatory power
#var.part.order.explained = list()
#var.part.order.explained[[1]] = 0
#var.part.order.explained[[1]] = "decreasing"
#var.part.order.explained[[2]] = c(2,1)
var.part.order.explained = lapply(1:nm, function(x) "decreasing")

#use var.part.order.raw to select which order species are shown in the explained variance partitioning
#same options apply as for var.part.order.explained
#var.part.order.raw = list()
#var.part.order.raw[[1]] = "decreasing"
#var.part.order.raw[[2]] = c(1,2)
var.part.order.raw = lapply(1:nm, function(x) "decreasing")


#use show.sp.names.beta to choose to show / not show species names in betaPlot
#if given, show.sp.names.beta should be a vector with length equalling number of models
show.sp.names.beta = rep(TRUE, nm)

#use plotTree to choose to plot / not plot the tree in betaPlot
#if given, plotTree should be a vector with length equalling number of models
#plotTree = c(FALSE,FALSE)

#use omega.order to select which order species are shown in omega plots
#omega.order should be a list of length the number of models. 
#for each element provide either 0 (use default);
#or a vector of species indices;
#or "AOE" if you wish to use the angular order of the eigenvectors.
#omega.order = list()
#omega.order[[1]] = "AOE"
#omega.order[[2]] = c(2,1)
#Default: species shown in the order they are in the model
#show.sp.names.omega = c(TRUE,FALSE) #Default: species names shown in beta plot if there are at most 30 species




# Fonction custom plotVP ----
plotVP_gg <- function(VP, main = "", cols = NULL, ylim = c(0, 1)) {
  vals <- VP$vals
  nr<- nrow(vals)
  
  if (is.null(cols)) {
    cols <- colorspace::qualitative_hcl(nr, palette = "Pastel1")
  }
  
  
  row_means <- rowMeans(vals, na.rm = TRUE)
  legend_labels <- paste0(rownames(vals), " (mean=", round(row_means, 2), ")")
  
  #Making data long
  df <- as.data.frame(vals)
  df$Variable <- rownames(vals)
  df_long <- tidyr::pivot_longer(df, -Variable, 
                                 names_to = "Species", 
                                 values_to = "Proportion")
  
  #species order is from VP$vals
  df_long$Species <- factor(df_long$Species, levels = colnames(vals))
  df_long$Variable <- factor(df_long$Variable, levels = rownames(vals))
  
  named_labels <- setNames(legend_labels, rownames(vals))
  
  ggplot(df_long, aes(x = Species, y = Proportion, fill = Variable)) +
    geom_bar(stat = "identity", color = "grey40", linewidth = 0.5) +
    scale_fill_manual(
      values = setNames(cols, rownames(vals)),
      labels = named_labels
    ) +
    scale_y_continuous(limits = ylim, expand = c(0, 0)) +
    labs(title = main, y = "Proportion of variance", x = NULL, fill = NULL) +
    theme_classic() +
    theme(
      axis.text.x  = element_text(angle = 90, hjust = 1, size = 7),
      legend.position = "bottom",
      legend.text  = element_text(size = 5),
      plot.title   = element_text(size = 7),
      plot.margin  = margin(5, 5, 5, 5)
    )
}



# Execute estimates visualization ----

if(is.null(support.level.beta)) support.level.beta = 0.95
if(is.null(support.level.gamma)) support.level.gamma =  0.95
if(is.null(support.level.omega)) support.level.omega =  0.9


if(file.exists(filename)){
  if(is.null(var.part.order.explained)){
    var.part.order.explained = list()
    for(j in 1:nm) var.part.order.explained[[j]] = 0
  }
  if(is.null(var.part.order.raw)){
    var.part.order.raw = list()
    for(j in 1:nm) var.part.order.raw[[j]] = 0
  }
  if(is.null(omega.order)){
    omega.order = list()
    for(j in 1:nm) omega.order[[j]] = 0
  }
  
  modelnames = names(models)
  
  
  
  pdf(file = file.path(resultDir, "parameter_estimates.pdf"))
  
  for(j in 1:nm){
    cat(c("\n", names(models)[j], "\n", "\n"), file = text.file, sep = "", append=TRUE)
    m = models[[j]]
    if(m$XFormula=="~."){
      covariates = colnames(m$X)[-1]
    } else {
      covariates = attr(terms(m$XFormula),"term.labels")
    }
    
    if(m$nr+length(covariates)>1 & m$ns>1){
      preds = computePredictedValues(m)
      VP = computeVariancePartitioning(m)
      vals = VP$vals
      mycols = rainbow(nrow(VP$vals))
      MF = evaluateModelFit(hM=m, predY=preds)
      R2 = NULL
      if(!is.null(MF$TjurR2)){
        TjurR2 = MF$TjurR2
        vals = rbind(vals,TjurR2)
        R2=TjurR2
      }
      if(!is.null(MF$R2)){
        R2=MF$R2
        vals = rbind(vals,R2)
      }
      if(!is.null(MF$SR2)){
        R2=MF$SR2
        vals = rbind(vals,R2)
      }
      filename = file.path(resultDir, paste("parameter_estimates_VP_",modelnames[j],".csv"))
      write.csv(vals,file=filename)
      if(!is.null(VP$R2T$Beta)){
        filename = file.path(resultDir,paste("parameter_estimates_VP_R2T_Beta",modelnames[j],".csv"))
        write.csv(VP$R2T$Beta,file=filename)
      }
      if(!is.null(VP$R2T$Y)){
        filename = file.path(resultDir, paste("parameter_estimates_VP_R2T_Y",modelnames[j],".csv"))
        write.csv(VP$R2T$Y,file=filename)
      }
      
      
      #Explained variance
      if(all(var.part.order.explained[[j]] == 0)){
        c.var.part.order.explained = 1:m$ns
      } else if(all(var.part.order.explained[[j]] == "decreasing")){
        c.var.part.order.explained = order(R2, decreasing = TRUE)
      } else {
        c.var.part.order.explained = var.part.order.explained[[j]]
      }
      
      VP_plot = VP
      VP_plot$vals = VP$vals[, c.var.part.order.explained]
      
      print(plotVP_gg(
        VP   = VP_plot,
        main = paste0("Proportion of explained variance, ", modelnames[j]),
        cols = NULL,
        ylim = c(0, 1.01)) 
      )
      
      #Raw variance
      if(!is.null(R2)){
        if(all(var.part.order.raw[[j]] == 0)){
          c.var.part.order.raw = 1:m$ns
        } else if(all(var.part.order.raw[[j]] == "decreasing")){
          c.var.part.order.raw = order(R2, decreasing = TRUE)
        } else {
          c.var.part.order.raw = var.part.order.raw[[j]]
        }
        
        VPr = VP
        for(k in 1:m$ns){
          VPr$vals[, k] = R2[k] * VPr$vals[, k]
        }
        VPr$vals = VPr$vals[, c.var.part.order.raw]
        
        print(plotVP_gg(
          VP   = VPr,
          main = paste0("Proportion of raw variance, ", modelnames[j]),
          cols = NULL)
        )
      }
    }
  }
  for(j in 1:nm){
    m = models[[j]]
    if(m$nc>1){
      postBeta = getPostEstimate(m, parName = "Beta")
      filename = file.path(resultDir, paste("parameter_estimates_Beta_",modelnames[j],".xlsx"))
      me = as.data.frame(t(postBeta$mean))
      me = cbind(m$spNames,me)
      colnames(me) = c("Species",m$covNames)
      po = as.data.frame(t(postBeta$support))
      po = cbind(m$spNames,po)
      colnames(po) = c("Species",m$covNames)
      ne = as.data.frame(t(postBeta$supportNeg))
      ne = cbind(m$spNames,ne)
      colnames(ne) = c("Species",m$covNames)
      vals = list("Posterior mean"=me,"Pr(x>0)"=po,"Pr(x<0)"=ne)
      writexl::write_xlsx(vals,path = filename)
      if(is.null(show.sp.names.beta)){
        c.show.sp.names = (is.null(m$phyloTree) && m$ns<=30) 
      } else {
        c.show.sp.names = show.sp.names.beta[j]
      }
      c.plotTree = !is.null(m$phyloTree)
      if(!is.null(plotTree)){
        c.plotTree = c.plotTree & plotTree[j]
      }
      plotBeta(m, post=postBeta, supportLevel = support.level.beta, param="Sign",
               plotTree = c.plotTree,
               covNamesNumbers = c(TRUE,FALSE),
               spNamesNumbers=c(c.show.sp.names,FALSE),
               cex=c(0.6,0.6,0.8))
      mymain = paste0("BetaPlot, ",modelnames[j])
      if(!is.null(m$phyloTree)){
        mpost = convertToCodaObject(m)
        rhovals = unlist(poolMcmcChains(mpost$Rho))
        mymain = paste0(mymain,", E[rho] = ",round(mean(rhovals),2),", Pr[rho>0] = ",round(mean(rhovals>0),2))
      }
      title(main=mymain, line=2.5, cex.main=0.8)
    }
  }
  for(j in 1:nm){
    m = models[[j]]      
    if(m$nt>1 & m$nc>1){
      postGamma = getPostEstimate(m, parName="Gamma")
      plotGamma(m, post=postGamma, supportLevel = support.level.gamma, param="Sign",
                covNamesNumbers = c(TRUE,FALSE),
                trNamesNumbers=c(m$nt<21,FALSE),
                cex=c(0.6,0.6,0.8))
      title(main=paste0("GammaPlot ",modelnames[j]), line=2.5,cex.main=0.8)
    }
  }
  for(j in 1:nm){
    m = models[[j]]
    if(m$nr>0 & m$ns>1){
      OmegaCor = computeAssociations(m)
      for (r in 1:m$nr){
        toPlot = ((OmegaCor[[r]]$support>support.level.omega) + (OmegaCor[[r]]$support<(1-support.level.omega))>0)*sign(OmegaCor[[r]]$mean)
        if(is.null(show.sp.names.omega)){
          c.show.sp.names = (m$ns<=30) 
        } else {
          c.show.sp.names = show.sp.names.omega[j]
        }
        if(!c.show.sp.names){
          colnames(toPlot)=rep("",m$ns)
          rownames(toPlot)=rep("",m$ns)
        }
        if(all(omega.order[[j]]==0)){
          plotOrder = 1:m$ns
        } else {
          if(all(omega.order[[j]]=="AOE")){
            plotOrder = corrMatOrder(OmegaCor[[r]]$mean,order="AOE")
          } else {
            plotOrder = omega.order[[j]]
          }
        }
        cat(c("\n","omega.order","\n","\n"),file=text.file,sep="",append=TRUE)
        cat(m$spNames[plotOrder],file=text.file,sep="\n",append=TRUE)
        mymain = paste0("Associations, ",modelnames[j], ": ",names(m$ranLevels)[[r]])
        if(m$ranLevels[[r]]$sDim>0){
          mpost = convertToCodaObject(m)
          alphavals = unlist(poolMcmcChains(mpost$Alpha[[1]][,1]))
          mymain = paste0(mymain,", E[alpha1] = ",round(mean(alphavals),2),", Pr[alpha1>0] = ",round(mean(alphavals>0),2))
        }
        corrplot(toPlot[plotOrder,plotOrder], method = "color",
                 col=colorRampPalette(c("blue","white","red"))(3),
                 mar=c(0,0,1,0),
                 main=mymain,cex.main=0.8)
        me = as.data.frame(OmegaCor[[r]]$mean)
        me = cbind(m$spNames,me)
        colnames(me)[1] = ""
        po = as.data.frame(OmegaCor[[r]]$support)
        po = cbind(m$spNames,po)
        colnames(po)[1] = ""
        ne = as.data.frame(1-OmegaCor[[r]]$support)
        ne = cbind(m$spNames,ne)
        colnames(ne)[1] = ""
        vals = list("Posterior mean"=me,"Pr(x>0)"=po,"Pr(x<0)"=ne)
        filename = file.path(resultDir, paste("parameter_estimates_Omega_",modelnames[j],"_",names(m$ranLevels)[[r]],".xlsx"))
        writexl::write_xlsx(vals,path = filename)
      }
    }
  }
  dev.off()
}



# Manual plots ----
#BEWARE BEWARE, model used for this plots may not be the right one if multiple models from "models" object were conserved
m$call

# coefficient extraction (inspired from Hmsc::getPostEstimate function)
#Pool MCMC chains

postList_pooled <- Hmsc:::poolMcmcChains(m$postList)


## Beta coefficient ----

#Extract Beta of each chain into a list
beta_list <- lapply(postList_pooled, function(a) a[["Beta"]])

#Bind all values from each element within the list
val <- abind::abind(beta_list, along = 0)  

#Quick renaming
dimnames(val) <- list(
  NULL,
  m$covNames,
  m$spNames
)

#Compute means and distribution
beta_mean <- colMeans(val)          
beta_q025 <- apply(val, c(2,3), quantile, 0.025)
beta_q975 <- apply(val, c(2,3), quantile, 0.975)
beta_q25 <- apply(val, c(2,3), quantile, 0.25)
beta_q75 <- apply(val, c(2,3), quantile, 0.75)


#Convert each matrix into long format
df_mean <- as.data.frame(beta_mean) %>%
  tibble::rownames_to_column("covariate") %>%
  tidyr::pivot_longer(-covariate, names_to = "species", values_to = "mean")

df_q025 <- as.data.frame(beta_q025) %>%
  tibble::rownames_to_column("covariate") %>%
  tidyr::pivot_longer(-covariate, names_to = "species", values_to = "lower")

df_q975 <- as.data.frame(beta_q975) %>%
  tibble::rownames_to_column("covariate") %>%
  tidyr::pivot_longer(-covariate, names_to = "species", values_to = "upper")

df_q25 <- as.data.frame(beta_q25) %>%
  tibble::rownames_to_column("covariate") %>%
  tidyr::pivot_longer(-covariate, names_to = "species", values_to = "q25")

df_q75 <- as.data.frame(beta_q75) %>%
  tibble::rownames_to_column("covariate") %>%
  tidyr::pivot_longer(-covariate, names_to = "species", values_to = "q75")

#Generate dataframe based on beta coefficient extraction
beta_coef <- df_mean %>%
  left_join(df_q025, by = c("covariate", "species")) %>%
  left_join(df_q975, by = c("covariate", "species")) %>%
  mutate(significant_bilateral = sign(lower) == sign(upper)) %>%
  left_join(df_q25, by = c("covariate", "species")) %>%
  left_join(df_q75, by = c("covariate", "species"))


#Now we want to add information obtained using Hmsc::getPostEstimate to our new dataframe (as a check and support provider)
#Extract beta coefficients using Hmsc::getPostEstimate
postBeta_bis <- Hmsc::getPostEstimate(m, parName = "Beta")
rownames(postBeta_bis$support)    <- m$covNames
rownames(postBeta_bis$supportNeg) <- m$covNames
rownames(postBeta_bis$mean) <- m$covNames


#Make a dataframe out of those
df_support <- as.data.frame(postBeta_bis$support) %>%
  tibble::rownames_to_column("covariate") %>%
  tidyr::pivot_longer(-covariate, names_to = "species", values_to = "support") %>%
  left_join(
    as.data.frame(postBeta_bis$supportNeg) %>%
      tibble::rownames_to_column("covariate") %>%
      tidyr::pivot_longer(-covariate, names_to = "species", values_to = "supportNeg"),
    by = c("covariate", "species")) %>%
  left_join(
    as.data.frame(postBeta_bis$mean) %>%
      tibble::rownames_to_column("covariate") %>%
      tidyr::pivot_longer(-covariate, names_to = "species", values_to = "mean_getestimate"),
    by = c("covariate", "species")) %>% 
  mutate(significant_support_unilateral = support > 0.95 | supportNeg > 0.95) %>% 
  mutate(significant_support_bilateral = support > 0.975 | supportNeg > 0.975)

#Join this df_support to beta_coef
beta_coef <- beta_coef %>%
  left_join(df_support %>% select(covariate, species, mean_getestimate, support, supportNeg, significant_support_unilateral, significant_support_bilateral),
            by = c("covariate", "species"))


### Exploration plots ----

#Making a forest plot
beta_coef %>%
  filter(covariate != "(Intercept)") %>%                          #REMOVING INTERCEPT FOR VISUALIZATION PURPOSES
  ggplot(aes(x = mean, y = covariate, color = species, group = species)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = lower, xmax = upper),
                 width = 0,
                 linewidth = 0.9,
                 position = position_dodge(width = 0.85)) +
  geom_errorbarh(aes(xmin = q25, xmax = q75),
                 width = 0, linewidth = 1.2,
                 position = position_dodge(width = 0.85)) +
  geom_point(aes(fill = significant_support_unilateral,
                 color = species),
             size = 2.5, shape = 21, stroke = 1.1,
             position = position_dodge(width = 0.85)) +
  scale_color_viridis_d(option = "magma", end = 0.9) +
  scale_fill_manual(
    values = c("TRUE" = "grey20", "FALSE" = "white"),
    labels = c("TRUE" = "P ≥ 0.95", "FALSE" = "P < 0.95")
  ) +
  labs(x = "β coefficients", y = NULL, color = NULL) +
  theme_bw() +
  theme(
    strip.text = element_text(color = "white", face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right"
  )


#Same but with facet and direction color
beta_coef %>%
  filter(covariate != "(Intercept)") %>%                           #REMOVING INTERCEPT FOR VISUALIZATION PURPOSES
  mutate(effect_direction = case_when(
    significant_support_unilateral == FALSE ~ "non_significant",
    mean  > 0 ~ "positive",
    mean  < 0 ~ "negative",
    TRUE ~ "problem")) %>% 
  ggplot(aes(x = mean, y = species, color = effect_direction)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), width = 0.3) +
  geom_point(size = 2.5) +
  facet_wrap(~ covariate, 
             #scales = "free_x"
  ) +
  scale_color_manual(values = effect_direction_palette) +
  labs(x = "Coefficient β", y = NULL, color = NULL) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey40"),
    strip.text = element_text(color = "white", face = "bold"),
    legend.position = "bottom"
  )


#Small look at the distribution of one beta coefficients
# plot(density(val[, "categoryhedgerows", "Neoehrlichia_mikurensis"]))




### Synthesis plot  ----


#### Tree panel ----
#Clean the tree
tree_clean <- m$phyloTree
tree_clean$edge.length <- NULL   

#Make the ggtree (discrete one)
p_tree <- ggtree::ggtree(tree_clean, 
                         branch.length = "none",
                         size = 0.4,
                         color = "grey30") +
  ggtree::theme_tree()




#### Beta coefficient panel ----

# Make labels for covariates
covariate_labels <- c(
  # "(Intercept)" = "Intercept"
  "seasonSpring" = "Season (Spring)",
  "maturitymature" = "Maturity (mature)",
  "categoryhedgerows" = "Habitat (Hedgerows)",
  "categorypine_edge" = "Habitat (Plantation)",
  "categoryhedgerows:seasonSpring" = "Hedgerows-Spring",
  "categorypine_edge:seasonSpring" = "Plantation-Spring"
)

#Apply labels to dataframe
beta_plot_df <- beta_coef %>%
  mutate(
    covariate_label = factor(
      covariate_labels[covariate],
      levels = covariate_labels
    ))


#Make the heatmap of β coefficients

#For legend max value:
lim_ab_beta <- max(abs(beta_plot_df %>% 
                      filter(!is.na(covariate_label)) %>%  #BEWARE this will remove covariates without label
                      pull(mean)
), na.rm = TRUE)

p_ab_beta <- beta_plot_df %>%
  filter(!is.na(covariate_label)) %>%                      #BEWARE this will remove covariates without label
  ggplot(aes(x = covariate_label, y = species)) +
  geom_tile(aes(fill = mean, alpha = support > 0.75 | supportNeg > 0.75), 
            color = "white", linewidth = 0.5) +
  geom_text(
    data = . %>% filter(significant_support_unilateral),
    aes(label = "*"),
    color = "black", size = 4, vjust = 0.75
  ) +
  scale_fill_gradient2(
    low = "#2C7796",
    mid = "#F7F7F7",
    high = "#D55763",
    midpoint = 0,
    name = "β coefficient",
    limits = c(-lim_ab_beta, lim_ab_beta),
  ) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 1), guide = "none") + # Decide if badly well supported tiles are transparent (to remove noise)
  labs(x = NULL, y = NULL) +
  scale_x_discrete(position = "bottom") +
  theme_bw(base_size = 8) +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 7, angle = 90, hjust = 0.9, vjust = 0, color = "grey40"),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "bottom"
  ) 





#### Variance partitioning panel ----

#Palette of covariates
vp_palette <- c(
  "category" = "#ff9aa2",
  "category:season" = "#feb7b1", 
  "season" = "#ffdac0",   
  "Random: site" = "#dbdbdb",
  "Random: location" = "#b6b6b6",
  "Random: mission" = "#929292"
)

#Covariates labels
vp_labels <- c(
  "category" = "Habitat",
  "category:season" = "Habitat x Season",
  "season"= "Season",
  "Random: site" = "Random: site",
  "Random: location" = "Random: Spatial",
  "Random: mission" = "Random: mission"
)

#Generate long format dataframe for plotting
vp_df <- as.data.frame(VP$vals) %>%
  tibble::rownames_to_column("var_level") %>%
  tidyr::pivot_longer(-var_level, names_to = "species", values_to = "proportion") %>%
  mutate(
    var_level = factor(var_level, levels = rev(names(vp_labels))) )

#Generate var part plot
p_ab_vp <- vp_df %>%
  ggplot(aes(x = proportion, y = species, fill = var_level)) +
  geom_bar(stat = "identity", color = "grey60", linewidth = 0.5) +
  scale_fill_manual(values = vp_palette,
                    name = NULL,
                    labels = vp_labels) +
  scale_x_continuous(
    limits = c(0, 1.01),
    expand = c(0, 0),
    breaks = c(0, 0.5, 1),
    labels = c("0", "0.5", "1")
  ) +
  labs(x = "Variance proportion", 
       y = NULL) +
  theme_bw(base_size = 8) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    # plot.margin = margin(0, 5, 0, 5),
    legend.key.size  = unit(0.4, "cm"),
  )

vp_df %>%
  group_by(var_level) %>% 
  summarise(mean_prop = mean(proportion)) %>% 
  arrange(desc(mean_prop))


#### Explanatory power panel ----

#Prepare SR² data
sr2_df <- data.frame(
  species = m$spNames,
  sr2 = as.numeric(MF$SR2)
)

#Prepare R² plot
p_ab_r2 <- sr2_df %>%
  ggplot(aes(x = 1, y = species, label = round(sr2, 2))) +
  geom_text(size = 2, color = "grey40") +
  labs(x = "SR²", y = NULL) +
  scale_x_continuous(limits = c(0, 2), expand = c(0, 0)) +
  theme_bw(base_size = 8) +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border= element_blank(),
  )


#### Final assemblage ----
#using aplot ensure everything is aligned

#The problem is that I can not make legend work
#so we have to use cowplot

#no legend version
p_ab_beta_noleg <- p_ab_beta + theme(legend.position = "none")


p_ab_vp_noleg <- p_ab_vp + theme(legend.position = "none")
p_ab_r2_noleg <- p_ab_r2 + theme(legend.position = "none")

#extracted legends
leg_ab_beta <- cowplot::get_legend(p_ab_beta + theme(legend.position = "bottom"))
leg_ab_vp   <- cowplot::get_legend(p_ab_vp+ theme(legend.position = "bottom"))



#Panel assemblage
p_assembled <- p_ab_beta_noleg +coord_equal()

p_assembled <- p_assembled  %>%
  aplot::insert_left(p_tree , width = 0.1) %>%
  aplot::insert_right(p_ab_vp_noleg, width = 0.6) %>%
  aplot::insert_right(p_ab_r2_noleg + coord_equal(), width = 0.09) %>%
  aplot::set_panel_spacing(spacing = 2) %>%
  ggplotify::as.ggplot() #needed for cowplot combination


#Combine legends 
p_ab_legends <- cowplot::plot_grid(
  leg_ab_beta, leg_ab_vp,
  nrow = 1,
  rel_widths = c(0.4, 0.6)
)

#Final assemblage
p_ab_final <- cowplot::plot_grid(
  p_assembled,
  p_ab_legends,
  ncol = 1,
  rel_heights = c(1, 0.15)
)
p_ab_final


## Rho parameter ----

#Extract Rho of each chain into a list
rho_list <- lapply(postList_pooled, function(a) a[["rho"]])
rhovals <- unlist(rho_list)


mean(rhovals)
mean(rhovals > 0)
quantile(rhovals, c(0.025, 0.975))

plot(density(rhovals))


## Gamma parameter ----
gamma_list <- lapply(postList_pooled, function(a) a[["Gamma"]])

#Bind all values from each element within the list
val_gamma <- abind::abind(gamma_list, along = 0)  

#Quick renaming
dimnames(val_gamma) <- list(
  NULL, 
  m$covNames,
  m$trNames
)

#Compute means, distribution quantile and support
gamma_mean <- colMeans(val_gamma)
gamma_q025 <- apply(val_gamma, c(2,3), quantile, 0.025)
gamma_q975 <- apply(val_gamma, c(2,3), quantile, 0.975)
gamma_support <- apply(val_gamma > 0, c(2,3), mean)
gamma_supportNeg <- apply(val_gamma < 0, c(2,3), mean)

# Long format
gamma_coef <- as.data.frame(gamma_mean) %>%
  tibble::rownames_to_column("covariate") %>%
  tidyr::pivot_longer(-covariate, names_to = "trait", values_to = "mean") %>%
  left_join(
    as.data.frame(gamma_q025) %>%
      tibble::rownames_to_column("covariate") %>%
      tidyr::pivot_longer(-covariate, names_to = "trait", values_to = "lower"),
    by = c("covariate", "trait")
  ) %>%
  left_join(
    as.data.frame(gamma_q975) %>%
      tibble::rownames_to_column("covariate") %>%
      tidyr::pivot_longer(-covariate, names_to = "trait", values_to = "upper"),
    by = c("covariate", "trait")
  ) %>%
  left_join(
    as.data.frame(gamma_support) %>%
      tibble::rownames_to_column("covariate") %>%
      tidyr::pivot_longer(-covariate, names_to = "trait", values_to = "support"),
    by = c("covariate", "trait")
  ) %>%
  left_join(
    as.data.frame(gamma_supportNeg) %>%
      tibble::rownames_to_column("covariate") %>%
      tidyr::pivot_longer(-covariate, names_to = "trait", values_to = "supportNeg"),
    by = c("covariate", "trait")
  ) %>%
  mutate(
    significant_support_unilateral = support > 0.95 | supportNeg > 0.95
  )


#Apply labels to dataframe
gamma_plot_df <- gamma_coef %>%
  mutate(
    covariate_label = factor(
      covariate_labels[covariate],
      levels = covariate_labels
    ))


#Make the heatmap of Γ coefficients

#For legend max value:
lim_gamma <- max(abs(gamma_plot_df %>% 
                       filter(!is.na(covariate_label)) %>%  #BEWARE this will remove covariates without label
                       pull(mean)
), na.rm = TRUE)

p_gamma <- gamma_plot_df %>%
  filter(!is.na(covariate_label)) %>%                      #BEWARE this will remove covariates without label
  ggplot(aes(x = covariate_label, y = trait)) +
  geom_tile(aes(fill = mean, alpha = support > 0.75 | supportNeg > 0.75), 
            color = "white", linewidth = 0.4) +
  geom_text(
    data = . %>% filter(significant_support_unilateral),
    aes(label = "*"),
    color = "black", size = 4, vjust = 0.75
  ) +
  scale_fill_gradient2(
    low = "#2C7796",
    mid = "#F7F7F7",
    high = "#D55763",
    midpoint = 0,
    name = "Γ coefficient",
    limits = c(-lim_gamma, lim_gamma),
  ) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0), guide = "none") + # to make not well supported things transparent
  labs(x = NULL, y = NULL) +
  scale_x_discrete(position = "bottom") +
  coord_fixed() +
  theme_bw(base_size = 8) +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 7, angle = 90, hjust = 0.9, vjust = 0, color = "grey40"),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "bottom"
  ) 
p_gamma





gamma_plot_df %>%
  # filter(!is.na(covariate_label)) %>%
  ggplot(aes(x = mean, y = covariate_label, fill = trait, group = trait)) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = lower, xmax = upper, color = trait),
                 width = 0,
                 linewidth = 0.9,
                 position = position_dodge(width = 0.65)) +
  geom_point(aes(fill = significant_support_unilateral,
                 color = trait),
             size = 2.5, shape = 21, stroke = 1.1,
             position = position_dodge(width = 0.65)) +
  scale_color_viridis_d(option = "magma", end = 0.9, begin = 0.2) +
  scale_fill_manual(
    values = c("TRUE" = "grey20", "FALSE" = "white"),
    labels = c("TRUE" = "P ≥ 0.95", "FALSE" = "P < 0.95")
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(color = "white", face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right"
  )



## Omega parameter ----

#Extract residual associations
Omega_assoc <- Hmsc::computeAssociations(m)

#Matrix to long format (for every random effects)
omega_ab_df <- purrr::map_dfr(1:length(Omega_assoc), function(r) {
  
  omega <- Omega_assoc[[r]]
  support_neg <- 1 - omega$support
  
  as.data.frame(omega$mean) %>%
    tibble::rownames_to_column("species_1") %>%
    tidyr::pivot_longer(-species_1, names_to = "species_2", values_to = "mean") %>%
    left_join(
      as.data.frame(omega$support) %>%
        tibble::rownames_to_column("species_1") %>%
        tidyr::pivot_longer(-species_1, names_to = "species_2", values_to = "support"),
      by = c("species_1", "species_2")
    ) %>%
    mutate(
      supportNeg  = 1 - support,
      significant = support > 0.95 | supportNeg > 0.95,
      random_level = names(m$ranLevels)[r]
    )
})

head(omega_ab_df)


#Extract species order (same as tree used for beta plot)
sp_order <- ggtree::ggtree(tr = tree_clean)$data %>%
  dplyr::filter(isTip) %>%
  dplyr::arrange(y) %>%
  dplyr::pull(label)

#reorder dataframe according to species order from tree
omega_ab_plot_df <- omega_ab_df %>%
  mutate(
    species_1 = factor(species_1, levels = sp_order),
    species_2 = factor(species_2, levels = sp_order)
  ) 

#Generate heatmaps
p_omega <- omega_ab_plot_df %>%
  filter(species_1 != species_2) %>%  # to remove diagonal
  # filter(as.integer(species_1) > as.integer(species_2)) %>%  # to remove bottom triangle
  # filter(support > 0.75 | supportNeg > 0.75) %>%  # to remove noise
  ggplot(aes(x = species_2, y = species_1, fill = mean)) +
  geom_tile(aes(fill = mean, alpha = support > 0.75 | supportNeg > 0.75), 
            color = "white", linewidth = 0.4) +
  geom_text(
    data = . %>% filter(significant),
    aes(label = "*"),
    color = "black", size = 3, vjust = 0.75
  ) +
  scale_fill_gradient2(
    low = "#2C7796",
    mid = "#F7F7F7",
    high= "#D55763",
    midpoint = 0,
    limits= c(-1, 1),
    name= "Residual association"
  ) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0), guide = "none") + # to make not well supported things transparent
  facet_wrap(~ random_level) +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 8) +
  theme(
    axis.text.x= element_text(angle = 45, hjust = 1, size = 7, color = "grey40"),
    axis.text.y = element_text(size = 7, color = "grey40"),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    strip.background = element_rect(fill = "grey30"),
    strip.text = element_text(color = "white", size = 7),
    legend.position  = "bottom"
  )

p_omega


# Save outputs ----
#Save figures and elements to make assembled plot in next script

#Save plots
if(!dir.exists(here::here("figures", "parameters"))) {
  dir.create(here::here("figures", "parameters"), recursive = TRUE)
}
ggsave(filename = here::here("figures", "parameters", "parameters_beta_ab_plot.png"),
       plot = p_ab_final,
       width = 3000, height = 2500, units = "px", dpi = 300, bg = NULL)

ggsave(filename = here::here("figures", "parameters", "parameters_gamma_ab_plot_.png"),
       plot = p_gamma,
       width = 2000, height = 1500, units = "px", dpi = 300, bg = NULL)

ggsave(filename = here::here("figures", "parameters", "parameters_omega_ab_plot_.png"),
       plot = p_omega,
       width = 3000, height = 2500, units = "px", dpi = 300, bg = NULL)

#Save elements for combined plots
save(p_ab_final, 
     lim_ab_beta,
     p_ab_beta,
     p_ab_vp,
     p_ab_r2,
     omega_ab_plot_df,
     file = here::here(resultDir, "parameters_beta_ab_plot.RData"))

