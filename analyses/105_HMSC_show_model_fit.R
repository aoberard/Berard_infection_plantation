# Global parameters ----


#	INPUT: Model fits.

#	OUTPUT. Model fits illustrated (for highest RUN of S4) in the file "results/model_fit.pdf".


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
resultDir <- file.path(dataDir, "results", "model_fit")

if(!dir.exists(resultDir)) {
  dir.create(resultDir, recursive = TRUE)
}

#Loading fitting parameters
load(file = here::here(modelDir, "fitting_parameters.RData"))

#Loading model fit parameters
load(here::here(modelDir, "nfolds_parameter.RData"))

#Loading original data used to specify HMSC models
load(here::here(dataDir, "allData.RData"))




# Execute model fit plotting ----

nst = length(thin_list)

for (Lst in nst:1) {
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  
  filename = file.path(modelDir,paste("MF_thin_", as.character(thin),
                                      "_samples_", as.character(samples),
                                      "_chains_",as.character(nChains),
                                      "_nfolds_", as.character(nfolds),
                                      ".Rdata",sep = ""))
  if(file.exists(filename)){break}
}

if(file.exists(filename)){
  
  load(filename)
  cat("Loading file : ", filename, "\n")
  
  nm = length(MF)
  modelnames = names(MF)
  pdf(file = file.path(resultDir,paste0("/model_fit_nfolds_",nfolds,".pdf")))
  for(j in 1:nm){
    cMF = MF[[j]]
    cMFCV = MFCV[[j]]
    if(!is.null(cMF$TjurR2)){
      plot(cMF$TjurR2,cMFCV$TjurR2,xlim=c(-1,1),ylim=c(-1,1),
           xlab = "explanatory power",
           ylab = "predictive power",
           main=paste0(modelnames[j],", thin = ",
                       as.character(thin),
                       ", samples = ",as.character(samples),
                       ": Tjur R2.\n",
                       "mean(MF) = ",as.character(mean(cMF$TjurR2,na.rm=TRUE)),
                       ", mean(MFCV) = ",as.character(mean(cMFCV$TjurR2,na.rm=TRUE))))
      abline(0,1)
      abline(v=0)
      abline(h=0)
    }
    if(!is.null(cMF$R2)){
      plot(cMF$R2,cMFCV$R2,xlim=c(-1,1),ylim=c(-1,1),
           xlab = "explanatory power",
           ylab = "predictive power",
           main=paste0(modelnames[[j]],", thin = ",as.character(thin),
                       ", samples = ",as.character(samples),
                       ": R2. \n",
                       "mean(MF) = ",as.character(mean(cMF$R2,na.rm=TRUE)),
                       ", mean(MFCV) = ",as.character(mean(cMFCV$R2,na.rm=TRUE))))
      abline(0,1)
      abline(v=0)
      abline(h=0)
    }
    if(!is.null(cMF$AUC)){
      plot(cMF$AUC,cMFCV$AUC,xlim=c(0,1),ylim=c(0,1),
           xlab = "explanatory power",
           ylab = "predictive power",
           main=paste0(modelnames[[j]],", thin = ",as.character(thin),
                       ", samples = ",as.character(samples),
                       ": AUC. \n",
                       "mean(MF) = ",as.character(mean(cMF$AUC,na.rm=TRUE)),
                       ", mean(MFCV) = ",as.character(mean(cMFCV$AUC,na.rm=TRUE))))
      abline(0,1)
      abline(v=0.5)
      abline(h=0.5)
    }
    if(FALSE && !is.null(cMF$O.TjurR2)){
      plot(cMF$O.TjurR2,cMFCV$O.TjurR2,xlim=c(-1,1),ylim=c(-1,1),
           xlab = "explanatory power",
           ylab = "predictive power",
           main=paste0(modelnames[[j]],", thin = ",as.character(thin),", samples = ",as.character(samples),": O.Tjur R2"))
      abline(0,1)
      abline(v=0)
      abline(h=0)
    }
    if(FALSE && !is.null(cMF$O.AUC)){
      plot(cMF$O.AUC,cMFCV$O.AUC,xlim=c(0,1),ylim=c(0,1),
           xlab = "explanatory power",
           ylab = "predictive power",
           main=paste0(modelnames[[j]],", thin = ",as.character(thin),", samples = ",as.character(samples),": O.AUC"))
      abline(0,1)
      abline(v=0.5)
      abline(h=0.5)
    }      
    if(!is.null(cMF$SR2)){
      plot(cMF$SR2,cMFCV$SR2,xlim=c(-1,1),ylim=c(-1,1),
           xlab = "explanatory power",
           ylab = "predictive power",
           main=paste0(modelnames[[j]],", thin = ",as.character(thin),
                       ", samples = ",as.character(samples),
                       ": SR2. \n",
                       "mean(MF) = ",as.character(mean(cMF$SR2,na.rm=TRUE)),
                       ", mean(MFCV) = ",as.character(mean(cMFCV$SR2,na.rm=TRUE))))
      abline(0,1)
      abline(v=0)
      abline(h=0)
    }    
    if(FALSE && !is.null(cMF$C.SR2)){
      plot(cMF$C.SR2,cMFCV$C.SR2,xlim=c(-1,1),ylim=c(-1,1),
           xlab = "explanatory power",
           ylab = "predictive power",
           main=paste0(modelnames[[j]],", thin = ",as.character(thin),", samples = ",as.character(samples),": C.SR2"))
      abline(0,1)
      abline(v=0)
      abline(h=0)
    }  
  }
  dev.off()
}


# Relationship between predictive power and prevalence ----

#Load fitted models
nst = length(thin_list)
for (Lst in nst:1) {
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  filename_models = file.path(modelDir, paste("models_thin_", as.character(thin),
                                              "_samples_", as.character(samples),
                                              "_chains_", as.character(nChains),
                                              ".Rdata", sep = ""))
  if(file.exists(filename_models)){break}
}

if(!file.exists(filename_models)) stop("models file not found: ", filename_models)
load(filename_models)

#Generating dataframe
df_spfit <- do.call(rbind, lapply(names(MFCV), function(model_name) {
  #Extract for each model_name :
  #fit performance
  cMFCV <- MFCV[[model_name]]
  cMF <- MF[[model_name]]
  
  #Hmsc object
  m <- models[[model_name]]
  
  #Then calculate occurence of each infectious agent
  occur <- colSums(m$Y == 1, na.rm = TRUE)
  
  data.frame(
    model = model_name,
    species = m$spNames,
    occur = occur,
    effectif = nrow(m$XData),  
    TjurR2_cv = if(!is.null(cMFCV$TjurR2)) cMFCV$TjurR2 else NA,
    AUC_cv = if(!is.null(cMFCV$AUC)) cMFCV$AUC else NA,
    TjurR2 = if(!is.null(cMF$TjurR2)) cMF$TjurR2 else NA,
    AUC = if(!is.null(cMF$AUC)) cMF$AUC else NA,
    row.names = NULL
  )
}))

#Add prevalence column to datafram
df_spfit <- df_spfit %>% 
  mutate(prevalence = occur / effectif)

#Plot relationship
ggplot(df_spfit, aes(x = prevalence, y = TjurR2_cv, color = model, label = species)) +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(size = 2.5, max.overlaps = 15) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap( ~ model) +
  labs(
    title = "Predictive power vs prevalence",
    x     = "Infectious agent prevalence",
    y     = "Tjur R² (cross-validated)") +
  theme_bw() +
  theme(legend.position = "none")

ggplot(df_spfit, aes(x = prevalence, y = AUC_cv, color = model, label = species)) +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(size = 2.5, max.overlaps = 15) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap( ~ model) +
  labs(
    title = "Predictive power vs prevalence",
    x     = "Infectious agent prevalence",
    y     = "AUC (cross-validated)") +
  theme_bw() +
  theme(legend.position = "none")




# Model selection ----

#Based on mean values of parameter
#using explanatory power (from MF, based on all data) and predictive power (from MFCV, based on cross validation)
#and waic


nst <- length(thin_list)

for (Lst in nst:1) {
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  
  filename = file.path(modelDir,paste("MF_thin_", as.character(thin),
                                      "_samples_", as.character(samples),
                                      "_chains_",as.character(nChains),
                                      "_nfolds_", as.character(nfolds),
                                      ".Rdata",sep = "")  )

if (file.exists(filename)) break
}

if (!file.exists(filename)) {
  stop("NO MF_*.Rdata found")
}

cat("Loading file : ", filename, "\n")
load(filename)


#Extract model names
model_names <- names(MF)

#Generate summary result table
performances <- data.frame(
  model = model_names,
  
  # Tjur R²
  MF_TjurR2 = sapply(MF,   function(x) mean(x$TjurR2)),
  MFCV_TjurR2 = sapply(MFCV, function(x) mean(x$TjurR2)),
  
  # AUC
  MF_AUC  = sapply(MF,   function(x) mean(x$AUC)),
  MFCV_AUC  = sapply(MFCV, function(x) mean(x$AUC)),
  
  # WAIC
  WAIC =  sapply(WAIC, function(x) mean(x)),
  
  thin = thin,
  samples = samples,
  nChains = nChains,
  nfolds = nfolds,
  
  row.names = NULL,
  check.names = FALSE
)

#Add overlearning diagnostic columns
performances$Delta_TjurR2 <- performances$MF_TjurR2 - performances$MFCV_TjurR2
performances$Delta_AUC <- performances$MF_AUC - performances$MFCV_AUC

#Sort models according to MFCV R²
performances <- performances[order(-performances$MFCV_TjurR2), ]

#Print
print(performances, digits = 4)

#Save performances summary
write.csv(performances,
          file = here::here(modelDir, paste0("model_performance_summary_nfolds_", nfolds,".csv")),
          row.names = FALSE)




