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
dataDir <- here::here("data", "derived-data", "hmsc_ab_out")
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
    if(!is.null(cMF$SR2)){
      plot(cMF$SR2,cMFCV$SR2,xlim=c(-1,1),ylim=c(-1,1),
           xlab = "explanatory power",
           ylab = "predictive power",
           main=paste0(modelnames[j],", thin = ",
                       as.character(thin),
                       ", samples = ",as.character(samples),
                       ": Tjur R2.\n",
                       "mean(MF) = ",as.character(mean(cMF$SR2,na.rm=TRUE)),
                       ", mean(MFCV) = ",as.character(mean(cMFCV$SR2,na.rm=TRUE))))
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
  total_occur <- colSums(m$Y, na.rm = TRUE)
  line_occur <- colSums(m$Y > 0, na.rm = TRUE)
  
  data.frame(
    model = model_name,
    species = m$spNames,
    line_occur = line_occur,
    total_occur = total_occur,
    effectif = nrow(m$XData),  
    SR2 = if(!is.null(cMF$SR2)) cMF$SR2 else NA,
    SR2_c = if(!is.null(cMF$C.SR2)) cMF$C.SR2 else NA,

    SR2_cv = if(!is.null(cMFCV$SR2)) cMFCV$SR2 else NA,
    SR2_cv_c = if(!is.null(cMFCV$C.SR2)) cMFCV$C.SR2 else NA,
    
    row.names = NULL
  )
}))

#Add site-session prevalence column to dataframe
df_spfit <- df_spfit %>% 
  mutate(site_prevalence = line_occur / effectif)

#Plot relationship
ggplot(df_spfit, aes(x = total_occur, y = SR2_cv, color = model, label = species)) +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(size = 2.5, max.overlaps = 15) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap( ~ model) +
  labs(
    title = "Predictive power vs occurence",
    x     = "Infectious agent total number of occurence",
    y     = "SR² (cross-validated)") +
  theme_bw() +
  theme(legend.position = "none")

ggplot(df_spfit, aes(x = total_occur, y = SR2_cv_c, color = model, label = species)) +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(size = 2.5, max.overlaps = 15) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap( ~ model) +
  labs(
    title = "Predictive power vs occurence",
    x     = "Infectious agent total number of occurence",
    y     = "SR² conditional counts (cross-validated)") +
  theme_bw() +
  theme(legend.position = "none")

#Save performance per taxa
data.table::fwrite(x = df_spfit, file = here::here(resultDir, "taxa_performance.csv"))


#Save occurence - explanatory perfomance relationship for model of choice
relationship_plot <- ggplot(df_spfit %>% filter(model == "habitat_season_int"),
                            aes(x = total_occur, y = SR2, color = model, label = species)) +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(size = 2.5, max.overlaps = 15) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  labs(
    x = "Infectious agent total number of occurence",
    y = "SR²") +
  theme_bw() +
  theme(legend.position = "none")

ggsave(filename = here::here("figures", "performance_relationship_ab_plot.png"),
       plot = relationship_plot,
       width = 1500, height = 1500, units = "px", dpi = 300, bg = NULL)



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
  
  # SR²
  MF_SR2 = sapply(MF, function(x) mean(x$SR2, na.rm = TRUE)),
  MFCV_SR2 = sapply(MFCV, function(x) mean(x$SR2, na.rm = TRUE)),
  
  # SR² on conditional counts
  MF_SR2_c = sapply(MF, function(x) mean(x$C.SR2, na.rm = TRUE)),
  MFCV_SR2_c = sapply(MFCV, function(x) mean(x$C.SR2, na.rm = TRUE)),
  
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
performances$Delta_SR2 <- performances$MF_SR2 - performances$MFCV_SR2
performances$Delta_SR2_c <- performances$MF_SR2_c - performances$MFCV_SR2_c

#Sort models according to MFCV SR²
performances <- performances[order(-performances$MF_SR2), ]

#Print
print(performances, digits = 4)

#Save performances summary
write.csv(performances,
          file = here::here(modelDir, paste0("model_performance_summary_nfolds_", nfolds,".csv")),
          row.names = FALSE)



# Save species ordered by CV SR2 ----
# for each model
species_order <- lapply(names(MFCV), function(model_name) {
  cMFCV <- MFCV[[model_name]]
  cMF <- MF[[model_name]]
  m <- models[[model_name]]
  
  res <- list()
  
  #CV SR
  if(!is.null(cMFCV$SR2)){
    sr_cv <- setNames(cMFCV$SR2, m$spNames)
    res$sr_cv <- sort(sr_cv, decreasing = TRUE)
  }
  
  #Full data sr R
  if(!is.null(cMF$SR2)){
    sr <- setNames(cMF$SR2, m$spNames)
    res$sr <- sort(sr, decreasing = TRUE)
  }
  
  #Species order based on CV metric
  if(!is.null(res$sr_cv)){
    res$species_order <- names(res$sr_cv)
  } else if(!is.null(res$sr)){
    res$species_order <- names(res$sr)
  }
  return(res)
  
})

names(species_order) <- names(MFCV)
species_ab_order <-species_order 
save(species_ab_order, file = here::here(resultDir, "species_order_srCV.RData"))




