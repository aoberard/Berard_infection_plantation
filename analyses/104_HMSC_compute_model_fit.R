# Script parameters ----


# INPUT: Fitted models.

#	OUTPUT: Model fits computed by the cross-validation, with fitting
# (which is part of cross-validation) done for multiple RUNs:
# first short MCMC chains (to provide some results fast), and then with increasingly long MCMC chains
# (up to the longest run performed in S2). The results are stored in the files
# "models/MF_thin_1_samples_5_chains_4.Rdata" (RUN 0),
# "models/MF_thin_1_samples_250_chains_4.Rdata" (RUN 1),
# "models/MF_thin_10_samples_250_chains_4.Rdata" (RUN 2), 
# "models/MF_thin_100_samples_250_chains_4.Rdata" (RUN 3), and so on.


# Remove saved elements (clean workspace)
rm(list = ls())
gc()

# Load libraries
library("dplyr")
library("Hmsc")

# Make script reproducible
set.seed(1)




## Prediction parameters ----

#Setting Default parameter (Null implies default)
nfolds = NULL #Default: two-fold cross-validation
nParallel = NULL #Default: nParallel = nChains

#Change default parameters by removing comments and setting value
#Typically only section of the script that nee to be modified     ▲

nfolds = 5 #change the number of CV-folds
#nParallel = 1 #Set to 1 to disable parallel computing, otherwise comment

#Use socket clusters in parallel processing (TRUE or FALSE), must be TRUE if launched in Windows, else FALSE is usually quicker
useSocket_param <- FALSE



# Load saved elements ----

#Setting directories
dataDir <- here::here("data", "derived-data", "hmsc_out")
modelDir <- file.path(dataDir, "models")


#Loading fitting parameters
load(file = here::here(modelDir, "fitting_parameters.RData"))



# Execute model fit computation ----

# Overwrite parameter: set to TRUE to recompute metrics even if the output file already exists.
overwrite <- TRUE

#Default number of folds for cross-validation
if(is.null(nfolds)) nfolds = 2

#Default number of cores used in parallel
if(is.null(nParallel)) nParallel = nChains

Lst = 1


while(Lst <= length(samples_list)){
  
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  
  cat("\n\n=============== NEW BATCH =================\n")
  cat("Model fit settings:\n")
  cat("  thin    = ", thin, "\n")
  cat("  samples = ", samples, "\n")
  cat("  nfolds  = ", nfolds, "\n")
  
  
  filename.in = file.path(modelDir,paste("models_thin_", as.character(thin),
                                         "_samples_", as.character(samples),
                                         "_chains_",as.character(nChains),
                                         ".Rdata",sep = ""))
  
  filename.out = file.path(modelDir,paste("MF_thin_", as.character(thin),
                                          "_samples_", as.character(samples),
                                          "_chains_",as.character(nChains),
                                          "_nfolds_", as.character(nfolds),
                                          ".Rdata",sep = ""))
  cat("Input file:\n")
  cat("  ", filename.in, "\n")
  cat("Output file:\n")
  cat("  ", filename.out, "\n")
  
  if (file.exists(filename.out) && !overwrite) {
    
    cat("Status: output file already exists and overwrite = FALSE.\n")
    cat("Skipping this batch.\n")
    
  } else {
    
    if (!file.exists(filename.in)) {
      cat("Status: input file not found.\n")
      cat("Skipping this batch.\n")
      
    } else {
      
      batch_start <- Sys.time()
      cat("Batch started at:", format(batch_start), "\n")
      
      load(file = filename.in) #models
      nm = length(models)
        
      MF = list()
      MFCV = list()
      WAIC = list()
        
      #Loop over models
        for(mi in 1:nm){
          
          model_name <- names(models)[mi]
          
          cat("\n-----------------------------------------------\n")
          cat("Starting model ", mi, "/", nm, ": ", model_name, "\n", sep = "")
          
          model_start <- Sys.time()
          cat("Start time:", format(model_start), "\n")
          
          m = models[[mi]]
          
          
          preds = pcomputePredictedValues(m,
                                         nParallel = nParallel,
                                         useSocket = useSocket_param
          )
          
          MF[[mi]] = evaluateModelFit(hM = m,
                                      predY = preds)
          
          partition = createPartition(m, nfolds = nfolds) #USE column = ...
          preds = pcomputePredictedValues(m,
                                         partition = partition, 
                                         nParallel = nParallel,
                                         useSocket = useSocket_param
          )
          MFCV[[mi]] = evaluateModelFit(hM = m,
                                        predY = preds)
          
          WAIC[[mi]] = computeWAIC(m)
          
          model_end <- Sys.time()
          cat("Finished at:", format(model_end), "\n")
          cat(
            "Duration   :",
            round(difftime(model_end, model_start, units = "mins"), 2),
            "minutes\n")
          
          }
        
        names(MF)=names(models)
        names(MFCV)=names(models)
        names(WAIC)=names(models)
        
        save(MF, MFCV, WAIC, file = filename.out)
        cat("\n - ✓ Batch output saved ✓ - \n")
        cat("Saved objects: `MF`, `MFCV`, and `WAIC`\n")
        cat("Saved to:", filename.out, "\n")
        
        batch_end <- Sys.time()
        cat("Batch finished at:", format(batch_end), "\n")
        cat(
          "Total duration   :",
          round(difftime(batch_end, batch_start, units = "mins"), 2),
          "mins\n")
    }
  }
  
  Lst = Lst + 1
}

#Save parameter for show model fit script
save(nfolds,
     file = here::here(modelDir, "nfolds_parameter.RData")
     )
