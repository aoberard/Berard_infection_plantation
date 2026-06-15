# Script parameters ----

#	INPUT : Unfitted models.

#	OUTPUT : Fitted models, with fitting done for multiple RUNs:
# first short MCMC chains (to provide some results fast), and then with increasingly long MCMC chains
# (until MCMC convergence or computational limit is reached): the files
# "models/models_thin_1_samples_5_chains_4.Rdata" (RUN 0),
# "models/models_thin_1_samples_250_chains_4.Rdata" (RUN 1),
# "models/models_thin_10_samples_250_chains_4.Rdata" (RUN 2), 
# "models/models_thin_100_samples_250_chains_4.Rdata" (RUN 3), and so on.


# Remove saved elements (clean workspace)
rm(list = ls())
gc()

# Load libraries
library("dplyr")
library("Hmsc")

# Make script reproducible
set.seed(1)

## Parallel chains parameter ----
nParallel = NULL #Default: nParallel = nChains

#Change default options by removing comment and setting value
#Only section of the script typically modified ▲

#nParallel = 1 #Set to 1 to disable parallel computing, otherwise comment

#Use socket clusters in parallel processing (TRUE or FALSE), must be TRUE if launched in Windows, else FALSE is usually quicker
useSocket_param <- FALSE



# Set directories 
dataDir <- here::here("data", "derived-data", "hmsc_out")
modelDir <- file.path(dataDir, "models")

# Load unfitted models ----
load(file = here::here(modelDir, "unfitted_models.RData"))


# Model fitting ----

## Fitting parameters ----
#Choose fitting parameters
nm <- length(models)
samples_list <- c(250, 250)    #250, 250
thin_list <- c(1, 100)          #100, 1000
nChains <- 4

if(is.null(nParallel)) nParallel = nChains
Lst = 1

#And save them for future scripts
save(samples_list, thin_list, nChains, file = here::here(modelDir, "fitting_parameters.RData"))


#Explanation of MCMC fitting parameters:
#nChains: number of independent Monte Carlo Markov Chain (MCMC)
#samples: number of posterior samples retained for each chain
#thin   : only one state every "thin" iterations is retained

#transient: "burn-in" : number of initial iteration used for each Monte Carlo. 
#Initial iterations are discarded because the chain has not yet reached the stationary posterior distribution.

#Default value in this script:
#transient = ceiling(0.5 * samples * thin)
for(i in seq_along(samples_list)) {
  transient <- ceiling(0.5 * samples_list[i] * thin_list[i])
  cat("Run", i, ": samples =", samples_list[i], ", thin =", thin_list[i], ", transient =", transient, "\n")
}


#For example, with samples = 250, thin = 100, nchains = 4, using script default value for transient
#
#Burn-in: There is initially 12500 iterations performed per chain (burn-in)
#Those iterations are discarded
#
#Sampling phase: Then there is 250 * 100 additional iterations performed (25000 after the 12500 first iterations)
#
#Thinning: One state is retained every 100 iterations from the 25000 performed just before : 250 values saved per chain
#
#Total: as 4 chains are processed independently, 250 * 4 = 1000 posterior samples
#
#Posterior samples are then used to estimate posterior means and distribution, credible intervals...
#A retained sample is not a single coefficient value, 
#but a complete realization of all model parameters at a given MCMC iteration.




## Fitting execution and saving ----

#For each combination of thinning and number of posterior samples (batch),
#all models contained in the list "models" are fitted sequentially.
#After all models have finished, the updated "models" list containing 
#fitted models is saved to disk.

#Overwrite parameter : set to TRUE to refit models even if the output file already exists
overwrite <- FALSE



while(Lst <= length(samples_list)){
  
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  
  cat("\n\n=============== NEW BATCH =================\n")
  cat("MCMC settings:\n")
  cat("  thin    = ", thin, "\n")
  cat("  samples = ", samples, "\n")  
  
  filename = file.path(modelDir,
                       paste("models_thin_", as.character(thin),
                             "_samples_", as.character(samples),
                             "_chains_", as.character(nChains),
                             ".Rdata", sep = ""))
  
  cat("Output file:\n")
  cat("  ", filename, "\n")
  
  if(file.exists(filename)  && !overwrite){
    cat("Status: file already exists, skipping this batch.\n")
    
    
  } else {
    
    batch_start <- Sys.time()
    cat("Batch started at:", format(batch_start), "\n")
    
    for (mi in 1:nm) {
      
      model_name <- names(models)[mi]
      
      cat("\n-------------------------------------------------\n")
      cat("Starting model ", mi, "/", nm, ": ", model_name, "\n", sep = "")
      model_start <- Sys.time()
      cat("Start time:", format(model_start), "\n")
      
      m = models[[mi]]
      m = sampleMcmc(m, 
                     samples = samples,
                     thin = thin,
                     adaptNf = rep(ceiling(0.4*samples*thin), m$nr), 
                     transient = ceiling(0.5*samples*thin),
                     nChains = nChains,
                     nParallel = nParallel,
                     useSocket = useSocket_param
                     ) 
      
      models[[mi]] = m
      
      model_end <- Sys.time()
      cat("Finished at:", format(model_end), "\n")
      cat("Duration   :", round(difftime(model_end, model_start, units = "mins"), 2),
          "minutes\n")
      
    }
    save(models, file = filename)
    cat("\n - ✓ Batch output saved ✓ - \n")
    cat("Saved object: `models` (named list of all fitted models)\n")
    cat("Saved to:", filename, "\n")
    
    batch_end <- Sys.time()
    cat("Batch finished at:", format(batch_end), "\n")
    cat("Total duration   :",
        round(difftime(batch_end, batch_start, units = "mins"), 2),
        "mins\n")
    
  }
  
  Lst = Lst + 1
}
