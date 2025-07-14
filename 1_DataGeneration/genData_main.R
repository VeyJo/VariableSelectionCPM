##----------------------------------------------------------------------------##
##  generate data for simulation study
##  this is the main script to generate the data for the simulation study
##----------------------------------------------------------------------------##


## import variances of random error of respective outcome-generating model (as calculated in protocol)
var_error <- read.csv("1_DataGeneration/variance_error_R2.csv", sep = ";")

## source functions for each complexity level A-D
source("1_DataGeneration/genData_functions.R")



#--------------------------#
## training data ###########
#--------------------------#

# set scenario parameters
reps <- 500       # number of repetitions
compl <- c("A", "B", "C", "D")   
R2 <- c(0.3, 0.5, 0.8) 
n <- c(100, 250, 500, 1000)

parameters <- expand.grid(compl, R2, n, reps)
colnames(parameters) <- c("compl", "R2", "n", "reps")
# saveRDS(parameters, "1_DataGeneration/parameters_scenarios.rds")
# parameters <- readRDS("1_DataGeneration/parameters_scenarios.rds")

# run parallelization to simulate data
library(doFuture)
library(iterators)
# parallel::detectCores()
plan(multisession, workers = 12)


set.seed(69120)
sim = foreach(para = iter(parameters, by='row'), .combine = rbind,
              .options.future = list(packages = "mvtnorm",
                                     seed = TRUE)) %dofuture% {

      ## import quantiles for particular complexity level
      quantiles <- readRDS(paste0("1_DataGeneration/quantiles/quantiles_",
                                  para$compl, ".rds"))                 
      
      ## set variance of epsilon for particular R^2 
      var_e <- var_error[which(var_error$R2 == para$R2), as.character(para$compl)]
      
      ## generate "rep" data sets
      dat_sim <- rep(list(NA), reps)
      ## use the particular function to generate data 
      if (para$compl == "A") {
        for (r in 1:reps) {
          dat_sim[[r]] <- generate_compl_A(
            n         = para$n, 
            var_e     = var_e, 
            quantiles = quantiles)
          }
      } else if(para$compl == "B"){
        for (r in 1:reps) {
          dat_sim[[r]] <- generate_compl_B(
            n         = para$n, 
            var_e     = var_e, 
            quantiles = quantiles)
        }
      } else if(para$compl == "C"){
        for (r in 1:reps) {
          dat_sim[[r]] <- generate_compl_C(
            n         = para$n, 
            var_e     = var_e, 
            quantiles = quantiles)
        }
      } else if(para$compl == "D"){
        for (r in 1:reps) {
          dat_sim[[r]] <- generate_compl_D(
            n         = para$n, 
            var_e     = var_e, 
            quantiles = quantiles)
        }
      }
      ## save generated data sets named by appropriate sub-scenario parameters
      saveRDS(dat_sim, paste0("Data/", "simdata_", 
                                "compl-", para$compl, "_R", 
                                gsub("\\.", "", para$R2), 
                                "_n", para$n, ".rds"))
      }





#----------------------#
## test data ###########
#----------------------#

# set scenario parameters to generate one large test data
compl <- c("A", "B", "C", "D")   
R2 <- c(0.3, 0.5, 0.8) 
n <- 100000

parameters <- expand.grid(compl, R2, n)
colnames(parameters) <- c("compl", "R2", "n")

# run parallelization to simulate data
library(doFuture)
library(iterators)
# parallel::detectCores()
plan(multisession, workers = 12)

# para <- parameters[1,]
set.seed(231024)
sim = foreach(para = iter(parameters, by='row'), .combine = rbind,
              .options.future = list(packages = "mvtnorm",
                                     seed = TRUE)) %dofuture% {
     
     ## import quantiles for particular complexity level
     quantiles <- readRDS(paste0("1_DataGeneration/quantiles/quantiles_",
                                 para$compl, ".rds"))                 
     
     ## set variance of epsilon for particular R^2 
     var_e <- var_error[which(var_error$R2 == para$R2), as.character(para$compl)]
     
     ## generate "rep" data sets
     dat_sim <- NA
     ## use the particular function to generate data 
     if (para$compl == "A") {
         dat_sim <- generate_compl_A(
           n         = para$n, 
           var_e     = var_e, 
           quantiles = quantiles)
         
     } else if(para$compl == "B"){
         dat_sim <- generate_compl_B(
           n         = para$n, 
           var_e     = var_e, 
           quantiles = quantiles)
       
     } else if(para$compl == "C"){
         dat_sim <- generate_compl_C(
           n         = para$n, 
           var_e     = var_e, 
           quantiles = quantiles)
       
     } else if(para$compl == "D"){
         dat_sim <- generate_compl_D(
           n         = para$n, 
           var_e     = var_e, 
           quantiles = quantiles)
     }
     ## save generated data sets named by appropriate sub-scenario parameters
     saveRDS(dat_sim, paste0("Data/", "testdata_", 
                             "compl-", para$compl, "_R", 
                             gsub("\\.", "", para$R2), 
                             ".rds"))
   }

