##----------------------------------------------------------------------------##
##  apply the methods to the simulated datasets
##----------------------------------------------------------------------------##

## load function to apply the methods
source("2_Analyze/functions_methods.R")


## load true predictor variables
truePred <- read.csv2("2_Analyze/truePredictorMatrix.csv")

parameters <- readRDS("1_DataGeneration/parameters_scenarios.rds")
parameters <- arrange(parameters, desc(n)) ## to start with longest run time

# run parallelization to simulate data
library(doFuture)
library(iterators)
# parallel::detectCores()
plan(multisession, workers = 20)

seed <- 1900

## set path where raw simulation results should be stored (and create folder)
## this folder needs to contain subfolders for each method with names
## (tom, lmss, enet, rfb, rfh, gbm, gbt, mfp)
res_path <- "ResultsRaw/"


## each method is applied separately

#==============================================================================#
#### LMSS ######################################################################
#==============================================================================#

startTime <- Sys.time()
set.seed(seed)
sim = foreach(para = iter(parameters, by='row'), .combine = rbind,
              .options.future = list(packages = c("dplyr", "mfp2"),
                                     seed = TRUE)) %dofuture% {
   # load data of scenario
   dat_list <- readRDS(paste0("Data/", "simdata_",
                              "compl-", para$compl, "_R",
                              gsub("\\.", "", para$R2),
                              "_n", para$n, ".rds"))
   dat_test <- readRDS(paste0("Data/", "testdata_",
                              "compl-", para$compl, "_R",
                              gsub("\\.", "", para$R2),
                              ".rds"))

   res <- list()
   ## apply method in each repetition
   for (r in 1:para$reps) {
     res[[r]] <- apply_stepLin(dat = dat_list[[r]],
                               dat_test = dat_test)
   }

   saveRDS(res, paste0(res_path, "lmss/",
                       "rr_lmss_compl", para$compl, 
                       "_R", gsub("\\.", "", para$R2), 
                       "_n", para$n, ".rds"))
}
stopTime <- Sys.time()
duration_stepLin <- difftime(stopTime, startTime, units = "mins")



#==============================================================================#
#### ENET ######################################################################
#==============================================================================#

startTime <- Sys.time()
set.seed(seed)
sim = foreach(para = iter(parameters, by='row'), .combine = rbind,
              .options.future = list(packages = c("dplyr", "mfp2", "glmnet"),
                                     seed = TRUE)) %dofuture% {
   # load data of scenario
   dat_list <- readRDS(paste0("Data/", "simdata_", 
                              "compl-", para$compl, "_R", 
                              gsub("\\.", "", para$R2), 
                              "_n", para$n, ".rds"))
   dat_test <- readRDS(paste0("Data/", "testdata_", 
                              "compl-", para$compl, "_R", 
                              gsub("\\.", "", para$R2), 
                              ".rds"))
   
   res <- list()
   # apply method in each repetition
   for (r in 1:para$reps) {
     res[[r]] <- apply_enet(dat = dat_list[[r]],
                            dat_test = dat_test)
   }

   saveRDS(res, paste0(res_path, "enet/",
                       "rr_enet_compl", para$compl,
                       "_R", gsub("\\.", "", para$R2),
                       "_n", para$n, ".rds"))
   ## remove to 
   rm(dat_list, dat_test, res)
}
stopTime <- Sys.time()
duration_enet <- difftime(stopTime, startTime, units = "mins")



#==============================================================================#
#### RFB #######################################################################
#==============================================================================#

startTime <- Sys.time()
set.seed(seed)
sim = foreach(para = iter(parameters, by='row'), .combine = rbind,
              .options.future = list(packages = c("dplyr", "Boruta", "ranger"),
                                     seed = TRUE)) %dofuture% {
   # load data of scenario
   dat_list <- readRDS(paste0("Data/", "simdata_", 
                              "compl-", para$compl, "_R", 
                              gsub("\\.", "", para$R2), 
                              "_n", para$n, ".rds"))
   dat_test <- readRDS(paste0("Data/", "testdata_", 
                              "compl-", para$compl, "_R", 
                              gsub("\\.", "", para$R2), 
                              ".rds"))
   
   res <- list()
   ## apply method in each repetition
   for (r in 1:para$reps) {
     res[[r]] <- apply_rf_Bo(dat = dat_list[[r]], 
                             dat_test = dat_test)
   }
   
   saveRDS(res, paste0(res_path, "rfb/",
                       "rr_rfb_compl", para$compl, 
                       "_R", gsub("\\.", "", para$R2), 
                       "_n", para$n, ".rds"))
}
stopTime <- Sys.time()
duration_rfBo <- difftime(stopTime, startTime, units = "mins")


#==============================================================================#
#### RFH #######################################################################
#==============================================================================#

startTime <- Sys.time()
set.seed(seed)
sim = foreach(para = iter(parameters, by='row'), .combine = rbind,
              .options.future = list(packages = c("dplyr", "rfvimptest", "ranger"),
                                     seed = TRUE)) %dofuture% {
   # load data of scenario
   dat_list <- readRDS(paste0("Data/", "simdata_", 
                              "compl-", para$compl, "_R", 
                              gsub("\\.", "", para$R2), 
                              "_n", para$n, ".rds"))
   dat_test <- readRDS(paste0("Data/", "testdata_", 
                              "compl-", para$compl, "_R", 
                              gsub("\\.", "", para$R2), 
                              ".rds"))
   
   res <- list()
   ## apply method in each repetition
   for (r in 1:para$reps) {
     res[[r]] <- apply_rf_Ha(dat = dat_list[[r]], 
                             dat_test = dat_test)
   }
   
   saveRDS(res, paste0(res_path, "rfh/",
                       "rr_rfh_compl", para$compl, 
                       "_R", gsub("\\.", "", para$R2), 
                       "_n", para$n, ".rds"))
}
stopTime <- Sys.time()
duration_rfHa <- difftime(stopTime, startTime, units = "mins")


#==============================================================================#
#### GBT #######################################################################
#==============================================================================#

startTime <- Sys.time()
set.seed(seed)
sim = foreach(para = iter(parameters, by='row'), .combine = rbind,
              .options.future = list(packages = c("dplyr", "gbm3"),
                                     seed = TRUE)) %dofuture% {
   # load data of scenario
   dat_list <- readRDS(paste0("Data/", "simdata_", 
                              "compl-", para$compl, "_R", 
                              gsub("\\.", "", para$R2), 
                              "_n", para$n, ".rds"))
   dat_test <- readRDS(paste0("Data/", "testdata_", 
                              "compl-", para$compl, "_R", 
                              gsub("\\.", "", para$R2), 
                              ".rds"))
   
   res <- list()
   ## apply method in each repetition
   for (r in 1:para$reps) {
     res[[r]] <- apply_gbt(dat = dat_list[[r]], 
                           dat_test = dat_test)
   }
   
   saveRDS(res, paste0(res_path, "gbt/",
                       "rr_gbt_compl", para$compl, 
                       "_R", gsub("\\.", "", para$R2), 
                       "_n", para$n, ".rds"))
}
stopTime <- Sys.time()
duration_gbt <- difftime(stopTime, startTime, units = "mins")


#==============================================================================#
#### GBM #######################################################################
#==============================================================================#

startTime <- Sys.time()
set.seed(seed)
sim = foreach(para = iter(parameters, by='row'), .combine = rbind,
              .options.future = list(packages = c("dplyr", "mfp2", "mboost"),
                                     seed = TRUE)) %dofuture% {
   # load data of scenario
   dat_list <- readRDS(paste0("Data/", "simdata_", 
                              "compl-", para$compl, "_R", 
                              gsub("\\.", "", para$R2), 
                              "_n", para$n, ".rds"))
   dat_test <- readRDS(paste0("Data/", "testdata_", 
                              "compl-", para$compl, "_R", 
                              gsub("\\.", "", para$R2), 
                              ".rds"))
   
   res <- list()
   ## apply method in each repetition
   for (r in 1:para$reps) {
     res[[r]] <- apply_gbm(dat = dat_list[[r]], 
                           dat_test = dat_test)
   }
   
   saveRDS(res, paste0(res_path, "gbm/",
                       "rr_gbm_compl", para$compl, 
                       "_R", gsub("\\.", "", para$R2), 
                       "_n", para$n, ".rds"))
}
stopTime <- Sys.time()
duration_gbm <- difftime(stopTime, startTime, units = "mins")



#==============================================================================#
#### MFP #######################################################################
#==============================================================================#

startTime <- Sys.time()
set.seed(seed)
sim = foreach(para = iter(parameters, by='row'), .combine = rbind,
              .options.future = list(packages = c("dplyr", "mfp2"),
                                     seed = TRUE)) %dofuture% {
   # load data of scenario
   dat_list <- readRDS(paste0("Data/", "simdata_", 
                              "compl-", para$compl, "_R", 
                              gsub("\\.", "", para$R2), 
                              "_n", para$n, ".rds"))
   dat_test <- readRDS(paste0("Data/", "testdata_", 
                              "compl-", para$compl, "_R", 
                              gsub("\\.", "", para$R2), 
                              ".rds"))
   res <- list()
   ## apply method in each repetition
   for (r in 1:para$reps) {
     res[[r]] <- apply_mfp(dat = dat_list[[r]], 
                           dat_test = dat_test)
   }
   
   saveRDS(res, paste0(res_path, "mfp/",
                       "rr_mfp_compl", para$compl, 
                       "_R", gsub("\\.", "", para$R2), 
                       "_n", para$n, ".rds"))
}
stopTime <- Sys.time()
duration_mfp <- difftime(stopTime, startTime, units = "mins")




#==============================================================================#
#### Oracle Models #############################################################
#==============================================================================#


startTime <- Sys.time()
set.seed(seed)
sim = foreach(para = iter(parameters, by='row'), .combine = rbind,
              .options.future = list(packages = c("dplyr", "mfp2", "glmnet",
                                                  "ranger", "mboost", "gbm3"),
                                     seed = TRUE)) %dofuture% {
                             
   # load data of scenario
   dat_list <- readRDS(paste0("Data/", "simdata_",
                              "compl-", para$compl, "_R",
                              gsub("\\.", "", para$R2),
                              "_n", para$n, ".rds"))
   dat_test <- readRDS(paste0("Data/", "testdata_",
                              "compl-", para$compl, "_R",
                              gsub("\\.", "", para$R2),
                              ".rds"))
   
   
   res <- list()
   ## calculate oracle models in each repetition
   for (r in 1:para$reps) {
     res[[r]] <- apply_oracles(dat = dat_list[[r]],
                               dat_test = dat_test,
                               truePred = truePred)
   }
   
   saveRDS(res, paste0(res_path, "oracle_models/",
                       "rr_oracles_compl", para$compl, 
                       "_R", gsub("\\.", "", para$R2), 
                       "_n", para$n, ".rds"))
}
stopTime <- Sys.time()
duration_oracles <- difftime(stopTime, startTime, units = "mins")



#==============================================================================#
#### True Oracle Model #########################################################
#==============================================================================#

startTime <- Sys.time()
set.seed(seed)
sim = foreach(para = iter(parameters, by='row'), .combine = rbind,
              .options.future = list(packages = c("dplyr", "mfp2"),
                                     seed = TRUE)) %dofuture% {
   # load data of scenario
   dat_list <- readRDS(paste0("Data/", "simdata_",
                              "compl-", para$compl, "_R",
                              gsub("\\.", "", para$R2),
                              "_n", para$n, ".rds"))
   dat_test <- readRDS(paste0("Data/", "testdata_",
                              "compl-", para$compl, "_R",
                              gsub("\\.", "", para$R2),
                              ".rds"))
   
   ### for particular complexity setting (DGM):
   ## take formula with true predictors and functional forms
   true_frmla <- ls_true_frmlas[[para$compl]]
   ## take formula with true predictors (and functional forms) and non-predictors
   frmla_all_vars <- ls_frmlas_all_vars[[para$compl]]
   
   res <- list()
   ## apply method in each repetition
   for (r in 1:para$reps) {
     res[[r]] <- apply_true_oracle(dat = dat_list[[r]],
                                   dat_test = dat_test,
                                   true_formula = true_frmla,
                                   formula_all_vars = frmla_all_vars)
   }
   
   saveRDS(res, paste0(res_path, "tom/",
                       "rr_tom_compl", para$compl, 
                       "_R", gsub("\\.", "", para$R2), 
                       "_n", para$n, ".rds"))
 }
stopTime <- Sys.time()
duration_oracle <- difftime(stopTime, startTime, units = "mins")



