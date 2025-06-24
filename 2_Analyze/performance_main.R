##----------------------------------------------------------------------------##
## calculate performance measures for raw output of methods
##----------------------------------------------------------------------------##

## load functions to calculate performance measures
source("2_Analyze/functions_perf_measures.R")
## load scenario parameters
parameters <- readRDS("1_DataGeneration/parameters_scenarios.rds")
## load true predictors
truePred <- read.csv2("2_Analyze/truePredictorMatrix.csv")


### set paths
## to load raw results objects
path_resraw <- "ResultsRaw/"


## to load test data
path_dat <- "Data/"

## to load true variable importance
path_R2 <- "2_Analyze/tables_deltaR2/"

## to save performance measure objects (create folder first)
path_res <- "Results/performance_raw/"


# run parallelization to analyze results
library(doFuture)
library(iterators)
# parallel::detectCores()
plan(multisession, workers = 12)


#### methods ####
## run the foreach function for each method one after the other 
methods <- c("lmss", "mfp", "enet", "rfb", "rfh", "gbt", "gbm")
meth <- methods[6]

sim = foreach(para = iter(parameters, by='row'), .combine = rbind,
              .options.future = list(packages = c(),
                                     seed = TRUE)) %dofuture% {
   ## load raw results of particular scenario
   res <- readRDS(paste0(path_resraw, meth, "/",
                         "rr_", meth, "_compl", para$compl, 
                         "_R", gsub("\\.", "", para$R2), 
                         "_n", para$n, ".rds"))
   
   ## load test data of respective scenario
   dat_test <- readRDS(paste0(path_dat, "testdata_", "compl-", para$compl, 
                              "_R", gsub("\\.", "", para$R2), ".rds"))
   
   ## load true variable importance (as delta R^2) of respective DGM complexity 
   trueImp_dR2 <- read.csv(paste0(path_R2, "R2_compl-", para$compl, ".csv"))
   
   
   ## create table for performance measures
   perf_table <- data.frame(method = rep(meth, length(res)))
   perf_table[, "compl"] <- para$compl
   perf_table[, "R2"] <- para$R2
   perf_table[, "n"] <- para$n
   
   ## and create list for variable specific measures stored in data frame
   perf_list <- list()
   
   
   ## calculate and insert performance measure for each repetition 
   for (r in 1:length(res)) {
     perf_table[r, "rep"] <- r
     perf_table[r, "converged"] <- res[[r]]$converged
     ### variable selection measures
     ## these performance measures con only be calculated if method converged
     if (isTRUE(perf_table[r, "converged"])) {
       perf_table[r, "TIF"] <- fun_TIF(selMat = res[[r]]$selMat, truePred = truePred)
       perf_table[r, "TEF"] <- fun_TEF(selMat = res[[r]]$selMat, truePred = truePred)
       perf_table[r, "TCR"] <- fun_TCR(selMat = res[[r]]$selMat, truePred = truePred)
       perf_table[r, "MSF"] <- fun_MSF(selMat = res[[r]]$selMat, truePred = truePred)
       perf_table[r, "MS"]  <- fun_MS(selMat = res[[r]]$selMat) 
       perf_table[r, "SV"]  <- fun_SV(selMat = res[[r]]$selMat)
     } 
     ### prediction measures
     perf_table[r, "RMSPE_final"] <- fun_RMSPE(preds = res[[r]]$predFinal, obs = dat_test$y)
     if (isTRUE(perf_table[r, "converged"])) {
       perf_table[r, "RMSPE_min"]  <- fun_RMSPE(preds = res[[r]]$predMin, obs = dat_test$y)
     }
     ## calibration measures
     perf_table[r, "calA_final"] <- fun_calA(preds = res[[r]]$predFinal, obs = dat_test$y)
     perf_table[r, "calB_final"] <- fun_calB(preds = res[[r]]$predFinal, obs = dat_test$y)
     if (isTRUE(perf_table[r, "converged"])) {
       perf_table[r, "calA_min"]  <- fun_calA(preds = res[[r]]$predMin, obs = dat_test$y)
       perf_table[r, "calB_min"]  <- fun_calB(preds = res[[r]]$predMin, obs = dat_test$y)
     }
     ## variable importance
     if (isTRUE(perf_table[r, "converged"])) {
       vimp <- fun_vimp(results = res[[r]], trueImp_dR2 = trueImp_dR2) 
       perf_table[r, "tau"] <- vimp$tau
       perf_list[[r]] <- vimp$imp
     }
     ## duration
     perf_table[r, "duration"] <- res[[r]]$duration
   }
   
   ## save table and list of performance measures
   saveRDS(perf_table, paste0(path_res, meth, "/", 
                              "PerfTable_", meth, "_compl", para$compl, 
                              "_R", gsub("\\.", "", para$R2), 
                              "_n", para$n, ".rds"))
   saveRDS(perf_list, paste0(path_res, meth, "/", 
                              "PerfList_", meth, "_compl", para$compl, 
                              "_R", gsub("\\.", "", para$R2), 
                              "_n", para$n, ".rds"))
}



#### true oracle model ####
meth <- "oracle"
# oracle model: only prediction performance measures
sim = foreach(para = iter(parameters, by='row'), .combine = rbind,
              .options.future = list(packages = c(),
                                     seed = TRUE)) %dofuture% {
   ## load raw results of particular scenario
   res <- readRDS(paste0(path_resraw, meth, "/",
                         "rr_", meth, "_compl", para$compl, 
                         "_R", gsub("\\.", "", para$R2), 
                         "_n", para$n, ".rds"))
   
   ## load test data of respective scenario
   dat_test <- readRDS(paste0(path_dat, "testdata_", "compl-", para$compl, 
                              "_R", gsub("\\.", "", para$R2), ".rds"))
   
   ## load true variable importance (as delta R^2) of respective DGM complexity 
   trueImp_dR2 <- read.csv(paste0(path_R2, "R2_compl-", para$compl, ".csv"))
   
   ## create table for performance measures
   perf_table <- data.frame(method = rep(meth, length(res)))
   perf_table[, "compl"] <- para$compl
   perf_table[, "R2"] <- para$R2
   perf_table[, "n"] <- para$n
   
   ## and create list for variable specific measures stored in data frame
   perf_list <- list()
   

   ## calculate and insert performance measure for each repetition 
   for (r in 1:length(res)) {
     perf_table[r, "rep"] <- r
     
     ### prediction measures
     perf_table[r, "RMSPE_final"] <- fun_RMSPE(preds = res[[r]]$predFinal, obs = dat_test$y)
     perf_table[r, "RMSPE_min"]  <- fun_RMSPE(preds = res[[r]]$predMin, obs = dat_test$y)

     ## calibration measures
     perf_table[r, "calA_final"] <- fun_calA(preds = res[[r]]$predFinal, obs = dat_test$y)
     perf_table[r, "calB_final"] <- fun_calB(preds = res[[r]]$predFinal, obs = dat_test$y)
     perf_table[r, "calA_min"]  <- fun_calA(preds = res[[r]]$predMin, obs = dat_test$y)
     perf_table[r, "calB_min"]  <- fun_calB(preds = res[[r]]$predMin, obs = dat_test$y)
     
     ## variable importance
     vimp <- fun_vimp(results = res[[r]], trueImp_dR2 = trueImp_dR2) 
     perf_table[r, "tau"] <- vimp$tau
     perf_list[[r]] <- vimp$imp
   }
   
   ## save table and list of performance measures
   saveRDS(perf_table, paste0(path_res, meth, "/", 
                              "PerfTable_", meth, "_compl", para$compl, 
                              "_R", gsub("\\.", "", para$R2), 
                              "_n", para$n, ".rds"))
   saveRDS(perf_list, paste0(path_res, meth, "/", 
                             "PerfList_", meth, "_compl", para$compl, 
                             "_R", gsub("\\.", "", para$R2), 
                             "_n", para$n, ".rds"))
}



#### oracle models ####
meth <- "oracle_models"
# oracle model: only prediction performance measures
sim = foreach(para = iter(parameters, by='row'), .combine = rbind,
              .options.future = list(packages = c(),
                                     seed = TRUE)) %dofuture% {
   ## load raw results of particular scenario
   res <- readRDS(paste0(path_resraw, meth, "/",
                         "rr_oracles_compl", para$compl, 
                         "_R", gsub("\\.", "", para$R2), 
                         "_n", para$n, ".rds"))
   
   ## load test data of respective scenario
   dat_test <- readRDS(paste0(path_dat, "testdata_", "compl-", para$compl, 
                              "_R", gsub("\\.", "", para$R2), ".rds"))
   
   methods <- names(res[[1]])
   ## loop over oracle models
   for (m in 1:length(methods)) {
     ## create table for performance measures
     perf_table <- data.frame(method = rep(paste0(methods[m], "_oracle"), length(res)))
     perf_table[, "compl"] <- para$compl
     perf_table[, "R2"] <- para$R2
     perf_table[, "n"] <- para$n
     

     ## calculate and insert performance measure for each repetition 
     for (r in 1:length(res)) {
       perf_table[r, "rep"] <- r
       
       ### prediction measures
       perf_table[r, "RMSPE"] <- fun_RMSPE(preds = res[[r]][[m]], obs = dat_test$y)
       
       ## calibration measures
       perf_table[r, "calA"] <- fun_calA(preds = res[[r]][[m]], obs = dat_test$y)
       perf_table[r, "calB"] <- fun_calB(preds = res[[r]][[m]], obs = dat_test$y)
     }
     
     ## save table and list of performance measures
     saveRDS(perf_table, paste0(path_res, meth, "/", 
                                "PerfTable_", methods[m], "Oracle_compl", para$compl, 
                                "_R", gsub("\\.", "", para$R2), 
                                "_n", para$n, ".rds"))
   }
}

