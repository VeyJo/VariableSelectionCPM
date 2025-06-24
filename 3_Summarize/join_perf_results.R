##----------------------------------------------------------------------------##
##  summarize results of performance measures and join them together
##----------------------------------------------------------------------------##

library(tidyr)
library(dplyr)

path_perf <- "Results/performance_raw/"
path_sum <- "Results/results_summarized/"

methods <- c("lmss", "mfp", "enet", "rfb", "rfh", "gbt", "gbm")
# meth <- methods[1]


#--------------------------------------------#
### tables containing performance results ####
#--------------------------------------------#

ls_res_t <- list()
## loop across methods
for (i in 1:length(methods)) {
  meth <- methods[i]
  ## read files containing results as tables
  file_names_t <- list.files(paste0(path_perf, meth, "/"), 
                             pattern = "^PerfTable.*\\.rds$", full.names = TRUE)
  l_tables <- lapply(file_names_t, readRDS)
  ## join all scenarios together
  ls_res_t[[i]] <- do.call(rbind, l_tables)
}

## join all methods together 
res_table_A <- do.call(rbind, ls_res_t)



### True Oracle
## results of true oracle does not contain all columns, therefore separate joining of results
file_names_t <- list.files(paste0(path_perf, "tom/"), 
                           pattern = "^PerfTable.*\\.rds$", full.names = TRUE)
l_tables <- lapply(file_names_t, readRDS)
## join all scenarios together
res_table_A_oracle <- do.call(rbind, l_tables)

res_table_A_all <- bind_rows(res_table_A, res_table_A_oracle)



saveRDS(res_table_A_all, 
        paste0(path_sum, "results_perf_measures_A.rds"))



### oracle models
## results of oracle does not contain all columns, therefore separate joining of results
file_names_t <- list.files(paste0(path_perf, "oracle_models/"), 
                           pattern = "^PerfTable.*\\.rds$", full.names = TRUE)
l_tables <- lapply(file_names_t, readRDS)
## join all scenarios together
res_table_A_oracles <- do.call(rbind, l_tables)

saveRDS(res_table_A_oracles, 
        paste0(path_sum, "results_perf_measures_A_oracles.rds"))



#-------------------------------------------#
### lists containing performance results ####
#-------------------------------------------#

## parallelize across methods because it is computationally intensive
library(doFuture)
library(iterators)
plan(multisession, workers = 7)


ls_res_l = foreach(i = 1:length(methods), #.combine = "c",
                   .options.future = list(packages = c("tidyr", "dplyr"))
                  ) %dofuture% {
  meth <- methods[i]
  ## read files containing results as list
  file_names_l <- list.files(paste0(path_perf, meth, "/"), 
                             pattern = "^PerfList.*\\.rds$", full.names = TRUE) 
  l_lists <- lapply(file_names_l, readRDS)
  ## name list entry to corresponding scenario parameters
  names(l_lists) <- gsub(".*_compl(.*?)\\.rds", "\\1", basename(file_names_l))
  
  ## convert the data frames of each list entry to the wide format 
  l_1 <- lapply(l_lists, function(x) lapply(x, function(d)
    ## if method did not converged, the list entry is null
    if(!is.null(d)) {
      d |> 
        mutate(var = rownames(d),
               converged = TRUE) |> 
        select(var, converged, selection, rank_m) |> 
        rename(sel = selection, rank = rank_m) |> 
        pivot_wider(names_from = var, 
                    values_from = c(sel, rank),
                    names_glue = "{var}_{.value}")
    } else {
      ## this are the names that result if repetition has converged (code above)
      ## according to the coding used for respective method 
      if(meth %in% c("lmss", "mfp", "enet", "gbm", "tom")) {
        d <- as.data.frame(matrix(NA, nrow = 1, ncol = 35))
        names(d) <- c("converged", "x1_sel", "x2_sel", "x3_sel", "x4_1_sel", "x4_2_sel",
                      "x5_sel", "x6_sel", "x7_sel", "x8_sel", "x9_1_sel", "x9_2_sel", 
                      "x10_sel", "x11_sel", "x12_sel", "x13_sel", "x14_sel", "x15_sel",
                      "x1_rank", "x2_rank", "x3_rank", "x4_1_rank", "x4_2_rank",
                      "x5_rank", "x6_rank", "x7_rank", "x8_rank", "x9_1_rank", "x9_2_rank",
                      "x10_rank", "x11_rank", "x12_rank", "x13_rank", "x14_rank", "x15_rank")
      } else if(meth %in% c("rfb", "rfh", "gbt")) {
        d <- as.data.frame(matrix(NA, nrow = 1, ncol = 31))
        names(d) <- c("converged", "x1_sel", "x2_sel", "x3_sel", "x4_sel", 
                      "x5_sel", "x6_sel", "x7_sel", "x8_sel", "x9_sel", 
                      "x10_sel", "x11_sel", "x12_sel", "x13_sel", "x14_sel", "x15_sel",
                      "x1_rank", "x2_rank", "x3_rank", "x4_rank",
                      "x5_rank", "x6_rank", "x7_rank", "x8_rank", "x9_rank", 
                      "x10_rank", "x11_rank", "x12_rank", "x13_rank", "x14_rank", "x15_rank")
      }
      d$converged <- FALSE
      return(d)
    }
  )) 
  
  ## rbind all repetitions of a scenario 
  l_2 <- lapply(l_1, function(x) do.call(rbind, x))
  l_3 <- l_2
  ## add columns with scenario parameters - loop across scenarios
  for (s in 1:length(l_2)) {
    scen <- names(l_2)[[s]]
    df <- l_2[[s]]
    df$method <- meth
    df$compl  <- sub("^(.)_.*", "\\1", scen)
    df$R2     <- sub("^.*_R0(\\d+)_.*", "0.\\1", scen)
    df$n      <- sub("^.*_n(\\d+)$", "\\1", scen)
    df$rep    <- 1:nrow(df)
    ## rearrange columns
    df <- df |> select(method, compl, R2, n, rep, everything())
    l_3[[s]] <- df
  }
  ## join all scenarios together and return resulting data.frame
  do.call(rbind, l_3)
}

## join the parametric methods
res_table_B1 <- do.call(rbind, ls_res_l[c(1,2,3,7)])
## join the non-parametric methods
res_table_B2 <- do.call(rbind, ls_res_l[c(4,5,6)])





### True Oracle Model
meth <- "tom"
## read files containing results as list
file_names_l <- list.files(paste0(path_perf, meth, "/"), 
                           pattern = "^PerfList.*\\.rds$", full.names = TRUE) 
l_lists <- lapply(file_names_l, readRDS)
## name list entry to corresponding scenario parameters
names(l_lists) <- gsub(".*_compl(.*?)\\.rds", "\\1", basename(file_names_l))

## convert the data frames of each list entry to the wide format 
l_1 <- lapply(l_lists, function(x) lapply(x, function(d)
  ## all repetitions converged 
    d |> 
      mutate(var = rownames(d),
             converged = TRUE) |> 
      select(var, converged, rank_m) |> 
      rename(rank = rank_m) |> 
      pivot_wider(names_from = var, 
                  values_from = rank,
                  names_glue = "{var}_{.value}")
  )) 

## rbind all repetitions of a scenario 
l_2 <- lapply(l_1, function(x) do.call(rbind, x))
l_3 <- l_2
## add columns with scenario parameters - loop across scenarios
for (s in 1:length(l_2)) {
  scen <- names(l_2)[[s]]
  df <- l_2[[s]]
  df$method <- meth
  df$compl  <- sub("^(.)_.*", "\\1", scen)
  df$R2     <- sub("^.*_R0(\\d+)_.*", "0.\\1", scen)
  df$n      <- sub("^.*_n(\\d+)$", "\\1", scen)
  df$rep    <- 1:nrow(df)
  ## rearrange columns
  df <- df |> select(method, compl, R2, n, rep, everything())
  ## fix column name in compl D
  if(df$compl[1] == "D"){
    df <- df |> rename(x11_rank = `c("x8", "x11")_rank`)
  }
  l_3[[s]] <- df
}
## join all scenarios together and return resulting data.frame
res_table_B3 <- do.call(rbind, l_3)





saveRDS(res_table_B1, 
        paste0(path_sum, "results_perf_measures_B_regression.rds"))
saveRDS(res_table_B2, 
        paste0(path_sum, "results_perf_measures_B_trees.rds"))
saveRDS(res_table_B3, 
        paste0(path_sum, "results_perf_measures_B_TrueOracle.rds"))
