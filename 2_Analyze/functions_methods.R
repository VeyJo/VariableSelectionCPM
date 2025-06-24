##----------------------------------------------------------------------------##
##  functions to apply the methods
##----------------------------------------------------------------------------##

#### results table contain
# selected variables (selMat: data.frame with 1 row)
# variable importance (varImp: named numeric of selected variables)
# predicted values of developed model
# computation time 

library(dplyr)
library(mfp2)
library(glmnet)
library(Boruta)
library(ranger)
library(rfvimptest)
## load adapted version of rfvimptest() where I added ... to the arguments to disable parallelization (which is missing in the function)
source("2_Analyze/rfvimptest_adapted.R")
library(mboost)
library(gbm3)




#==============================================================================#
#### LM ########################################################################
#==============================================================================#


apply_stepLin <- function(dat, dat_test){
  ## Record start time
  startTime <- Sys.time()
  
  #### coding of categorical data 
  ## ordinal dummy coding of x4
  dat <- dat |>
    mutate(x4 = factor(x4, ordered = TRUE, levels = c(0, 1, 2)),
           x9 = factor(x9, ordered = FALSE, levels = c(0, 1, 2))) |>
    ## ordinal variable is dummy ordinal coded
    create_dummy_variables(drop_variables = TRUE,
                           var_ordinal = "x4", var_nominal = "x9") |>
    mutate(across(c(x4_1, x4_2), as.numeric)) |>
    relocate(c(x4_1, x4_2), .after = x3) |>
    relocate(c(x91, x92), .after = x8)|>
    rename(x9_1 = x91,
           x9_2 = x92)

  dat_test <- dat_test |>
    mutate(x4 = factor(x4, ordered = TRUE, levels = c(0, 1, 2)),
           x9 = factor(x9, ordered = FALSE, levels = c(0, 1, 2))) |>
    ## ordinal variable is dummy ordinal coded
    create_dummy_variables(drop_variables = TRUE,
                           var_ordinal = "x4", var_nominal = "x9") |>
    mutate(across(c(x4_1, x4_2), as.numeric)) |>
    relocate(c(x4_1, x4_2), .after = x3) |>
    relocate(c(x91, x92), .after = x8)|>
    rename(x9_1 = x91,
           x9_2 = x92)
  
  x_names <- colnames(dat)[colnames(dat) != "y"]
  
  #### variable selection
  ## track if an error occurs
  A <- tryCatch({
    ## perform stepwise forward if EPV < 10
    if (nrow(dat)/(ncol(dat)-1) < 10) {
      ## stepwise forward selection 
      fit_step <- step(
        # null model
        lm(y ~ 1, data = dat),
        # full model
        scope = formula(lm(y ~ ., data = dat)),
        direction = "both", trace = FALSE) 
    } else {
      fit_step <- step(
        # full model
        lm(y ~ ., data = dat), 
        direction = "both", trace = FALSE) 
      }
    ## collect selected variables and remove intercept
    selection <- names(fit_step[["coefficients"]])
    selection <- selection[selection != "(Intercept)"]
    },
    error = function(cond) {
      ## if error - set selection to NULL and return error object
      selection <<- NULL
      return(cond)
    } 
    )
  
  ## proceed if no error occurred and at least 1 variable was selected
  if(!any(class(A) == "error") & length(selection) > 0) {
    converged <- TRUE
    ## create selection matrix
    selMat <- data.frame(matrix(data = ifelse(x_names %in% selection, 1,0), 
                                ncol = length(x_names)))
    names(selMat) <- x_names
    
    #### variable importance
    ## standardize selected variables and y and calculate variable importance
    dat_sel <- dat[names(dat) %in% selection]
    dat_sel[,"y"] <- dat$y
    dat_sel_stand <- as.data.frame(scale(dat_sel))
    ## take absolute values of standardized regression coefficients
    varImp <- abs(coef(lm(y ~., dat_sel_stand)))
    varImp <- varImp[names(varImp) != "(Intercept)"]
    ## get the 4 most important predictors
    var4 <- names(sort(varImp, decreasing = TRUE)[1:4])
    ## only take selected variables if less than 4 variables were selected 
    var4 <- na.omit(var4)
    
    #### prediction
    ## predict test data using final and minimal model
    pred_final <- as.numeric(predict(fit_step, dat_test))

    ## minimal model with the 4 most important predictors
    fit_min <- lm(formula(paste("y~", paste(var4, collapse = "+"))), dat)
    pred_min <- as.numeric(predict(fit_min, dat_test))
    
    ## Record stop time and determine difference
    stopTime <- Sys.time()
    duration <- difftime(stopTime, startTime, units = "mins")
    ## return final list
    return(list("methName" = 'stepLin', "converged" = converged,
                "selMat" = selMat, "varImp" = varImp,
                "predFinal" = pred_final,
                "predMin" = pred_min,
                "duration" = duration))
  
    
  } else {    ### if an error occurred (method failed): 
    converged <- FALSE
    ## save corresponding error message
    if(any(class(A) == "error")) {
      error_message <- A$message
    } else if(length(selection) < 1) {
      error_message <- "no variable selected"
    } else {
      error_message <- "any other error"
    }
    ## if the method failed, only the most important predictor x5 is taken
    ## predict with model including x5 only as specified in protocol
    pred_final <- predict(lm(y ~ x5, data = dat), newdata = dat_test)
    ## all other measures are set to NA
    selMat <- NA
    varImp <- NA
    pred_min <- NA
    duration <- NA
    
    return(list("methName" = 'stepLin', 
                "converged" = converged, "errorMessage" = error_message,
                "selMat" = selMat, "varImp" = varImp,
                "predFinal" = pred_final,
                "predMin" = pred_min,
                "duration" = duration))
  }
}

# res_lm <- apply_stepLin(dat, dat_test)


#==============================================================================#
#### ENET ######################################################################
#==============================================================================#


apply_enet <- function(dat, dat_test){
  ## Record start time
  startTime <- Sys.time()
  
  ## coding of categorical data 
  dat <- dat |>
    mutate(x4 = factor(x4, ordered = TRUE, levels = c(0, 1, 2)),
           x9 = factor(x9, ordered = FALSE, levels = c(0, 1, 2))) |>
    ## ordinal variable is dummy ordinal coded
    create_dummy_variables(drop_variables = TRUE,
                           var_ordinal = "x4", var_nominal = "x9") |>
    mutate(across(c(x4_1, x4_2), as.numeric)) |>
    relocate(c(x4_1, x4_2), .after = x3) |>
    relocate(c(x91, x92), .after = x8)|>
    rename(x9_1 = x91,
           x9_2 = x92)

  dat_test <- dat_test |>
    mutate(x4 = factor(x4, ordered = TRUE, levels = c(0, 1, 2)),
           x9 = factor(x9, ordered = FALSE, levels = c(0, 1, 2))) |>
    ## ordinal variable is dummy ordinal coded
    create_dummy_variables(drop_variables = TRUE,
                           var_ordinal = "x4", var_nominal = "x9") |>
    mutate(across(c(x4_1, x4_2), as.numeric)) |>
    relocate(c(x4_1, x4_2), .after = x3) |>
    relocate(c(x91, x92), .after = x8)|>
    rename(x9_1 = x91,
           x9_2 = x92)
  

  X <- as.matrix(dat[names(dat) != "y"])
  Y <- as.numeric(dat$y)
  X_test <- as.matrix(dat_test[names(dat_test) != "y"])
  Y_test <- as.numeric(dat_test$y)
  
  x_names <- colnames(X)
  
  #### variable selection
  ## track if an error occurs
  A <- tryCatch({
    ## conduct 5-fold cross-validation to determine best alpha and lambda
    seq_alpha <- seq(from = 0.1, to = 1, by = 0.1)
    list_cv <- list()
    for(i in 1:length(seq_alpha)){
      list_cv[[i]] <- glmnet::cv.glmnet(X, Y, alpha = seq_alpha[i], 
                                        nfolds = 5, type.measure = "mse")
    }
    
    ## extract and identify combination of alpha/lambda with lowest CV mse
    cvm_min <- unlist(lapply(list_cv, function(x) 
      mse.min <- x$cvm[x$lambda == x$lambda.min]))
    k <- which.min(cvm_min) ## seq_alpha[k] is optimal alpha
    
    ## take final model from CV
    enet_cv_final <- list_cv[[k]]
    ## fit final with optimal alpha and lambda
    opt_a <- seq_alpha[k]
    opt_l <- enet_cv_final$lambda.min

    ## extract selected variables from final model
    c <- matrix(coef(enet_cv_final, s = "lambda.min"))
    ## remove intercept
    c <- c[-1, , drop = FALSE]
    ## create selection matrix
    selMat <- ifelse(c == 0, 0, 1)
    selMat <- data.frame(t(selMat))
    names(selMat) <- x_names
    selection <- names(selMat[, selMat[1,] == 1])
  }, 
  error = function(cond) {
    ## if error - set selection to NULL and return error object
    selection <<- NULL
    return(cond)
  } 
  )
  
  ## proceed if no error occurred and at least 1 variable was selected
  if(!any(class(A) == "error") & length(selection) > 0) {
    converged <- TRUE
    
    #### variable importance
    ### get standardized coefs for variable importance
    coefs_raw <- data.frame(t(c))
    names(coefs_raw) <- x_names
    ## only take selected variables
    coefs_raw <- coefs_raw[, names(coefs_raw) %in% selection]
    dat_sel <- dat[, names(dat) %in% selection]
    ## calculate sd of selected variables and calculate standardized coefs
    sd <- apply(dat_sel, 2, sd)
    ## multiply coefs by sd and convert to named numeric
    varImp <- as.numeric(abs(coefs_raw * sd))
    names(varImp) <- colnames(dat_sel)
    ## get the 4 most important predictors
    var4 <- names(sort(varImp, decreasing = TRUE)[1:4])
    ## only take selected variables if less than 4 variables were selected 
    var4 <- na.omit(var4)
    
    #### prediction
    ## predict test data using final and minimal modl
    pred_final <- as.numeric(predict(enet_cv_final, s = "lambda.min", 
                                     newx = X_test, type = "response"))

    ## minimal model with the 4 most important predictors (ridge)
    X_min <- X[, x_names %in% var4]
    enet_min <- glmnet::cv.glmnet(X_min, Y, alpha = 0)
    pred_min <- as.numeric(predict(enet_min, s = "lambda.min", 
                                   newx = X_test[, x_names %in% var4], 
                                   type = "response"))
    
    ## Record stop time and determine difference
    stopTime <- Sys.time()
    duration <- difftime(stopTime, startTime, units = "mins")
    ## return final list (when everything worked out)
    return(list("methName" = 'enet', "converged" = converged,
                "selMat" = selMat, "varImp" = varImp,
                "predFinal" = pred_final,
                "predMin" = pred_min,
                "duration" = duration))
  

  } else {    ### if an error occurred (method failed):
    converged <- FALSE
    ## save corresponding error message
    if(any(class(A) == "error")) {
      error_message <- A$message
    } else if(length(selection) < 1) {
      error_message <- "no variable selected"
    } else {
      error_message <- "any other error"
    }
    ## if the method failed, only the most important predictor x5 is taken
    ## predict with model including x5 only as specified in protocol
    ## lm is used since glmnet requires at least 2 variables
    pred_final <- predict(lm(y ~ x5, data = dat), newdata = dat_test)
    ## all other measures are set to NA
    selMat <- NA
    varImp <- NA
    pred_min <- NA
    duration <- NA
    return(list("methName" = 'enet', 
                "converged" = converged, "errorMessage" = error_message,
                "selMat" = selMat, "varImp" = varImp,
                "predFinal" = pred_final,
                "predMin" = pred_min,
                "duration" = duration))
  }
}

# res_enet <- apply_enet(dat, dat_test)


#==============================================================================#
#### RF ########################################################################
#==============================================================================#

#------------------------------------------------------------------------------#
##### Boruta ####

apply_rf_Bo <- function(dat, dat_test){
  ## Record start time
  startTime <- Sys.time()
  
  ## nominal variables are handled as factor variables
  ## ordinal variable x4 is handled as numeric to keep ordinal scale
  dat <- dat |> 
    mutate(across(.cols = c(x2, x8, x9, x14, x15), as.factor))
  x_names <- colnames(dat)[colnames(dat) != "y"]
  

  #### variable selection
  ## track if an error occurs
  A <- tryCatch({
    boruta <- Boruta::Boruta(
      x = dat[,-16], y = dat[,16],
      pValue = 0.01, 
      mcAdj = TRUE,
      maxRuns = 500,
      doTrace = 0,
      holdHistory = TRUE,
      getImp = getImpRfZ,
      ntree = 1000,
      mtry = function(x){ceiling(x/3)},   # default mtry (x/3) for regression
      num.threads = 1)                    # disable parallelization of ranger
    
    ## consider only confirmed as selected, tentative as not selected
    selection <- getSelectedAttributes(boruta) 
  }, 
  error = function(cond) {
    ## if error - set selection to NULL and return error object
    selection <<- NULL
    return(cond)
  } 
  )

    
  ## proceed if no error occurred and at least 1 variable was selected
  if(!any(class(A) == "error") & length(selection) > 0) {
    converged <- TRUE
    ## create selection matrix
    selMat <- data.frame(matrix(data = ifelse(x_names %in% selection, 1,0), 
                                ncol = length(x_names)))
    names(selMat) <- x_names
    
    ## determine best mtry to fit final model using selected variables
    dat_sel <- dat[, names(dat) %in% c(selection, "y")]
    mtry_seq <- 1:length(selection)
    list_rf <- lapply(mtry_seq, function(m) {
      ranger::ranger(y ~ ., data = dat_sel, num.trees = 1000, mtry = m,
                     importance = "permutation", num.threads = 1)
    })
    ## extract and identify mtry with lowest OOB error
    oob_error <- unlist(lapply(list_rf, function(x) x$prediction.error))
    k <- which.min(oob_error)
    mtryFin <- mtry_seq[k]
    ## take final model
    rf_Bo_final <- list_rf[[k]]
    
    #### variable importance 
    ### for all variables from Boruta procedure
    imp <- as.data.frame(boruta$ImpHistory) # scaled permutation importance from ranger
    imp <- imp[, names(imp) %in% x_names] # %in% selection
    ## calculate median importance across the runs
    varImp_all <- apply(imp, 2, function(x) median(x[is.finite(x)]))
    ### for selected variables from final model
    varImp <- rf_Bo_final$variable.importance
    ## get the 4 most important predictors
    var4 <- names(sort(varImp, decreasing = TRUE)[1:4])
    ## only take selected variables if less than 4 variables were selected 
    var4 <- na.omit(var4)
    
    
    #### prediction
    ## final developed model
    pred_final <- predict(rf_Bo_final, dat_test)[["predictions"]]
    
    ### minimal model
    mtry_seq <- 1:length(var4)
    list_rf <- lapply(mtry_seq, function(m) {
      ranger::ranger(formula(paste("y~", paste(var4, collapse = "+"))),
                     data = dat, num.trees = 1000, mtry = m, 
                     num.threads = 1)
    })
    ## extract and identify mtry with lowest OOB error
    oob_error <- unlist(lapply(list_rf, function(x) x$prediction.error))
    k <- which.min(oob_error)
    ## take final model
    rf_Bo_min <- list_rf[[k]]
    pred_min <- predict(rf_Bo_min, dat_test)[["predictions"]]
    
    ## Record stop time and determine difference
    stopTime <- Sys.time()
    duration <- difftime(stopTime, startTime, units = "mins")
    ## return final list
    return(list("methName" = 'boruta', "converged" = converged,
                "selMat" = selMat, "varImp" = varImp, "varImpBo" = varImp_all,
                "predFinal" = pred_final,
                "predMin" = pred_min,
                "duration" = duration,
                "mtryFinal" = mtryFin))
  
    
    } else {    ### if an error occurred (method failed): 
      converged <- FALSE
      ## save corresponding error message
      if(any(class(A) == "error")) {
        error_message <- A$message
      } else if(length(selection) < 1) {
        error_message <- "no variable selected"
      } else {
        error_message <- "any other error"
      }
      ## if the method failed, only the most important predictor x5 is taken
      ## predict with model including x5 only as specified in protocol
      rf_x5 <- ranger::ranger(y ~ x5,
                              data = dat, num.trees = 1000, mtry = 1, 
                              num.threads = 1)
      pred_final <- predict(rf_x5, dat_test)[["predictions"]]
      ## all other measures are set to NA
      selMat <- NA
      varImp <- NA
      pred_min <- NA
      duration <- NA
      
      return(list("methName" = 'boruta', 
                  "converged" = converged, "errorMessage" = error_message,
                  "selMat" = selMat, "varImp" = varImp,
                  "predFinal" = pred_final,
                  "predMin" = pred_min,
                  "duration" = duration))
    }
}




#------------------------------------------------------------------------------#
##### Hapfelmeier ####


apply_rf_Ha <- function(dat, dat_test){
  ## Record start time
  startTime <- Sys.time()
  
  ## nominal variables are handled as factor variables
  ## ordinal variable x4 is handled as numeric to keep ordinal scale
  dat <- dat |> 
    mutate(across(.cols = c(x2, x8, x9, x14, x15), as.factor))
  x_names <- colnames(dat)[colnames(dat) != "y"]
  
  #### variable selection
  ## track if an error occurs
  A <- tryCatch({
    ## use adapted version to disable parallelization (until this is fixed in the package)
    imptest <- rfvimptest(
      dat, yname = "y", condinf = FALSE,
      type = "SPRT", ntree = 1000, 
      mtry = ceiling(ncol(dat[,-1])/3),       # default mtry (x/3) for regression
      num.threads = 1)                        # disable parallelization of ranger 
    
    ## variables are selected if their H1 is accepted
    selection <- names(imptest$testres[imptest$testres == "accept H1"])
  }, 
  error = function(cond) {
    ## if error - set selection to NULL and return error object
    selection <<- NULL
    return(cond)
  } 
  )
  
  ## proceed if no error occurred and at least 1 variable was selected
  if(!any(class(A) == "error") & length(selection) > 0) {
    converged <- TRUE
    ## create selection matrix
    selMat <- data.frame(matrix(data = ifelse(x_names %in% selection, 1,0),
                                ncol = length(x_names)))
    names(selMat) <- x_names
    
    ## determine best mtry to fit final model using selected variables
    dat_sel <- dat[, names(dat) %in% c(selection, "y")]
    mtry_seq <- 1:length(selection)
    list_rf <- lapply(mtry_seq, function(m) {
      ranger::ranger(y ~ ., data = dat_sel, num.trees = 1000, mtry = m,
                     importance = "permutation", num.threads = 1)
    })
    ## extract and identify mtry with lowest OOB error
    oob_error <- unlist(lapply(list_rf, function(x) x$prediction.error))
    k <- which.min(oob_error)
    mtryFin <- mtry_seq[k]
    ## take final model
    rf_Ha_final <- list_rf[[k]]
                                
    
    #### variable importance
    ### for all variables from Hapfelmeier procedure 
    varImp_all <- imptest[["varimp"]]
    ### for selected variables from final model
    varImp <- rf_Ha_final$variable.importance
    ## get the 4 most important predictors
    var4 <- names(sort(varImp, decreasing = TRUE)[1:4])
    ## only take selected variables if less than 4 variables were selected 
    var4 <- na.omit(var4)
    
    #### prediction
    ## final developed model
    pred_final <- predict(rf_Ha_final, dat_test)[["predictions"]]
    
    ### minimal model
    mtry_seq <- 1:length(var4)
    list_rf <- lapply(mtry_seq, function(m) {
      ranger::ranger(formula(paste("y~", paste(var4, collapse = "+"))), 
                     data = dat, num.trees = 1000, mtry = m,
                     num.threads = 1)
    })
    ## extract and identify mtry with lowest OOB error
    oob_error <- unlist(lapply(list_rf, function(x) x$prediction.error))
    k <- which.min(oob_error)
    ## take final model
    rf_Ha_min <- list_rf[[k]]
    pred_min <- predict(rf_Ha_min, dat_test)[["predictions"]]
    
    ## Record stop time and determine difference
    stopTime <- Sys.time()
    duration <- difftime(stopTime, startTime, units = "mins")
    ## return final list
    return(list("methName" = 'Hapfelmeier', "converged" = converged,
                "selMat" = selMat, "varImp" = varImp, "varImpHa" = varImp_all,
                "predFinal" = pred_final,
                "predMin" = pred_min,
                "duration" = duration,
                "mtryFinal" = mtryFin))
    
    
  } else{    ### if an error occurred (method failed): 
    converged <- FALSE
    ## save corresponding error message
    if(any(class(A) == "error")) {
      error_message <- A$message
    } else if(length(selection) < 1) {
      error_message <- "no variable selected"
    } else {
      error_message <- "any other error"
    }
    ## if the method failed, only the most important predictor x5 is taken
    ## predict with model including x5 only as specified in protocol
    rf_x5 <- ranger::ranger(y ~ x5,
                            data = dat, num.trees = 1000, mtry = 1, 
                            num.threads = 1)
    pred_final <- predict(rf_x5, dat_test)[["predictions"]]
    ## all other measures are set to NA
    selMat <- NA
    varImp <- NA
    pred_min <- NA
    duration <- NA
    
    return(list("methName" = 'Hapfelmeier', 
                "converged" = converged, "errorMessage" = error_message,
                "selMat" = selMat, "varImp" = varImp,
                "predFinal" = pred_final,
                "predMin" = pred_min,
                "duration" = duration))
  }
}



#==============================================================================#
#### Boosting ##################################################################
#==============================================================================#

#------------------------------------------------------------------------------#
##### mboost ####


apply_gbm <- function(dat, dat_test){
  ## Record start time
  startTime <- Sys.time()
  
  ## ordinal coding of x4
  dat <- dat |>
    mutate(x4 = factor(x4, ordered = TRUE, levels = c(0, 1, 2)),
           x9 = factor(x9, ordered = FALSE, levels = c(0, 1, 2))) |>
    ## ordinal variable is dummy ordinal coded
    create_dummy_variables(drop_variables = TRUE,
                           var_ordinal = "x4", var_nominal = "x9") |>
    mutate(across(c(x4_1, x4_2), as.numeric)) |>
    relocate(c(x4_1, x4_2), .after = x3) |>
    relocate(c(x91, x92), .after = x8)|>
    rename(x9_1 = x91,
           x9_2 = x92)

  dat_test <- dat_test |>
    mutate(x4 = factor(x4, ordered = TRUE, levels = c(0, 1, 2)),
           x9 = factor(x9, ordered = FALSE, levels = c(0, 1, 2))) |>
    ## ordinal variable is dummy ordinal coded
    create_dummy_variables(drop_variables = TRUE,
                           var_ordinal = "x4", var_nominal = "x9") |>
    mutate(across(c(x4_1, x4_2), as.numeric)) |>
    relocate(c(x4_1, x4_2), .after = x3) |>
    relocate(c(x91, x92), .after = x8)|>
    rename(x9_1 = x91,
           x9_2 = x92)
  
  
  x_names <- colnames(dat)[colnames(dat) != "y"]
  
  #### variable selection
  ## track if an error occurs
  A <- tryCatch({
    fit_gbm <- mboost::glmboost(y~., data = dat,
                                 family = Gaussian(),
                                 control = boost_control(mstop = 10000,
                                                         nu = 0.01,
                                                         risk = "inbag"))
    ## calculate optimal stopping iteration using 5-fold CV
    cv <- mboost::cvrisk(fit_gbm, 
                         folds = mboost::cv(model.weights(fit_gbm), 
                                            type = "kfold", B = 5),
                         papply = lapply) # perform CV not in parallel
    ## take model at optimal iteration 
    m_stop <- mboost::mstop(cv)
    gbm_final <- fit_gbm[m_stop]
    
    #### variable importance
    ## extract variable importance (reduction) from final model at mstop 
    imp <- as.data.frame(mboost::varimp(gbm_final))
    ## remove intercept
    imp <- imp[-1, c("variable", "reduction")]
    ## take selected variables (reduction != 0)
    varImp <- imp[imp$reduction != 0, ]
    ## transform to named numeric
    varImp <- setNames(varImp[[2]], varImp[[1]])
    ## save selected variables to create selection matrix
    selection <- x_names[x_names %in% names(varImp)]
  }, 
  error = function(cond) {
    ## if error - set selection to NULL and return error object
    selection <<- NULL
    return(cond)
  } 
  )

  ## proceed if no error occurred and at least 1 variable was selected
  if(!any(class(A) == "error") & length(selection) > 0) {
    converged <- TRUE
    ## create selection matrix
    selMat <- data.frame(matrix(data = ifelse(x_names %in% selection, 1,0), 
                                ncol = length(x_names)))
    names(selMat) <- x_names
    ## get the 4 most important predictors
    var4 <- names(sort(varImp, decreasing = TRUE)[1:4])
    ## only take selected variables if less than 4 variables were selected 
    var4 <- na.omit(var4)
    
    #### prediction
    ## final developed model
    pred_final <- as.numeric(predict(gbm_final, dat_test))
    
    ## minimal model
    gbm_min <- mboost::glmboost(formula(paste("y~", paste(var4, collapse = "+"))), 
                                data = dat,
                                family = Gaussian(),
                                control = boost_control(mstop = 10000,
                                                        nu = 0.01,
                                                        risk = "inbag"))
    ## calculate optimal stopping iteration using 5-fold CV
    cv_min <- mboost::cvrisk(gbm_min, 
                             folds = mboost::cv(model.weights(fit_gbm), 
                                                type = "kfold", B = 5),
                             papply = lapply) # perform CV not in parallel
    ## take model at optimal iteration 
    pred_min <- as.numeric(predict(gbm_min[mboost::mstop(cv_min)], dat_test))
    
    ## Record stop time and determine difference
    stopTime <- Sys.time()
    duration <- difftime(stopTime, startTime, units = "mins")
    ## return final list
    return(list("methName" = 'gbm', "converged" = converged,
                "selMat" = selMat, "varImp" = varImp,
                "predFinal" = pred_final,
                "predMin" = pred_min,
                "duration" = duration,
                "mStop" = m_stop))
    
    
  } else{    ### if an error occurred (method failed): 
    converged <- FALSE
    ## save corresponding error message
    if(any(class(A) == "error")) {
      error_message <- A$message
    } else if(length(selection) < 1) {
      error_message <- "no variable selected"
    } else {
      error_message <- "any other error"
    }
    ## if the method failed, only the most important predictor x5 is taken
    ## predict with model including x5 only as specified in protocol
    fit_gbm <- mboost::glmboost(y ~ x5, data = dat,
                                family = Gaussian(),
                                control = boost_control(mstop = 10000,
                                                        nu = 0.01,
                                                        risk = "inbag"))
    ## calculate optimal stopping iteration using 5-fold CV
    cv <- mboost::cvrisk(fit_gbm, 
                         folds = mboost::cv(model.weights(fit_gbm), 
                                            type = "kfold", B = 5),
                         papply = lapply) # perform CV not in parallel
    ## take model at optimal iteration 
    m_stop <- mboost::mstop(cv)
    gbm_x5 <- fit_gbm[m_stop]
    pred_final <- as.numeric(predict(gbm_x5[mboost::mstop(cv)], dat_test))
    ## all other measures are set to NA
    selMat <- NA
    varImp <- NA
    pred_min <- NA
    duration <- NA
    
    return(list("methName" = 'gbm', 
                "converged" = converged, "errorMessage" = error_message,
                "selMat" = selMat, "varImp" = varImp,
                "predFinal" = pred_final,
                "predMin" = pred_min,
                "duration" = duration))
  }
}



#------------------------------------------------------------------------------#
##### gbm3 ####

apply_gbt <- function(dat, dat_test){
  ## Record start time
  startTime <- Sys.time()
  
  ## nominal variables are handled as factor variables
  ## ordinal variable x4 is handled as numeric to keep ordinal scale
  dat <- dat |> 
    mutate(across(.cols = c(x2, x8, x9, x14, x15), as.factor))
  x_names <- colnames(dat)[colnames(dat) != "y"]
  
  #### variable selection
  ## track if an error occurs
  A <- tryCatch({
    ## fit regression trees as base learners with gbmt function
    fit_gtb <- gbm3::gbmt(y ~ .,
                          data = dat,
                          distribution = gbm_dist("Gaussian"),
                          train_params = 
                            training_params(
                              num_trees = 10000,                   # default: 2000
                              interaction_depth = 3,               # default: 3
                              min_num_obs_in_node = 10,            # default: 10
                              shrinkage = 0.01,                    # default: 0.001
                              bag_fraction = 0.5,                  # default: 0.5
                              num_train = round(0.5 * nrow(dat)),  # default(0.5n) 
                              num_features = ncol(dat) - 1),       # default: all variables
                 cv_folds = 5,
                 par_details = gbmParallel(num_threads = 1)) # perform CV not in parallel
    ## calculate optimal stopping iteration using 5-fold CV
    best_iter <- gbm3::gbmt_performance(fit_gtb, method = "cv")
    
    #### variable importance
    ## extract varimp by the relative influence from final model at mstop 
    imp <- summary(fit_gtb, num_trees = best_iter, plot_it = FALSE)
    ## order the data by variable name
    imp <- imp[order(as.numeric(sub("x", "", imp$var))), ]
    ## take selected variables (reduction != 0)
    varImp <- imp[imp$rel_inf != 0, ]
    ## transform to named numeric
    varImp <- setNames(varImp[[2]], varImp[[1]])
    ## save selected variables to create selection matrix
    selection <- x_names[x_names %in% names(varImp)]
  }, 
  error = function(cond) {
    ## if error - set selection to NULL and return error object
    selection <<- NULL
    return(cond)
  } 
  )

  ## proceed if no error occurred and at least 1 variable was selected
  if(!any(class(A) == "error") & length(selection) > 0) {
    converged <- TRUE
    
    ## create selection matrix
    selMat <- data.frame(matrix(data = ifelse(x_names %in% selection, 1,0), 
                                ncol = length(x_names)))
    names(selMat) <- x_names
    ## get the 4 most important predictors
    var4 <- names(sort(varImp, decreasing = TRUE)[1:4])
    ## only take selected variables if less than 4 variables were selected 
    var4 <- na.omit(var4)
    
    
    #### prediction
    ## final developed model at optimal iteration
    pred_final <- as.numeric(predict(fit_gtb, n.trees = best_iter, 
                                     newdata = dat_test))

    ## minimal model
    gbt_min <- gbm3::gbmt(formula(paste("y~", paste(var4, collapse = "+"))),
                          data = dat,
                          distribution = gbm_dist("Gaussian"),
                          train_params = 
                            training_params(
                              num_trees = 10000,                   # default: 2000
                              interaction_depth = 3,               # default: 3
                              min_num_obs_in_node = 10,            # default: 10
                              shrinkage = 0.01,                    # default: 0.001
                              bag_fraction = 0.5,                  # default: 0.5
                              num_train = round(0.5 * nrow(dat)),  # default(0.5n) 
                              num_features = length(var4)),        # default: all variables
                 cv_folds = 5,
                 par_details = gbmParallel(num_threads = 1)) # perform CV not in parallel
    ## predict test data using model at stopping iteration
    pred_min <- as.numeric(predict(gbt_min, newdata = dat_test,
                                   n.trees = gbm3::gbmt_performance(gbt_min, method = "cv")))
    
  
    ## Record stop time and determine difference
    stopTime <- Sys.time()
    duration <- difftime(stopTime, startTime, units = "mins")
    ## return final list
    return(list("methName" = 'gbt', "converged" = converged,
                "selMat" = selMat, "varImp" = varImp,
                "predFinal" = pred_final,
                "predMin" = pred_min,
                "duration" = duration,
                "mStop" = as.numeric(best_iter)))
    
    
  } else{    ### if an error occurred (method failed): 
    converged <- FALSE
    ## save corresponding error message
    if(any(class(A) == "error")) {
      error_message <- A$message
    } else if(length(selection) < 1) {
      error_message <- "no variable selected"
    } else {
      error_message <- "any other error"
    }
    ## if the method failed, only the most important predictor x5 is taken
    ## predict with model including x5 only as specified in protocol
    gbt_x5 <- gbm3::gbmt(y ~ x5,
                         data = dat,
                         distribution = gbm_dist("Gaussian"),
                         train_params = 
                           training_params(
                             num_trees = 10000,                   
                             interaction_depth = 1,                # only 1 variable   
                             min_num_obs_in_node = 10,            
                             shrinkage = 0.01,                    
                             bag_fraction = 0.5,                  
                             num_train = round(0.5 * nrow(dat)),  
                             num_features = 1),       
                         cv_folds = 5,
                         par_details = gbmParallel(num_threads = 1)) # perform CV not in parallel
    ## calculate optimal stopping iteration using 5-fold CV
    ## predict test data using model at stopping iteration
    pred_final <- as.numeric(predict(gbt_x5, newdata = dat_test,
                                     n.trees = gbm3::gbmt_performance(gbt_x5, method = "cv")))
    
    ## all other measures are set to NA
    selMat <- NA
    varImp <- NA
    pred_min <- NA
    duration <- NA
    
    return(list("methName" = 'gbt', 
                "converged" = converged, "errorMessage" = error_message,
                "selMat" = selMat, "varImp" = varImp,
                "predFinal" = pred_final,
                "predMin" = pred_min,
                "duration" = duration))
  }
}




#==============================================================================#
#### MFP #######################################################################
#==============================================================================#

apply_mfp <- function(dat, dat_test){
  # Record start time
  startTime <- Sys.time()
  
  ## coding of categorical data 
  dat <- dat |>
    mutate(x4 = factor(x4, ordered = TRUE, levels = c(0, 1, 2)),
           x9 = factor(x9, ordered = FALSE, levels = c(0, 1, 2))) |>
    ## ordinal variable is dummy ordinal coded
    create_dummy_variables(drop_variables = TRUE,
                           var_ordinal = "x4", var_nominal = "x9") |>
    mutate(across(c(x4_1, x4_2), as.numeric)) |>
    relocate(c(x4_1, x4_2), .after = x3) |>
    relocate(c(x91, x92), .after = x8) |>
    rename(x9_1 = x91,
           x9_2 = x92)

  dat_test <- dat_test |>
    mutate(x4 = factor(x4, ordered = TRUE, levels = c(0, 1, 2)),
           x9 = factor(x9, ordered = FALSE, levels = c(0, 1, 2))) |>
    ## ordinal variable is dummy ordinal coded
    create_dummy_variables(drop_variables = TRUE,
                           var_ordinal = "x4", var_nominal = "x9") |>
    mutate(across(c(x4_1, x4_2), as.numeric)) |>
    relocate(c(x4_1, x4_2), .after = x3) |>
    relocate(c(x91, x92), .after = x8) |>
    rename(x9_1 = x91,
           x9_2 = x92)
  
  x_names <- colnames(dat)[colnames(dat) != "y"]
  
  #### variable selection
  ## track if an error occurs
  A <- tryCatch({
    fit_mfp <- mfp2::mfp2(y ~ ., data = dat,
                           center = TRUE, criterion = "aic", 
                           verbose = FALSE, xorder = "ascending")
    
    ## collect selected variables and remove intercept
    cfn <- names(fit_mfp[["coefficients"]])
    cfn <- cfn[cfn != "(Intercept)"]
    ## remove the dot indications (of the power functions) to get the original variables
    cfn_c <- unique(sub("\\..*", "", cfn))
    selection <- cfn_c[order(match(cfn_c, x_names))]
    selMat <- data.frame(matrix(data = ifelse(x_names %in% selection, 1,0), 
                                ncol = length(x_names)))
    names(selMat) <- x_names
  },
  error = function(cond) {
    ## if error - set selection to NULL and return error object
    selection <<- NULL
    return(cond)
  } 
  )

  ## check if an error occurred:
  if(!any(class(A) == "error")) {
    ### if algorithm converged (no error) execute the following code
    converged <- TRUE
    
    #### variable importance
    ## calculate change in AIC when variable is removed
    ## first, calculate AIC of model with selected variables
    fit_mfp_sel <- mfp2::mfp2(formula(paste("y~", paste(selection, collapse = "+"))),
                              center = TRUE, criterion = "aic",
                              verbose = FALSE, xorder = "ascending",
                              keep = selection,
                              data = dat)
    
    sel_AIC <- AIC(fit_mfp_sel)
    ### calculate AIC of model when variable i is dropped
    ## if more than 1 variable was selected
    if(length(selection)>1) {
      ls_AIC <- lapply(seq_along(selection), function(i) {
        ## drop variable
        vars_left <- selection[-i]
        ## refit mfp model without variable i
        fit_mfp_i <- mfp2::mfp2(formula(paste("y~", paste(vars_left, collapse = "+"))),
                                center = TRUE, verbose = FALSE, xorder = "ascending",
                                criterion = "aic", keep = vars_left,
                                data = dat)
        ## calculate change in AIC: greater values implies more importance
        (sel_AIC-AIC(fit_mfp_i)) *-1
      })
      ls_AIC <- setNames(ls_AIC, selection)
      varImp <- as.vector(do.call(cbind, ls_AIC))
      names(varImp) <- selection
    } else {
      varImp <- rep(1, length(selection))
      names(varImp) <- selection
    }

  
    ## get the 4 most important predictors
    var4 <- names(sort(varImp, decreasing = TRUE)[1:4])
    ## only take selected variables if less than 4 variables were selected
    var4 <- na.omit(var4)
    
    #### prediction
    ## final developed model
    pred_final <- as.numeric(predict(fit_mfp, dat_test))
    
    ## minimal model
    fit_mfp_min <- mfp2::mfp2(formula(paste("y~", paste(var4, collapse = " + "))), data = dat,
                              center = TRUE, verbose = FALSE, xorder = "ascending",
                              criterion = "aic", keep = var4)
    pred_min <- as.numeric(predict(fit_mfp_min, dat_test))

    
    ## Record stop time and determine difference
    stopTime <- Sys.time()
    duration <- difftime(stopTime, startTime, units = "mins")
    ## return final list
    return(list("methName" = 'mfp', "converged" = converged,
                "selMat" = selMat, "varImp" = varImp,
                "predFinal" = pred_final,
                "predMin" = pred_min,
                "duration" = duration,
                "fp" = fit_mfp[["fp_terms"]]))
    
    
  } else{    ### if an error occurred (method failed): 
    converged <- FALSE
    ## if the method failed, only the most important predictor x5 is taken
    ## predict with model including x5 only as specified in protocol
    mfp_x5 <- mfp2::mfp2(y ~ x5, data = dat,
                         center = TRUE, criterion = "aic", 
                         verbose = FALSE)
    pred_final <- as.numeric(predict(mfp_x5, dat_test))
    ## all other measures are set to NA
    selMat <- NA
    varImp <- NA
    pred_min <- NA
    duration <- NA
    
    return(list("methName" = 'mfp', 
                "converged" = converged, "errorMessage" = A$message,
                "selMat" = selMat, "varImp" = varImp,
                "predFinal" = pred_final,
                "predMin" = pred_min,
                "duration" = duration))
  }
}




#==============================================================================#
#### Oracle Models of methods ##################################################
#==============================================================================#
apply_oracles <- function(dat, dat_test, truePred){
  ## define predictor variables
  predictors <- names(truePred)[which(truePred[1, ] == 1)]
  ## add x4_2 due to ordinal coding for regression-based methods
  predictors_reg <- ifelse(predictors == "x4", "x4_1", predictors)
  predictors_reg <- c(predictors_reg, "x4_2")[c(1,2,3,9,4,5,6,7,8)]
  
  ### data preparation for regression-based methods
  ## train data
  dat_reg <- dat |> 
    mutate(x4 = factor(x4, ordered = TRUE, levels = c(0, 1, 2))) |> 
    create_dummy_variables(drop_variables = TRUE, var_ordinal = "x4") |>
    mutate(across(c(x4_1, x4_2), as.numeric)) |> 
    select(any_of(predictors_reg), y)
  ## test data
  dat_reg_test <- dat_test |> 
    mutate(x4 = factor(x4, ordered = TRUE, levels = c(0, 1, 2))) |> 
    create_dummy_variables(drop_variables = TRUE, var_ordinal = "x4") |>
    mutate(across(c(x4_1, x4_2), as.numeric)) |> 
    select(any_of(predictors_reg), y)
  
  
  ### data preparation for tree-based methods
  ## train data
  dat_tree <- dat |> 
    mutate(across(.cols = c(x2, x8, x9, x14, x15), as.factor)) |> 
    select(any_of(predictors), y)
  ## test data
  dat_tree_test <- dat_test |> 
    select(any_of(predictors), y)
  
  
  #### LM 
  lm_oracle <- lm(y ~ ., data = dat_reg)
  preds_lm <- as.numeric(predict(lm_oracle, dat_reg_test))
  
  
  ### ENET
  X <- as.matrix(dat_reg[names(dat_reg) != "y"])
  Y <- as.numeric(dat_reg$y)
  X_test <- as.matrix(dat_reg_test[names(dat_reg_test) != "y"])
  ## 5-fold CV of enet (ridge to have all variables in model)
  enet_oracle <- glmnet::cv.glmnet(X, Y, alpha = 0, 
                                   nfolds = 5, type.measure = "mse")
  ## predict test data
  preds_enet <- as.numeric(predict(enet_oracle, s = "lambda.min", 
                                   newx = X_test, type = "response"))
  
  
  ### RF
  ## determine best mtry to fit oracle random forest based on predictor variables
  mtry_seq <- 1:length(predictors)
  list_rf <- lapply(mtry_seq, function(m) {
    ranger::ranger(y ~ ., data = dat_tree, num.trees = 1000, mtry = m,
                   importance = "permutation", num.threads = 1)
  })
  ## extract and identify mtry with lowest OOB error
  oob_error <- unlist(lapply(list_rf, function(x) x$prediction.error))
  k <- which.min(oob_error)
  ## take best model
  rf_oracle <- list_rf[[k]]
  ## predict test data
  preds_rf <- predict(rf_oracle, dat_tree_test)[["predictions"]]
  
  
  ### GBM
  fit_gbm <- mboost::glmboost(y ~ . , data = dat_reg,
                              family = Gaussian(),
                              control = boost_control(mstop = 10000,
                                                      nu = 0.01,
                                                      risk = "inbag"))
  ## calculate optimal stopping iteration using 5-fold CV
  cv <- mboost::cvrisk(fit_gbm, 
                       folds = mboost::cv(model.weights(fit_gbm), 
                                          type = "kfold", B = 5),
                       papply = lapply) # perform CV not in parallel
  ## take model at optimal iteration 
  m_stop <- mboost::mstop(cv)
  gbm_oracle <- fit_gbm[m_stop]
  ## predict test data
  preds_gbm <- as.numeric(predict(gbm_oracle, dat_reg_test))
  
  
  ### GBT
  fit_gtb <- gbm3::gbmt(y ~ .,
                        data = dat_tree,
                        distribution = gbm_dist("Gaussian"),
                        train_params = 
                          training_params(
                            num_trees = 10000,                  
                            interaction_depth = 3,               
                            min_num_obs_in_node = 10,            
                            shrinkage = 0.01,                    
                            bag_fraction = 0.5,                 
                            num_train = round(0.5 * nrow(dat_tree)),  
                            num_features = ncol(dat_tree) - 1),       
                        cv_folds = 5,
                        par_details = gbmParallel(num_threads = 1)) # perform CV not in parallel
  ## calculate optimal stopping iteration using 5-fold CV
  best_iter <- gbm3::gbmt_performance(fit_gtb, method = "cv")
  ## predict test data
  preds_gbt <- as.numeric(predict(fit_gtb, n.trees = best_iter, 
                                  newdata = dat_tree_test))
  
  
  ### MFP
  fit_mfp <- mfp2::mfp2(y ~ ., data = dat_reg,
                        center = TRUE, criterion = "aic", keep = predictors_reg,
                        verbose = FALSE, xorder = "ascending")
  
  ## predict test data
  preds_mfp <- as.numeric(predict(fit_mfp, dat_reg_test))
  
  
  ## return predicted values
  return(list("lm"   = preds_lm,
              "enet" = preds_enet,
              "rf"   = preds_rf,
              "gbm"  = preds_gbm,
              "gbt"  = preds_gbt,
              "mfp"  = preds_mfp))
}

# res <- apply_oracles(dat, dat_test, truePred)



#==============================================================================#
#### True Oracle Model #########################################################
#==============================================================================#

## load true formulas
source("2_Analyze/formulas_true_mods.R")

apply_true_oracle <- function(dat, dat_test, true_formula, formula_all_vars){

  ## fit oracle model with true functional forms of the predictors
  fit_lm <- lm(true_formula, data = dat)
  
  
  #### variable importance
  ## extract predictors with true functional form
  true_preds <- names(coef(fit_lm))[-1]
  ## calculate R^2 of oracle model
  R2 <- summary(fit_lm)$r.square
  ## calculate delta R^2 of model when variable i is dropped
  ls_R2 <- lapply(seq_along(true_preds), function(i) {
    ## drop variable
    preds_left <- true_preds[-i]
    ## refit mfp model without variable i
    fit_lm_i <- lm(formula(paste("y ~", paste(preds_left, collapse = " + "))),
                   data = dat)
    ## calculate change in R^2: greater values implies more importance
    R2 - summary(fit_lm_i)$r.square
  })
  varImp <- unlist(ls_R2)
  
  ## extract variable names from the formulas to get names for varImp
  p_names <- sapply(regmatches(true_preds, gregexpr("x\\d+", true_preds)), unique)
  ## replace x4 by x4_1 to make it consistent with all other methods
  p_names_dummy <- ifelse(p_names == "x4", "x4_1", p_names)
  names(varImp) <- p_names_dummy
  ## get the 4 most important predictors
  var4 <- names(sort(varImp, decreasing = TRUE)[1:4])
  ## only take selected variables if less than 4 variables were selected 
  var4 <- na.omit(var4)
  
  #### prediction
  ### final model
  ## predict test data using final and minimal modl
  ## the final model is the true oracle model
  pred_final <- as.numeric(predict(fit_lm, dat_test))
  
  
  ### minimal model
  ## create formula with true functional forms of 4 most important preds
  ## extract true functions, store as named list, and select those included in var4 
  terms <- as.list(attr(terms(true_formula), "term.labels"))
  names(terms) <- p_names
  terms_var4 <- terms[p_names %in% var4]
  ## create formula and fit lm
  formula_var4 <- formula(paste("y ~", paste(unlist(terms_var4), collapse = " + ")))
  fit_min <- lm(formula_var4, dat)
  pred_min <- as.numeric(predict(fit_min, dat_test))

  ## return final list
  return(list("methName" = 'trueOracle', 
              "varImp" = varImp,
              "predFinal" = pred_final, #"modOracle" = fit_lm, 
              "predMin" = pred_min))#, "modMin" = fit_min))
}










