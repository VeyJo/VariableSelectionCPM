##----------------------------------------------------------------------------##
##  functions to calculate performance measures
##----------------------------------------------------------------------------##

# variable selection ####

## TIF #####
## True Inclusion Frequency
fun_TIF <- function(selMat, truePred) { # truePred: true predictors; selMat: selected variables
  ## combine dummy variables to one variable if dummies exist
  if (length(selMat) > 15) {
    ## state variable as selected if one of the dummy variables was selected
    selMat$x4 <- pmax(selMat$x4_1, selMat$x4_2)
    selMat$x9 <- pmax(selMat$x9_1, selMat$x9_2)
    selMat <- selMat[, !names(selMat) %in% c("x4_1", "x4_2", "x9_1", "x_9_2")]
    ## reorder selMat
    selMat <- selMat[, paste0("x", 1:15)]
  }
  
  TI <- vector()
  for (i in 1:length(truePred)) {
    if (truePred[1, i] == 1 & selMat[1, i] == 1) {
      TI[i] <- 1 
    } else {
      TI[i] <- 0
    }
  }
  # calculate proportion of correctly selected true variables 
  TIF <- sum(TI) / sum(truePred == 1)
  TIF
}
# fun_TIF(selMat, truePred)



## TEF ####
## True Exclusion Frequency
fun_TEF <- function(selMat, truePred) { # truePred: true predictors; selMat: selected variables
  ## combine dummy variables to one variable if dummies exist
  if (length(selMat) > 15) {
    ## state variable as selected if one of the dummy variables was selected
    selMat$x4 <- pmax(selMat$x4_1, selMat$x4_2)
    selMat$x9 <- pmax(selMat$x9_1, selMat$x9_2)
    selMat <- selMat[, !names(selMat) %in% c("x4_1", "x4_2", "x9_1", "x_9_2")]
    ## reorder selMat
    selMat <- selMat[, paste0("x", 1:15)]
  }
  
  TE <- vector()
  for (i in 1:length(truePred)) {
    if (truePred[1, i] == 0 & selMat[1, i] == 0) {
      TE[i] <- 1 
    } else {
      TE[i] <- 0
    }
  }
  # calculate proportion of correctly excluded true variables  
  TEF <- sum(TE) / sum(truePred == 0)
  TEF
}
# fun_TEF(selMat, truePred)


## TCR ####
# True Classification Rate = Accuracy
fun_TCR <- function(selMat, truePred){
  ## combine dummy variables to one variable if dummies exist
  if (length(selMat) > 15) {
    ## state variable as selected if one of the dummy variables was selected
    selMat$x4 <- pmax(selMat$x4_1, selMat$x4_2)
    selMat$x9 <- pmax(selMat$x9_1, selMat$x9_2)
    selMat <- selMat[, !names(selMat) %in% c("x4_1", "x4_2", "x9_1", "x_9_2")]
    ## reorder selMat
    selMat <- selMat[, paste0("x", 1:15)]
  }
  
  TC <- vector()
  for(i in 1:length(truePred)){
    if (truePred[1, i] == selMat[1, i]){
      TC[i] <- 1 
    } else{
      TC[i] <- 0
    }
  }
  # calculate accuracy rate
  TCR <- sum(TC) / ncol(selMat)
  TCR
}


## MSF ####
## Model Selection Frequency
fun_MSF <- function(selMat, truePred) {
  ## combine dummy variables to one variable if dummies exist
  if (length(selMat) > 15) {
    ## state variable as selected if one of the dummy variables was selected
    selMat$x4 <- pmax(selMat$x4_1, selMat$x4_2)
    selMat$x9 <- pmax(selMat$x9_1, selMat$x9_2)
    selMat <- selMat[, !names(selMat) %in% c("x4_1", "x4_2", "x9_1", "x_9_2")]
    ## reorder selMat
    selMat <- selMat[, paste0("x", 1:15)]
  }
  ## because truePred contains integer 
  selMat <- data.frame(lapply(selMat, as.integer))
  TM <- if (identical(truePred, selMat)){
      TM <- 1 
    } else{
      TM <- 0
    }
  TM
}
# fun_MSF(selMat, truePred)



## MS ####
## Model Size as total number of selected variables
fun_MS <- function(selMat) {
  ## calculate number of selected variables
  sum(selMat)
}

## SV ####
## Proportion of Selected Variables
fun_SV <- function(selMat) {
  # calculate % of selected variables, as different numbers of variables (15, 17)
  sum(selMat) / ncol(selMat)
}


# variable importance ####

## tau  & varImp ####
## input is: results <- res[[r]]
fun_vimp <- function(results, trueImp_dR2) {
  m <- results$methName
  varImp <- results$varImp
  trueImp <- trueImp_dR2
  
  ## adapt names of true varimp vector in dependence of coding of categorical variables
  if(m %in% c("stepLin", "mfp", "enet", "gbm", "trueOracle")) {
    ## dummy coding of parametric methods
    trueImp$predictor <- ifelse(trueImp$predictor == "x4a", "x4_1", trueImp$predictor)
    trueImp <- setNames(as.numeric(trueImp[, 2]), trueImp$predictor)
  } else {                                             
    ## dummy coding of non-parametric methods
    trueImp$predictor <- ifelse(trueImp$predictor == "x4a", "x4", trueImp$predictor)
    trueImp <- setNames(as.numeric(trueImp[, 2]), trueImp$predictor)
  }
  
  ## create data.frame with selection status, variable importance and true varimp
  if(m != "trueOracle"){
    tab <- data.frame(t(results$selMat))
    colnames(tab) <- "selection"
    tab[c("varImp", "trueImp")] <- 0
    tab[names(varImp), "varImp"] <- varImp
    tab[names(trueImp), "trueImp"] <- trueImp
  } else if(m == "trueOracle") {
    tab <- data.frame(cbind(varImp, trueImp))
  }
  
  ## add ranks
  tab["rank_m"] <- rank(-tab[, "varImp"], ties.method = "average")
  tab["rank_true"] <- rank(-tab[, "trueImp"], ties.method = "average")
  
  ## calculate Kendall's tau
  tau <- cor(tab$rank_m, tab$rank_true, method = "kendall")
  ## return varimp table and tau
  return(list("imp" = tab, 
              "tau" = tau))
}
# fun_vimp(results, trueImp_dR2)




# prediction ####

## RMSPE ####
fun_RMSPE <- function(preds, obs) {
  sqrt(mean((obs - preds)^2))
}


## Calibration intercept ####
fun_calA <- function(preds, obs) {
  unname(coef(lm(obs ~ offset(preds)))[1])
}

## Calibration slope ####
fun_calB <- function(preds, obs) {
  unname(coef(lm(obs ~ preds))[2])
}


## for MFP ####
## RMSPE 
fun_RMSPE <- function(preds, obs) {
  preds <- ifelse(is.infinite(preds), NA, preds)
  sqrt(mean((obs - preds)^2, na.rm = TRUE))
}

## Calibration intercept
fun_calA <- function(preds, obs) {
  preds <- ifelse(is.infinite(preds), NA, preds)
  unname(coef(lm(obs ~ offset(preds)))[1])
}

## Calibration slope 
fun_calB <- function(preds, obs) {
  preds <- ifelse(is.infinite(preds), NA, preds)
  unname(coef(lm(obs ~ preds))[2])
}





