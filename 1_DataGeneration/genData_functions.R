##----------------------------------------------------------------------------##
##    funcitons to generate data for complexity levels A-D
##----------------------------------------------------------------------------##

library(mvtnorm)


#==============================================================================#
# complexity level A 
#==============================================================================#

generate_compl_A <- function(n, var_e, quantiles){
  ## create correlation matrix
  sigma <- matrix(0, 15, 15)
  diag(sigma) <- 1
  
  sigma[1, 2] <- 0.8
  sigma[1, 9] <- 0.3
  sigma[3, 9] <- -0.5
  sigma[8, 9] <- -0.3
  sigma[3, 5] <- 0.3
  sigma[8, 11] <- 0.3
  sigma[4, 7] <- -0.3
  sigma[6, 9] <- 0.3
  sigma[4, 6] <- -0.5
  sigma[5, 6] <- -0.3
  sigma[5, 12] <- 0.5
  sigma[6, 7] <- 0.5
  sigma[7, 11] <- 0.3
  sigma[6, 11] <- 0.5
  sigma[6, 14] <- 0.3
  sigma[7, 14] <- 0.3
  sigma[11, 14] <- 0.5
  
  for (i in 1:14) {
    for (j in (i + 1):15)
      sigma[j, i] <- sigma[i, j]
  }
  
  # draw from multivariate normal
  z <- rmvnorm(n, mean = rep(0, 15), sigma)
  
  ## transformation of variables
  x <- z * 0
  x[, 1]  <- floor(10 * z[, 1] + 55)
  x[, 2]  <- (z[, 2] < 0.6) * 1
  x[, 3]  <- exp(0.4 * z[, 3] + 3)
  x[, 4]  <- (z[, 4] > -1.2) + (z[, 4] >= 0.75)
  x[, 5]  <- exp(0.5 * z[, 5] + 1.5)
  x[, 6]  <- floor(apply(cbind(100 * exp(z[, 6]) - 20, 0), 1, max))
  x[, 7]  <- floor(apply(cbind(80 * exp(z[, 7]) - 20, 0), 1, max))
  x[, 8]  <- (z[, 8] < -0.35)
  x[, 9]  <- ((z[, 9] >= 0.5) & (z[, 9] < 1.5)) + 2 * (z[, 9] >= 1.5)
  x[, 10] <- 0.01 * (100 * (z[, 10] + 4) ** 2)
  x[, 11] <- floor(10 * z[, 11] + 55)
  x[, 12] <- floor(10 * z[, 12] + 55)
  x[, 13] <- floor(10 * z[, 13] + 55)
  x[, 14] <- (z[, 14] < 0)
  x[, 15] <- (z[, 15] < 0)
  

  ## remove extreme values and draw new value
  for (i in c(1,3,5,6,7,10,11,12,13)) { 
    for(j in 1:nrow(x)) {
      ## upper bound
      while(x[j,i] > quantiles[i, "third_plus_5iqr"]){
        x[j,i] <- rnorm(1, quantiles[i, "75%"], (quantiles[i, "50%"]/5))
      }
      ## lower bound not necessary
    }
  }
  
  ## create true predictors 
  pred <- x * 0
  pred[, 1]  <- -0.02 * x[, 1]
  pred[, 3]  <- 0.7 * ((x[, 3] + 10) / 25)
  pred[, 4]  <- -0.3 * (x[, 4] == 1)
  pred[, 5]  <- -0.12 * x[, 5]
  pred[, 6]  <- 0.0015 * x[, 6]
  pred[, 8]  <- 0.2 * x[, 8]
  pred[, 10] <- 0.015 * x[, 10]
  pred[, 11] <- 0.02 * x[, 11]
  
  ## calculate outcome by adding up predictor effects
  y <- apply(pred, 1, sum)
  ## add random error (noise) to y 
  y = y + rnorm(n, mean = 0, sd = sqrt(var_e))
  
  df <- data.frame(x, y)
  colnames(df)[1:15] <- paste0("x", 1:15)
  
  return(df)
}



#==============================================================================#
# complexity level B 
#==============================================================================#

generate_compl_B <- function(n, var_e, quantiles){
  ## create correlation matrix
  sigma <- matrix(0, 15, 15)
  diag(sigma) <- 1
  
  sigma[1, 2] <- 0.8
  sigma[1, 9] <- 0.3
  sigma[3, 9] <- -0.5
  sigma[8, 9] <- -0.3
  sigma[3, 5] <- 0.3
  sigma[8, 11] <- 0.3
  sigma[4, 7] <- -0.3
  sigma[6, 9] <- 0.3
  sigma[4, 6] <- -0.5
  sigma[5, 6] <- -0.3
  sigma[5, 12] <- 0.5
  sigma[6, 7] <- 0.5
  sigma[7, 11] <- 0.3
  sigma[6, 11] <- 0.5
  sigma[6, 14] <- 0.3
  sigma[7, 14] <- 0.3
  sigma[11, 14] <- 0.5
  
  for (i in 1:14) {
    for (j in (i + 1):15)
      sigma[j, i] <- sigma[i, j]
  }
  
  # draw from multivariate normal
  z <- rmvnorm(n, mean = rep(0, 15), sigma)
  
  ## transformation of variables
  x <- z * 0
  x[, 1]  <- floor(10 * z[, 1] + 55)
  x[, 2]  <- (z[, 2] < 0.6) * 1
  x[, 3]  <- exp(0.4 * z[, 3] + 3)
  x[, 4]  <- (z[, 4] > -1.2) + (z[, 4] >= 0.75)
  x[, 5]  <- exp(0.5 * z[, 5] + 1.5)
  x[, 6]  <- floor(apply(cbind(100 * exp(z[, 6]) - 20, 0), 1, max))
  x[, 7]  <- floor(apply(cbind(80 * exp(z[, 7]) - 20, 0), 1, max))
  x[, 8]  <- (z[, 8] < -0.35)
  x[, 9]  <- ((z[, 9] >= 0.5) & (z[, 9] < 1.5)) + 2 * (z[, 9] >= 1.5)
  x[, 10] <- 0.01 * (100 * (z[, 10] + 4) ** 2)
  x[, 11] <- floor(10 * z[, 11] + 55)
  x[, 12] <- floor(10 * z[, 12] + 55)
  x[, 13] <- floor(10 * z[, 13] + 55)
  x[, 14] <- (z[, 14] < 0)
  x[, 15] <- (z[, 15] < 0)
  
  ## remove extreme values and draw new value
  for (i in c(1,3,5,6,7,10,11,12,13)) { 
    for(j in 1:nrow(x)) {
      ## upper bound
      while(x[j,i] > quantiles[i, "third_plus_5iqr"]){
        x[j,i] <- rnorm(1, quantiles[i, "75%"], (quantiles[i, "50%"]/5))
      }
      ## lower bound not necessary
    }
  }
  
  # create predictor variables
  pred <- x * 0
  pred[, 1] <- 3.5 * x[, 1] ^ 0.5 - 0.28 * x[, 1]
  pred[, 3] <- 5 * (log((x[, 3] + 10) / 27)) ^ 2
  pred[, 4] <- -0.8 * (x[, 4] == 1)
  pred[, 5] <- -0.37 * x[, 5] - 1.8 * exp(-((log(x[, 5]) - 1.5) ^ 2) / 0.4)
  pred[, 6] <- 0.42 * log(x[, 6] + 1)
  pred[, 8] <- 0.55 * x[, 8]
  pred[, 10] <- 0.04 * x[, 10]
  pred[, 11] <- 0.07 * x[, 11]
  
  ## calculate outcome by adding up predictor effects
  y <- apply(pred, 1, sum)
  ## add random error (noise) to y 
  y = y + rnorm(n, mean = 0, sd = sqrt(var_e))
  
  df <- data.frame(x, y)
  colnames(df)[1:15] <- paste0("x", 1:15)
  
  return(df)
} 



#==============================================================================#
# complexity level C 
#==============================================================================#

generate_compl_C <- function(n, var_e, quantiles){
  ## create correlation matrix
  sigma <- matrix(0, 15, 15)
  diag(sigma) <- 1
  
  sigma[1, 2] <- 0.8
  sigma[1, 9] <- 0.3
  sigma[3, 9] <- -0.5
  sigma[8, 9] <- -0.3
  sigma[3, 5] <- 0.3
  sigma[8, 11] <- 0.3
  sigma[4, 7] <- -0.3
  sigma[6, 9] <- 0.3
  sigma[4, 6] <- -0.5
  sigma[5, 6] <- -0.3
  sigma[5, 12] <- 0.5
  sigma[6, 7] <- 0.5
  sigma[7, 11] <- 0.3
  sigma[6, 11] <- 0.5
  sigma[6, 14] <- 0.3
  sigma[7, 14] <- 0.3
  sigma[11, 14] <- 0.5
  
  for (i in 1:14) {
    for (j in (i + 1):15)
      sigma[j, i] <- sigma[i, j]
  }
  
  # draw from multivariate normal
  z <- rmvnorm(n, mean = rep(0, 15), sigma)
  
  ## transformation of variables
  x <- z * 0
  x[, 1]  <- floor(10 * z[, 1] + 55)
  x[, 2]  <- (z[, 2] < 0.6) * 1
  x[, 3]  <- exp(0.4 * z[, 3] + 3)
  x[, 4]  <- (z[, 4] > -1.2) + (z[, 4] >= 0.75)
  x[, 5]  <- exp(0.5 * z[, 5] + 1.5)
  x[, 6]  <- floor(apply(cbind(100 * exp(z[, 6]) - 20, 0), 1, max))
  x[, 7]  <- floor(apply(cbind(80 * exp(z[, 7]) - 20, 0), 1, max))
  x[, 8]  <- (z[, 8] < -0.35)
  x[, 9]  <- ((z[, 9] >= 0.5) & (z[, 9] < 1.5)) + 2 * (z[, 9] >= 1.5)
  x[, 10] <- 0.01 * (100 * (z[, 10] + 4)^2)
  ## x11-x13 generated below
  x[, 14] <- (z[, 14] < 0)
  x[, 15] <- (z[, 15] < 0)
  
  ## remove extreme values and draw new value
  for (i in c(1,3,5,6,7,10)) { 
    for(j in 1:nrow(x)) {
      ## upper bound
      while(x[j,i] > quantiles[i, "third_plus_5iqr"]){
        x[j,i] <- rnorm(1, quantiles[i, "75%"], (quantiles[i, "50%"]/5))
      }
      ## lower bound not necessary
    }
  }
  
  ## generating variables from others to introduce non-linear correlation
  # x11 correlation with: x6 & x14 0.5;   x7 & x8 0.3
  x[, 11] <- 8 + 7*log(x[, 6]+1) + 2*x[, 14] + 0.05*(x[, 7])^0.5 - 2*x[, 8]
  ## x12 correlation with: x5 0.5
  x[, 12] <- 180-0.4*(x[, 5]-7)^2   ## quadratic relation
  ## x13 non-monotonic relation with x3 (mimicking e.g. age-bmi)
  x[, 13] <- 46 - 3.91*x[, 3] + 0.14*x[, 3]^2 - 0.0012*x[, 3]^3
  
  ## adding noise
  x[, 11] <- apply(cbind(x[, 11] + rnorm(n = n, mean = 0, sd = 0.2*median(x[, 11])), 0), 1, max)
  x[, 12] <- x[, 12] + rnorm(n = n, mean = 0, sd = 0.05*median(x[, 12]))
  x[, 13] <- x[, 13] + rnorm(n = n, mean = 0, sd = 0.2*median(x[, 13]))
  
  ## remove extreme values and draw new value
  for (i in c(11,12,13)) {
    for (j in 1:nrow(x)) {
      ## upper bound
      while(x[j,i] > quantiles[i, "third_plus_5iqr"]){
        x[j,i] <- rnorm(1, quantiles[i, "75%"], (quantiles[i, "50%"]/5))
      }
      ## lower bound
      while(x[j,i] < quantiles[i, "first_minus_5iqr"]){
        x[j,i] <- rnorm(1, quantiles[i, "25%"], (quantiles[i, "50%"]/5))
      }
    }
  }
  
  ## create predictor variables
  pred <- x * 0
  pred[, 1] <- 3.5 * x[, 1] ^ 0.5 - 0.28 * x[, 1]
  pred[, 3] <- 5 * (log((x[, 3] + 10) / 27)) ^ 2
  pred[, 4] <- -0.8 * (x[, 4] == 1)
  pred[, 5] <- -0.37 * x[, 5] - 1.8 * exp(-((log(x[, 5]) - 1.5) ^ 2) / 0.4)
  pred[, 6] <- 0.42 * log(x[, 6] + 1)
  pred[, 8] <- 0.55 * x[, 8]
  pred[, 10] <- 0.04 * x[, 10]
  pred[, 11] <- 0.07 * x[, 11]
  
  ## calculate outcome by adding up predictor effects
  y <- apply(pred, 1, sum)
  ## add random error (noise) to y 
  y = y + rnorm(n, mean = 0, sd = sqrt(var_e))
  
  df <- data.frame(x, y)
  colnames(df)[1:15] <- paste0("x", 1:15)
  
  return(df)
} 




#==============================================================================#
# complexity level D
#==============================================================================#

generate_compl_D <- function(n, var_e, quantiles){
  ## create correlation matrix
  sigma <- matrix(0, 15, 15)
  diag(sigma) <- 1
  
  sigma[1, 2] <- 0.8
  sigma[1, 9] <- 0.3
  sigma[3, 9] <- -0.5
  sigma[8, 9] <- -0.3
  sigma[3, 5] <- 0.3
  sigma[8, 11] <- 0.3
  sigma[4, 7] <- -0.3
  sigma[6, 9] <- 0.3
  sigma[4, 6] <- -0.5
  sigma[5, 6] <- -0.3
  sigma[5, 12] <- 0.5
  sigma[6, 7] <- 0.5
  sigma[7, 11] <- 0.3
  sigma[6, 11] <- 0.5
  sigma[6, 14] <- 0.3
  sigma[7, 14] <- 0.3
  sigma[11, 14] <- 0.5
  
  for (i in 1:14) {
    for (j in (i + 1):15)
      sigma[j, i] <- sigma[i, j]
  }
  
  # draw from multivariate normal
  z <- rmvnorm(n, mean = rep(0, 15), sigma)
  
  ## transformation of variables
  x <- z * 0
  x[, 1]  <- floor(10 * z[, 1] + 55)
  x[, 2]  <- (z[, 2] < 0.6) * 1
  x[, 3]  <- exp(0.4 * z[, 3] + 3)
  x[, 4]  <- (z[, 4] > -1.2) + (z[, 4] >= 0.75)
  x[, 5]  <- exp(0.5 * z[, 5] + 1.5)
  x[, 6]  <- floor(apply(cbind(100 * exp(z[, 6]) - 20, 0), 1, max))
  x[, 7]  <- floor(apply(cbind(80 * exp(z[, 7]) - 20, 0), 1, max))
  x[, 8]  <- (z[, 8] < -0.35)
  x[, 9]  <- ((z[, 9] >= 0.5) & (z[, 9] < 1.5)) + 2 * (z[, 9] >= 1.5)
  x[, 10] <- 0.01 * (100 * (z[, 10] + 4)^2)
  ## x11-x13 generated below
  x[, 14] <- (z[, 14] < 0)
  x[, 15] <- (z[, 15] < 0)
  
  ## remove extreme values and draw new value
  for (i in c(1,3,5,6,7,10)) { 
    for(j in 1:nrow(x)) {
      ## upper bound
      while(x[j,i] > quantiles[i, "third_plus_5iqr"]){
        x[j,i] <- rnorm(1, quantiles[i, "75%"], (quantiles[i, "50%"]/5))
      }
      ## lower bound not necessary
    }
  }
  
  ## generating variables from others to introduce non-linear correlation
  # x11 correlation with: x6 & x14 0.5;   x7 & x8 0.3
  x[, 11] <- 8 + 7*log(x[, 6]+1) + 2*x[, 14] + 0.05*(x[, 7])^0.5 - 2*x[, 8]
  ## x12 correlation with: x5 0.5
  x[, 12] <- 180-0.4*(x[, 5]-7)^2   ## quadratic relation
  ## x13 non-monotonic relation with x3 (mimicking e.g. age-bmi)
  x[, 13] <- 46 - 3.91*x[, 3] + 0.14*x[, 3]^2 - 0.0012*x[, 3]^3
  
  ## adding noise
  x[, 11] <- apply(cbind(x[, 11] + rnorm(n = n, mean = 0, sd = 0.2*median(x[, 11])), 0), 1, max)
  x[, 12] <- x[, 12] + rnorm(n = n, mean = 0, sd = 0.05*median(x[, 12]))
  x[, 13] <- x[, 13] + rnorm(n = n, mean = 0, sd = 0.2*median(x[, 13]))

  ## remove extreme values and draw new value
  for (i in c(11,12,13)) {
    for (j in 1:nrow(x)) {
      ## upper bound
      while(x[j,i] > quantiles[i, "third_plus_5iqr"]){
        x[j,i] <- rnorm(1, quantiles[i, "75%"], (quantiles[i, "50%"]/5))
      }
      ## lower bound
      while(x[j,i] < quantiles[i, "first_minus_5iqr"]){
        x[j,i] <- rnorm(1, quantiles[i, "25%"], (quantiles[i, "50%"]/5))
      }
    }
  }
  
  ## create predictor variables
  pred <- x * 0
  pred[, 1] <- 3.5 * x[, 1] ^ 0.5 - 0.28 * x[, 1]
  pred[, 3] <- 5 * (log((x[, 3] + 10) / 27)) ^ 2
  pred[, 4] <- -0.8 * (x[, 4] == 1)
  pred[, 5] <- -0.37 * x[, 5] - 1.8 * exp(-((log(x[, 5]) - 1.5) ^ 2) / 0.4)
  pred[, 6] <- 0.42 * log(x[, 6] + 1)
  pred[, 8] <- 3 * x[, 8]
  pred[, 10] <- 0.04 * x[, 10]
  ## interaction term
  condition_x8 <- x[, 8] == 0
  pred[condition_x8, 11] <- 0.92 * log(x[condition_x8, 11] + 1)
  pred[!condition_x8, 11] <- 1 - 0.04 * x[!condition_x8, 11]
  
  ## calculate outcome by adding up predictor effects
  y <- apply(pred, 1, sum)
  ## add random error (noise) to y 
  y = y + rnorm(n, mean = 0, sd = sqrt(var_e))
  
  df <- data.frame(x, y)
  colnames(df)[1:15] <- paste0("x", 1:15)
  
  return(df)
} 

