##----------------------------------------------------------------------------##
##  true formulas for complexity level A-D
##----------------------------------------------------------------------------##


### containing true predictors ####
#### A -------------------------------------------------------------------------
frmla_true_pred_A <- as.formula(
  y ~ 
    I(-0.02 * x1) +
    I(0.7 * ((x3 + 10) / 25)) +
    I(-0.3*(x4==1)) + 
    I(-0.12 * x5) + 
    I(0.0015 * x6) +
    I(0.2 * x8) + 
    I(0.015 * x10) + 
    I(0.02*x11)
)

#### B -------------------------------------------------------------------------
frmla_true_pred_B <- as.formula(
  y ~ 
    I(3.5*x1^0.5-0.28*x1) +
    I(5*(log((x3+10)/27))^2) +
    I(-0.8*(x4==1)) + 
    I(-0.37*x5-1.8*exp(-((log(x5)-1.5)^2)/0.4)) + 
    I(0.42*log(x6+1)) +
    I(0.55*x8) + 
    I(0.04*x10) + 
    I(0.07*x11)
)


#### C -------------------------------------------------------------------------
frmla_true_pred_C <- as.formula(
  y ~
    I(3.5*x1^0.5-0.28*x1) +
    I(5*(log((x3+10)/27))^2) +
    I(-0.8*(x4==1)) + 
    I(-0.37*x5-1.8*exp(-((log(x5)-1.5)^2)/0.4)) + 
    I(0.42*log(x6+1)) +
    I(0.55*x8) + 
    I(0.04*x10) + 
    I(0.07*x11)
)



#### D -------------------------------------------------------------------------
frmla_true_pred_D <- as.formula(
  y ~
    I(3.5*x1^0.5-0.28*x1) +
    I(5*(log((x3+10)/27))^2) +
    I(-0.8*(x4==1)) +
    I(-0.37*x5-1.8*exp(-((log(x5)-1.5)^2)/0.4)) +
    I(0.42*log(x6+1)) +
    I(3*x8) +
    I(0.04*x10) +
    I(ifelse(x8==0, 0.92*log(x11+1), 1-0.04*x11))
)


## store as list
ls_true_frmlas <- list(A = frmla_true_pred_A,
                       B = frmla_true_pred_B,
                       C = frmla_true_pred_C,
                       D = frmla_true_pred_D)

rm(frmla_true_pred_A, frmla_true_pred_B, frmla_true_pred_C, frmla_true_pred_D)

