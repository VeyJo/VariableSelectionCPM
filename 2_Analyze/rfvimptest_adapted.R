##' Testing the statistical significance of predictors in random forests using sequential permutation testing
##'
##' Implements several strategies for testing the statistical significance of predictors in random forests using sequential permutation testing procedures based on the permutation variable importance measure.
##' See Hapfelmeier et al. (2022) for details.
##'
##' Only the general permutation test (\code{test="general"}) controls the type I error. In contrast, the two-sample permutation test (\code{test="twosample"})
##' is associated with inflated type I error, which can lead to false positive findings. An advantage of the two-sample permutation test is that it is
##' very fast. Therefore, this experimental approach may be used as an informal screening tool for finding informative variables.
##' It is, however, not a valid testing procedure. Note also that
##' the paper of Coleman et al. (2019) on which the two-sample test is based has not yet been published in a peer-reviewed journal and that
##' the theory underlying this procedure might thus still need further review.
##'
##' SRPT (\code{type="SRPT"}) and SAPT (\code{type="SAPT"}) are similar sequential procedures, where SRPT is faster with respect to accepting H0, that is, detecting non-informative variables,
##' whereas SAPT is faster with respect to accepting H1, that is, detecting informative variables. Therefore, SRPT may be preferred for
##' datasets with only few informative variables, whereas SAPT is preferable for datasets with many informative variables.
##' The Monte Carlo p-value based testing procedure (\code{type="pval"}) should be used, when p-values are required.
##' The choice \code{type="complete"} offers a conventional permutation test (that is, without sequential testing) (Hapfelmeier and Ulm, 2013). This choice
##' is computationally the most intensive. Lastly, the choice \code{type="certain"} is similar to \code{type="complete"}, but performs
##' early stopping by ending the permutation iterations as soon as it is certain which outcome the conventional permutation test would
##' take. That is, \code{type="certain"} can be considered as a computationally more effective version of \code{type="complete"}.
##'
##' @param data A \code{data.frame} containing the variables in the model.
##' @param yname Name of outcome variable.
##' @param Mmax Maximum number of permutations used in each permutation test. Default is 500.
##' @param varnames Optional. Names of the variables for which testing should be performed. By default all variables in \code{data} with the exception of the outcome variable are used.
##' @param p0 The value of the p-value in the null hypothesis (H0: p = p0) of SPRT and SAPT. Default is 0.06.
##' @param p1 The value of the p-value in the alternative hypothesis (H1: p = p1) of SPRT and SAPT. Default is 0.04.
##' @param alpha The significance level of SPRT when p = p0. Also known as type I error. Default is 0.05.
##' @param beta One minus the power of SPRT when p = p1. Also known as type II error. Default is 0.2.
##' @param A The quantity A in the formula of SAPT. Default is 0.1 for a type I error of 0.05. Usually not changed by the user.
##' @param B The quantity B in the formula of SAPT. Default is 10 (1/A) for a type I error of 0.05. Usually not changed by the user.
##' @param h The quantity h in the formula for the sequential Monte Carlo p-value. The default value for h is 8. Larger values lead to more precise p-value estimates,
##' but are computationally more expensive.
##' @param nperm The numbers of permutations of the out-of-bag observations over which the results are averaged, when calculating the variable importance measure values. Default is 1. Larger values than 1 can only be considered when \code{condinf=TRUE}, that is, when using random forests
##' with conditional inference trees (Hothorn et al., 2006) as base learners.
##' @param ntree Number of trees per forest. Default is 500.
##' @param progressbar Output the current progress of the calculations for each variable to the console? Default is TRUE.
##' @param test Type of the permutation test to perform. This can be either "general" or "twosample", where "general" refers to the usual (sequential) permutation
##' test and "twosample" refers to the two-sample (sequential) permutation test. For the latter, see also Coleman et al. (2019).
##' Note, however, that "twosample" is experimental and should not be used for formal testing. See the details section below.
##' @param type Type of the sequential method to use in the permutation tests. The choices are: "SPRT", "SAPT", "pval", "certain", and "complete". See the 'Details' section below for details.
##' @param condinf Set this value to \code{TRUE} if random forests using conditional inference trees (Hothorn et al., 2006) should
##' be used and to \code{FALSE} if classical random forests using CART trees should be used. Default is \code{FALSE}.
##' @param ... Further arguments passed to \code{ranger::ranger} (if \code{condinf=FALSE}) or \cr \code{party::cforest_unbiased()} (if \code{condinf=TRUE}).
##' @return Object of class \code{rfvimptest} with elements
##'   \item{\code{testtype}}{Type of the permutation test performed and sequential method used.}
##'   \item{\code{varimp}}{Variable importance for each considered independent variable.}
##'   \item{\code{testres}}{The results ("keep H0" vs. "accept H1") of the tests for each considered independent variable.}
##'   \item{\code{pvalues}}{The p-values of the tests for each considered independent variable. Note that p-values are only obtained for the
##'   method types "pval" and "complete".}
##'   \item{\code{stoppedearly}}{For each independent variable, whether the calculations stopped early ("yes") or the maximum of \code{Mmax} permutations was reached ("no").}
##'   \item{\code{perms}}{The number of permutations performed for each independent variable.}
##'   \item{\code{Mmax}}{Maximum number of permutations used in each permutation test.}
##'   \item{\code{ntree}}{Number of trees per forest.}
##'   \item{\code{comptime}}{The time the computations needed.}
##' @examples
##' \donttest{
##'
##' ## Load package:
##' library("rfvimptest")
##'
##' ## Set seed to obtain reproducible results:
##' set.seed(1234)
##'
##' # Load example data:
##' data(hearth2)
##'
##' # NOTE: For illustration purposes a very small number (Mmax=20) of maximum
##' # permutations is considered. This number would be much too small for actual
##' # applications. The default number is Max=500.
##'
##' # By default, SPRT is performed:
##' (ptest_sprt <- rfvimptest(data=hearth2, yname="Class", Mmax=20))
##' ptest_sprt$varimp
##' ptest_sprt$testres
##'
##' # Calculation of p-values using the Monte Carlo p-value based testing procedure:
##' (ptest_pval <- rfvimptest(data=hearth2, yname="Class", type="pval", Mmax=20))
##' ptest_pval$pvalues
##'
##' # If the frequency of informative variables is expected to be high SAPT can be used:
##' (ptest_sapt <- rfvimptest(data=hearth2, yname="Class", type="SAPT", Mmax=20))
##' ptest_sapt$testres
##'
##'
##' # If it is only of interest to test specific variables in the dataset these variables
##' # should be passed to rfvimptest() vias the argument 'varnames' because this
##' # reduces the computational burden considerably:
##'
##' (ptest_twovar <- rfvimptest(data=hearth2, yname="Class", varnames=c("age", "sex"), Mmax=20))
##' ptest_twovar$varimp
##' ptest_twovar$testres
##'
##'
##' # Two-sample permutation test procedures:
##'
##' # NOTE: These should be used only for informal screening for informative variables.
##' # They are not valid statistical tests.
##'
##' # Here, the maximum number of permutations can be much higher because it is necessary
##' # here to construct a new forest for each permutation:
##' rfvimptest(data=hearth2, yname="Class", test="twosample", condinf=TRUE, Mmax=1000)
##'
##' rfvimptest(data=hearth2, yname="Class", test="twosample", type="pval", condinf=TRUE, Mmax=1000)
##'
##' rfvimptest(data=hearth2, yname="Class", test="twosample", type="SAPT", condinf=TRUE, Mmax=1000)
##'
##' }
##'
##' @author Alexander Hapfelmeier, Roman Hornung
##' @references
##' \itemize{
##'   \item Breiman, L. (2001). Random forests. Mach Learn, 45:5-32, <\doi{10.1023/A:1010933404324}>.
##'   \item Coleman, T., Peng, W., Mentch, L. (2019). Scalable and efficient hypothesis testing with random forests. arXiv preprint arXiv:1904.07830, <\doi{10.48550/arXiv.1904.07830}>.
##'   \item Hapfelmeier, A., Hornung, R., Haller, B. (2022). Sequential Permutation Testing of Random Forest Variable Importance Measures. arXiv preprint arXiv:2206.01284, <\doi{10.48550/arXiv.2206.01284}>.
##'   \item Hapfelmeier, A., Ulm, K. (2013). A new variable selection approach using Random Forests. CSDA 60:50–69, <\doi{10.1016/j.csda.2012.09.020}>.
##'   \item Hapfelmeier, A., Hothorn, T., Ulm, K., Strobl, C. (2014). A new variable importance measure for random forests with missing data. Stat Comput 24:21–34, <\doi{10.1007/s11222-012-9349-1}>.
##'   \item Hothorn, T., Hornik, K., Zeileis, A. (2006). Unbiased Recursive Partitioning: A Conditional Inference Framework. J Comput Graph Stat 15(3):651–674, <\doi{10.1198/106186006X133933}>.
##'   \item Wright, M. N., Ziegler, A. (2017). ranger: A fast implementation of random forests for high dimensional data in C++ and R. J Stat Softw 77:1-17, <\doi{10.18637/jss.v077.i01}>.
##'   }
##' @encoding UTF-8
##' @importFrom stats as.formula
##' @importFrom utils setTxtProgressBar txtProgressBar
##' @export
rfvimptest <- function(data, yname, Mmax = 500, varnames = NULL, p0 = 0.06, p1 = 0.04, alpha = 0.05, beta = 0.2, A = 0.1, B = 10, h = 8, nperm = 1,
                       ntree = 500,
                       progressbar = TRUE,
                       test = c("general", "twosample")[1],
                       type = c("SPRT", "SAPT", "pval", "certain", "complete")[1], condinf=FALSE, ...) {
  starttime <- Sys.time()
  # @seealso \code{\link{predict.divfor}}

  if(any(is.na(data))) {
    missvariables <- paste(names(data)[apply(data, 2, function(x) any(is.na(x)))], collapse = ", ")
    stop(paste0("Missing data in columns: ", missvariables, ". Please provide complete data or consider setting condinf=TRUE."))
  }

  if(!condinf & test == "twosample")
    stop("'twosample' approach only useable if condinf=TRUE.")

  if(!condinf & nperm > 1)
    stop("Values of 'nperm' different than 1 are only useable if condinf=TRUE.")

  if (progressbar) pb <- txtProgressBar(min = 1, max = Mmax, initial = 1, width = 10, style = 3, char = "|")

  if (type == "SPRT") {
    A <- beta / (1 - alpha)
    B <- (1 - beta) / alpha
  }
  if (type %in% c("SPRT", "SAPT")) {
    logA <- log(A)
    logB <- log(B)
    help1 <- log((1 - p0) / (1 - p1))
    help2 <- log((p1 * (1 - p0)) / (p0 * (1 - p1)))
  }

  stop_crits <- switch(type,
                       SPRT = list((logA + 1:Mmax * help1) / help2, (logB + 1:Mmax * help1) / help2),
                       SAPT = list((logA + 1:Mmax * help1) / help2, (logB + 1:Mmax * help1) / help2),
                       pval = list(rep(h, times = Mmax), rep(h, times = Mmax)),
                       certain = list(rep(alpha*Mmax, times = Mmax), Mmax*alpha - Mmax + 1:Mmax),
                       complete = NULL)

  if (is.null(varnames))
    varnames <- names(data)[names(data)!=yname]

  if (!condinf) {
    rfmod <- ranger::ranger(data = data, dependent.variable.name=yname, num.tree=ntree, importance = "permutation", ...)
    vimp_orig <- list()
    vimp_orig$values <- rfmod$variable.importance[varnames]
  }
  else {
    rfmod <- party::cforest(as.formula(paste(yname, " ~ .", sep="")), data = data, controls = party::cforest_unbiased(ntree=ntree, ...))
    vimp_orig <- permimp::permimp(rfmod, whichxnames = varnames, nperm = nperm, asParty = TRUE, progressBar = FALSE)
  }


  if (test == "general") {
    testresult <-
      lapply(varnames, function(v) {
        permdata <- data
        permvimps <- c()
        for (m in 1:Mmax) {
          if (progressbar) {setTxtProgressBar(pb, m)}
          if (m == 1 & progressbar) cat(" of variable", v)
          permdata[, v] <- sample(permdata[, v])

          if (!condinf) {
            permmod <- ranger::ranger(data = permdata, dependent.variable.name=yname, num.tree=ntree, importance="permutation", ...)
            permvimps <- c(permvimps, permmod$variable.importance[v])
          }
          else {
            permmod <- party::cforest(as.formula(paste(yname, " ~ .", sep="")), data = permdata, controls = party::cforest_unbiased(ntree=ntree, ...))
            permvimps <- c(permvimps, permimp::permimp(permmod, whichxnames = v, nperm = nperm, asParty = TRUE, progressBar = FALSE)$values)
          }
          d <- sum(permvimps >= vimp_orig$values[v])
          if (type == "certain") {
            if (d > stop_crits[[1]][m]) {result <- "keep H0"; pvalue <- NA; break} else if (d <= stop_crits[[2]][m]) {result <- "accept H1"; pvalue <- NA; break}
          } else if (type %in% c("SPRT", "SAPT")) {
            if (d >= stop_crits[[1]][m]) {result <- "keep H0"; pvalue <- NA; break} else if (d <= stop_crits[[2]][m]) {result <- "accept H1"; pvalue <- NA; break}
          } else if (type == "pval") {
            if (d == h) {
              pvalue <- d/m
              result <- ifelse(pvalue > 0.05, "keep H0", "accept H1")
              break
            }
          }
        }
        if (m == Mmax) {
          if (type == "pval") {
            if (d < h) {
              pvalue <- (d + 1) / (Mmax + 1)
              result <- ifelse(pvalue > 0.05, "keep H0", "accept H1")
            }
            else {
              pvalue <- d / Mmax
              result <- ifelse(pvalue > 0.05, "keep H0", "accept H1")
            }
          } else  if (type == "complete") {
            pvalue <- d / Mmax
            result <- ifelse(pvalue > 0.05, "keep H0", "accept H1")
          } else {
            pvalue <- NA
            result <- ifelse(d / Mmax > 0.05, "keep H0", "accept H1")
          }
        }
        if (progressbar) cat(" - finished", "\n")
        list(testres=result, pvalue=pvalue, stoppedearly=ifelse(m < Mmax, "yes", "no"), permperf=m)
      })
  } else if (test == "twosample") {
    testresult <-
      lapply(varnames, function(v) {
        permdata <- data
        permdata[, v] <- sample(permdata[, v])
        permmod <- party::cforest(as.formula(paste(yname, " ~ .", sep="")), data = permdata, controls = party::cforest_unbiased(ntree=ntree, ...))
        permmodvimps <- permimp::permimp(permmod, whichxnames = v, nperm = nperm, asParty = TRUE, progressBar = FALSE)
        permvimps <- c()
        for (m in 1:Mmax) {
          if (progressbar) {setTxtProgressBar(pb, m)}
          if (m == 1 & progressbar) cat(" of variable", v)
          permvimps <- c(permvimps, mean(c(vimp_orig$perTree[, v], permmodvimps$perTree[, v])[sample(1:(2*ntree), ntree)]))
          d <- sum(permvimps >= vimp_orig$values[v])
          if (type == "certain") {
            if (d > stop_crits[[1]][m]) {result <- "keep H0"; pvalue <- NA; break} else if (d <= stop_crits[[2]][m]) {result <- "accept H1"; pvalue <- NA; break}
          } else if (type %in% c("SPRT", "SAPT")) {
            if (d >= stop_crits[[1]][m]) {result <- "keep H0"; pvalue <- NA; break} else if (d <= stop_crits[[2]][m]) {result <- "accept H1"; pvalue <- NA; break}
          } else if (type == "pval") {
            if (d == h) {
              pvalue <- d/m
              result <- ifelse(pvalue > 0.05, "keep H0", "accept H1")
              break
            }
          }
        }
        if (m == Mmax) {
          if (type == "pval") {
            if (d < h) {
              pvalue <- (d + 1) / (Mmax + 1)
              result <- ifelse(pvalue > 0.05, "keep H0", "accept H1")
            }
            else {
              pvalue <- d / Mmax
              result <- ifelse(pvalue > 0.05, "keep H0", "accept H1")
            }
          } else  if (type == "complete") {
            pvalue <- d / Mmax
            result <- ifelse(pvalue > 0.05, "keep H0", "accept H1")
          } else {
            pvalue <- NA
            result <- ifelse(d / Mmax > 0.05, "keep H0", "accept H1")
          }
        }
        if (progressbar) cat(" - finished", "\n")
        list(testres=result, pvalue=pvalue, stoppedearly=ifelse(m < Mmax, "yes", "no"), permperf=m)
      })
  }
  stoptime <- Sys.time()
  time_elapsed <- paste0(round(as.numeric(difftime(stoptime, starttime, units = "secs")), 1), " seconds")
  testres <- sapply(testresult, function(x) x$testres)
  pvalues <- sapply(testresult, function(x) x$pvalue)
  stoppedearly <- sapply(testresult, function(x) x$stoppedearly)
  perms <- sapply(testresult, function(x) x$permperf)
  names(testres) <- names(pvalues) <- names(stoppedearly) <- names(perms) <- names(vimp_orig$values)

  result <- list(testtype = paste(test, type, sep=", "), varimp=vimp_orig$values, testres = testres, pvalues = pvalues,
       stoppedearly = stoppedearly, perms = perms, Mmax=Mmax, ntree=ntree, comptime = time_elapsed)

  class(result) <- "rfvimptest"

  return(result)

}
