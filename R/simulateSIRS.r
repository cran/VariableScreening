#' @title Simulate a dataset for demonstrating the performance of screenIID with the SIRS method
#'
#' @description Simulates a dataset that can be used to demonstrate variable screening
#' for ultrahigh-dimensional regression with the SIRS option in screenIID.
#' The simulated dataset has p numerical predictors X and a categorical Y-response.
#' The data-generating scenario is a simplified version of Example 1 of  Zhu, Li, Li and Zhu (2011).
#' Specifically, the X covariates are normally distributed with mean zero and variance one, and
#' may be correlated if the argument rho is set to a nonzero value.  The response Y is generated as
#'        Y = c*X1 + 0.8*c*X2 + 0.6*c*X3 + 0.4*c*X4 + 0.5*c*X5 + sigma*e.
#' where c is the argument SignalStrength, e is either a standard normal
#' distribution (if HeavyTailedResponse==FALSE) or t distribution with 1 degree of
#' freedom (if HeavyTailedResponse==TRUE). sigma is either sqrt(6.83) if heteroskedastic==FALSE,
#' or else exp(X20+X21+X22) if  heteroskedastic=TRUE.
#'
#' @references Zhu, L.-P., Li, L., Li, R., & Zhu, L.-X. (2011). Model-free feature screening for
#' ultrahigh-dimensional data. Journal of the American Statistical Association, 106,
#' 1464-1475. <DOI:10.1198/jasa.2011.tm10563>
#' @param n Number of subjects in the dataset to be simulated.  It will also equal to the number of rows
#' in the  dataset to be simulated, because it is assumed that each row represents a different independent and
#' identically distributed subject.
#' @param p Number of predictor variables (covariates) in the simulated dataset.  These covariates
#' will be the features screened by DC-SIS.
#' @param rho The correlation between adjacent covariates in the simulated matrix X.  The within-subject covariance
#' matrix of X is assumed to has the same form as an AR(1) autoregressive covariance matrix,
#' although this is not meant to imply that the X covariates for each subject are in fact a time
#' series.  Instead, it is just used as an example of a parsimonious but nontrivial covariance
#' structure.  If rho is left at the default of zero, the X covariates will be independent
#' and the simulation will run faster.
#' @param HeavyTailedResponse If this is true, Y residuals will be generated to have much heavier
#' tails (more unusually high or low values) then a normal distribution would have.
#' @param heteroskedastic Whether the error variance should be allowed to depend on one of the predictor variables.
#' @param SignalStrength A constant used in the simulation to increase or decrease the signal-to-noise ratio;
#' it was set to 0.5, 1, or 2 for weaker, medium or stronger signal.
#' @return A list with following components:
#'      X Matrix of predictors to be screened. It will have n rows and p columns.
#'      Y Vector of responses.  It will have length n.
#' @importFrom stats rt
#' @export simulateSIRS

#' @examples
#' set.seed(12345678)
#' results <- simulateSIRS()

simulateSIRS <- function(n=200,
                         p=5000,
                         rho=0,
                         HeavyTailedResponse=TRUE,
                         heteroskedastic=TRUE,
                         SignalStrength=1) {
  if ((p<30)|(p>100000)) {stop("Please select a number p of predictors between 30 and 100000.")}
  if ((rho<0)|(rho>=1)) {stop("Please select a rho parameter in the interval [0,1).")}
  if (rho != 0) {
    SigmaMatrix <- matrix(0,p,p)
    SigmaMatrix <- rho^abs(row(SigmaMatrix)-col(SigmaMatrix));
    # This is the within-subject covariance matrix of the X covariates.  It has the same form
    # as an AR(1) autoregressive covariance matrix, although this is not meant to imply that the X
    # covariates for each subject are in fact a time series.  Instead, it is just used as an
    # example of a parsimonious but nontrivial covariance structure.
    ev <- eigen(sigma, symmetric = TRUE)
    # This eigenvalue decomposition takes some time if p is large, and is not necessary
    # in the rho=0 case, so
    if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
      warning("sigma is numerically not positive definite")
    }
    SqrtSigmaMatrix <- ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*% t(ev$vectors)
    # We will generate the X predictors as independent normal, and then transform them
    # to have covariance SigmaMatrix by premultiplying them by SqrtSigmaMatrix.
    X <- matrix(rnorm(n*p),n,p)%*%SqrtSigmaMatrix
  } else {
    # We can save computational time by just treating all the X as independent.
    X <- matrix(rnorm(n*p),n,p)
  }
  if (HeavyTailedResponse) {
    error <- rt(n,1)
  } else {
    error<-rnorm(n,0,1)  # Residuals for Y
  }
  if (heteroskedastic) {
    sigma <- exp(X[,20]+X[,21]+X[,22])
  } else {
    sigma <- sqrt(6.83)
  }
  Y <- SignalStrength*1*X[,1] +
    SignalStrength*0.8*X[,2] +
    SignalStrength*0.6*X[,3] +
    SignalStrength*0.4*X[,4] +
    SignalStrength*0.2*X[,5] +
    sigma*error;
  results <- list()
  results$X <- X
  results$Y <- Y
  return(results)
}
