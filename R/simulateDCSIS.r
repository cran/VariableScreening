#' @title Simulate a dataset for demonstrating the performance of screenIID with the DC-SIS method
#'
#' @description Simulates a dataset that can be used to demonstrate variable screening
#' for ultrahigh-dimensional regression with the DC-SIS option in screenIID.
#' The simulated dataset has p numerical predictors X and a categorical Y-response.
#' The data-generating scenario is a simplified version of Example 3.1a (homoskedastic) or 3.1d
#' (heteroskedastic) of  Li, Zhong & Zhu (2012).
#' Specifically, the X covariates are normally distributed with mean zero and variance one, and
#' may be correlated if the argument rho is set to a nonzero value.  The response Y is generated as
#' either
#'        Y = 6*X1 + 1.5*X2  + 9*1{X12 < 0} + exp(2*X22)*e
#' if heteroskedastic=TRUE, or
#'        Y = 6*X1 + 1.5*X2  + 9*1{X12 < 0} + 6*X22 + e
#' if heteroskedastic=FALSE, where e is a standard normal error term and 1{} is a zero-one
#' indicator function for the truth of the statement contained.
#'  Special thanks are due to Wei Zhong for providing some of the code upon which this function is based.
#'
#' @references Li, R., Zhong, W., & Zhu, L. (2012) Feature screening via distance correlation learning. Journal of
#' the American Statistical Association, 107: 1129-1139. <DOI:10.1080/01621459.2012.695654>
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
#' @param heteroskedastic Whether the error variance should be allowed to depend on one of the predictor variables.
#' @param sigma The error standard deviation of the response
#' @return A list with following components:
#'      \item{X:}{Matrix of predictors to be screened. It will have n rows and p columns.}
#'      \item{Y:}{Vector of responses.It will have length n.}
#' @export simulateDCSIS

#' @examples
#' set.seed(12345678)
#' results <- simulateDCSIS()

simulateDCSIS <- function(n=200,
                          p=5000,
                          rho=0,
                          heteroskedastic=TRUE,
                          sigma = 1) {
  if ((p<30)|(p>100000)) {stop("Please select a number p of predictors between 30 and 100000.")}
  if ((rho<0)|(rho>=1)) {stop("Please select a rho parameter in the interval [0,1).")}
  if (rho != 0) {
    X <- matrix(NA,n,p)
    X[,1] <- rnorm(n)
    for (i in 2:p)  {
      X[,i]=rho*X[,i-1]+sqrt(1-rho^2)*rnorm(n)
    }
  } else {
    # We can save computational time by just treating all the X as independent.
    X <- matrix(rnorm(n*p),n,p)
  }
  # Models loosely based on Example 1a (homoskedastic) and Example1d (heteroskedastic)
  # in section 3 of Li, Zhong & Zhu (2012)
  error<-rnorm(n,0,sigma)  # Residuals for Y
  if (heteroskedastic) {
    Y <- 6*X[,1] + 1.5*X[,2] + 9*(X[,12]<0) + exp(2*X[,22])*error;
  } else {
    Y <- 6*X[,1] + 1.5*X[,2] + 9*(X[,12]<0) + 6*X[,22] + error;
  }
  results <- list()
  results$X <- X
  results$Y <- Y
  return(results)
}
