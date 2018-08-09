#' @title Simulate a dataset for demonstrating the performance of screenIID with the MV-SIS method with numeric outcome Y
#'
#' @description Simulates a dataset that can be used to demonstrate variable screening
#' for ultrahigh-dimensional regression with categorical predictors and numerical outcome
#' variable using the MV-SIS-NY option in screenIID.
#' The simulated dataset has p numerical predictors X and a categorical response Y.
#' The X covariates are generated as binary with success probability 0.5 each.  The
#' response Y is generated as
#'        Y = 5*X1 + 5*X2  + 5*X12 + 5*X22 + e
#' if heteroskedastic=FALSE, where e is a standard normal error term and 1{} is a zero-one
#' indicator function for the truth of the statement contained.
#'  Special thanks are due to Wei Zhong for providing some of the code upon which this function is based.
#' @references Cui, H., Li, R., & Zhong, W. (2015). Model-free feature screening
#' for ultrahigh dimensional discriminant analysis. Journal of the American
#' Statistical Association, 110: 630-641. <DOI:10.1080/01621459.2014.920256>

#' @param n Number of subjects in the dataset to be simulated.  It will also equal to the number of rows
#' in the  dataset to be simulated, because it is assumed that each row represents a different independent and
#' identically distributed subject.
#' @param p Number of predictor variables (covariates) in the simulated dataset.  These covariates
#' will be the features screened by DC-SIS.
#' @return A list with following components:
#'      \item{X:}{Matrix of predictors to be screened. It will have n rows and p columns.}
#'      \item{Y:}{Vector of responses.  It will have length n.}
#' @export simulateMVSISNY

#' @examples
#' set.seed(12345678)
#' results <- simulateMVSISNY()

simulateMVSISNY <- function(n=500,
                          p=1000) {
  if ((p<30)|(p>100000)) {stop("Please select a number p of predictors between 30 and 100000.")}
  X <- matrix(rbinom(n*p,1,.5),n,p)
  error<-rnorm(n,0,1)  # Residuals for Y
  Y <- 5*X[,1] + 5*X[,2] + 5*X[,12] + 5*X[,22] + error
  results <- list()
  results$X <- X
  results$Y <- Y
  return(results)
}
