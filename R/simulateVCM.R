#' @title Simulate a dataset for testing the performance of screenVCM
#'
#' @description Simulates a dataset that can be used to test the screenVCM
#' function, and to test the performance of the proposed method under different
#' scenarios. The simulated dataset has a single U-covariate and p X-predictors,
#' only a few of which have nonzero effect.
#'
#    The function code is adapted from the relevant authors' code. Special thanks are due to
#'   Jingyuan Liu for providing some of the code upon which this function is based.
#'
#' @param n Number of subjects in the simulated dataset
#' @param rho The correlation matrix of columns of X.
#' @param p The total number of features to be screened from
#' @param trueIdx The indexes for the active features in the simulated X matrix.
#'                This should be a vector, and the values should be a subset of 1:p.
#' @param betaFun A list of functions of U, one function for each entry in trueIdx,
#'                giving the varying effects of each active predictor in the simulated
#'                X matrix.
#' @importFrom MASS mvrnorm
#' @importFrom stats rbinom rnorm runif pnorm
#'
#' @return A list with following components:
#'      X Matrix of predictors to be screened. It will have n rows and p columns.
#'      Y Vector of responses.  It will have length of n.
#'      U A vector representing a covariate with which the coefficient functions vary.
#' @export simulateVCM
#' @examples
#' set.seed(12345678)
#' results <- simulateVCM(p=1000)


simulateVCM <- function(n = 200,
                       rho = 0.4,
                       p = 1000,
                       trueIdx = c(2, 100, 400, 600, 1000),
                       betaFun = NULL) {
  if ((p<30)|(p>100000)) {stop("Please select a number p of predictors between 30 and 100000.")}
  if ((rho<0)|(rho>=1)) {stop("Please select a rho parameter in the interval [0,1).")}
  if( p < 1000 && is.null(betaFun)){stop("If p is set to be less than 1000, then betaFun has to be user specified.")}
  if (is.null(betaFun)) {
    betaFun <- function(U) {
      beta2 <- 2*I(U>0.4)

      beta100 <- 1+U

      beta400 <- (2-3*U)^2

      beta600 <- 2*sin(2*pi*U)

      beta1000 <- exp(U/(U+1))

      return(c(beta2, beta100, beta400, beta600, beta1000))

    }
  }

  d0<-floor((n^.8)/(log(n^.8)))                                                                   ## screened submodel size
  sigma<-1                                                                                        ## variance of xj's
  covm<-diag(p+1)                                                                                 ## construct covariance matrix of X
  covm<-sigma^2*rho^abs(row(covm)-col(covm))
  ev<-eigen(covm, symmetric = TRUE)
  if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
    warning("sigma is numerically not positive definite")
  }
  covmroot<- ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*% t(ev$vectors)            ## squre root of covariance matrix of X


    predictor<-matrix(rnorm(n*(p+1)),n,(p+1))%*%covmroot;                                           ## jointly draw u_star and xj's from Multivariate Normal
    U<-pnorm(predictor[,1]);                                                                        ## generate U from u_star, then U~Unif(0,1)
    coefMat <- t(sapply(U, betaFun))                                                               ## generate nonzero beta(U)'s as functions of U:
    X<-predictor[,2:(p+1)]                                                                          ## X matrix (correlated with U)
    eps<-rnorm(n,0,1)                                                                                    ## generate noise epsilon
    xMat <- X[, trueIdx]
    Y<-rowSums(coefMat*xMat)+eps          ## construct response Y using Xj's, U and epsilon

  results <- list()
  results$X <- X
  results$Y <- Y
  results$U <- U

  return(results)

}
