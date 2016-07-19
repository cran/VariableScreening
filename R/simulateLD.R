#' Simulate a dataset for testing the performance of screenlong
#'
#' Simulates a dataset that can be used to test the screenlong function,
#' and to test the performance of the proposed method under different scenarios.
#' The simulated dataset has two z-covariates and p x-covariates, only a
#' few of which have nonzero effect.  There are n subjects in the simulated
#' dataset, each having J observations, which are not necessarily evenly timed,
#' we randomly draw a subset to create an unbalanced dataset. The within-subject
#' correlation is assumed to be AR-1.
#'
#' @param n Number of subjects in the simulated dataset
#' @param J Number of observations per subject
#' @param rho The correlation parameter for the AR-1 correlation structure.
#' @param p The total number of features to be screened from
#' @param trueIdx The indexes for the active features in the simulated x matrix.
#' This should be a vector, and the values should be a subset of 1:p.
#' @param beta0Fun The time-varying intercept for the data-generating model, as a function of
#' time.  If left as null, it will default to \code{f(t) 2 * t^2 - 1}. Time is assumed to be scaled
#' to the interval [0,1].
#' @param betaFun The time-varying coefficients for z in the data-generating model, as a
#' function of time.  If left as null, it will be specified as two functions. The first is
#' \code{f(t) exp(t + 1)/2}. The second is \code{f(t) t^2 + 0.5}. Time is assumed
#' to be scaled to the interval [0,1].
#' @param gammaFun A list of functions of time, one function for each entry in trueIdx,
#' giving the time-varying effects of each active feature in the simulated x matrix.
#' If left as null, it will be specified as four functions.  The first is a step function
#' \code{f(t)=(t > 0.4)}.  The second is \code{f(t)=- cos(2 * pi * t)}.  the third is \code{f(t)=(2 - 3 * t)^2/2 - 1}.
#' The fourth is \code{f(t)=sin(2 * pi * t)}.
#' @param varFun A function of time telling the marginal variance of the error function at a
#' given time.  If left as null, it will be specified as \code{function(t) 0.5 + 3 * t^3}.
#' @importFrom MASS mvrnorm
#' @importFrom stats rbinom rnorm runif
#'
#' @return A list with following components:
#'      x Matrix of features to be screened.  It will have n*J rows and p columns.
#'      y Vector of responses.  It will have length of n*J.
#'      z A matrix representing covariates to be included in each of the screening models.
#'        The first column will be all ones, representing the intercept.  The second will
#'        consist of random ones and zeros, representing simulated genders.
#'      id Vector of integers identifying the subject to which each observation belongs.
#'      time Vector of real numbers identifying observation times. It should have the same
#'           length as the number of rows of x.
#' @export
#' @examples
#' set.seed(12345678)
#' results <- simulateLD(p=1000)

simulateLD <- function(n = 100,
                       J = 10,
                       rho = 0.6,
                       p = 500,
                       trueIdx = c(5, 100, 200, 400),
                       beta0Fun = NULL,
                       betaFun = NULL,
                       gammaFun = NULL,
                       varFun = NULL) {
  N <- n * J
  # total number of observations
  id <- rep(1:n, each = J)


  if (is.null(beta0Fun))
    beta0Fun <- function(t)
      2 * t ^ 2 - 1

  if (is.null(betaFun)) {
    betaFun <- function(t) {
      beta1 <- exp(t + 1) / 2
      beta2 <- t ^ 2 + 0.5
      return(c(beta1, beta2))
    }
  }


  if (is.null(gammaFun)) {
    gammaFun <- function(t) {
      gamma1 <- (t > 0.4)

      gamma2 <- -cos(2 * pi * t)

      gamma3 <- (2 - 3 * t) ^ 2 / 2 - 1

      gamma4 <- sin(2 * pi * t)

      return(c(gamma1, gamma2, gamma3, gamma4))

    }
  }

  trueIdx <- c(5, 100, 200, 400)
  if (is.null(varFun))
    varFun <- function(t)
      0.5 + 3 * t ^ 3


  ## Correlation matrix for error term in simulation
  corMat <- diag(1, J)

  for (i in 1:J) {
    for (j in 1:J)
      corMat[i, j] <- rho ^ abs(i - j)

  }

  ## Covariance matrix for error term for each subject  in simulation
  arCovFun <- function(tVec) {
    varVec <- sapply(tVec, varFun)

    covMat <- diag(varVec)

    for (i in 1:J) {
      for (j in 1:J)
        covMat[i, j] <-
          varVec[i] ^ 0.5 * varVec[j] ^ 0.5 * corMat[i, j]

    }
    return(covMat)

  }

  ## Generate x-predictors
  zVec1 <- rnorm(N)
  zVec2 <- rbinom(N, size = 1, prob = 0.5)
  zMat <- matrix(c(zVec1, zVec2), ncol = 2)   # z-predictor

  xMat <- matrix(rnorm(p * N), ncol = p)
  # x-predictor
  tVec <- sort(runif(N))
  # time points for all observations

  ## Generate correlated error with time-varying variance

  errVec <- c()

  for (i in 1:n) {
    idx <- ((i - 1) * J + 1):(i * J)

    covMat_i <- arCovFun(tVec[idx])

    errVec <- c(errVec, mvrnorm(1, rep(0, J), covMat_i))

  }

  coefMat <-
    cbind(t(sapply(tVec, betaFun)), t(sapply(tVec, gammaFun)))


  predMat <- cbind(zMat, xMat[, trueIdx])


  ## Generate response
  yVec <-
    sapply(tVec, beta0Fun) + apply(cbind(coefMat, predMat), 1, function(x)
      sum(x[1:ncol(coefMat)] * x[(ncol(coefMat) + 1):(ncol(coefMat) + ncol(predMat))])) + errVec


  ########################################
  ## Take a subset to create unbalanced data set ##
  ########################################

  obs_index <- sort(sample(1:N, size = 0.8 * N))

  yVec <- yVec[obs_index]
  zMat <- zMat[obs_index, ]
  xMat <- xMat[obs_index,]
  id <- id[obs_index]
  tVec <- tVec[obs_index]

  results <- list()

  results$x <- xMat

  results$y <- yVec

  results$z <- zMat

  results$id <- id

  results$time <- tVec


  return(results)

}
