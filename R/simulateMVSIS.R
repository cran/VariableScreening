#' @title Simulate a dataset for demonstrating the performance of screenIID with the MV-SIS option with categorical outcome variable
#'
#' @description Simulates a dataset that can be used to test screenIID for ultrahigh-dimensional discriminant analysis with the MV-SIS option.
#' The simulation is based on the balanced scenarios in Example 3.1 of Cui, Li & Zhong (2015).
#' The simulated dataset has p numerical X-predictors and a categorical Y-response.
#'  Special thanks are due to Wei Zhong for providing some of the code upon which this function is based.
#'
#' @references Cui, H., Li, R., & Zhong, W. (2015). Model-free feature screening
#' for ultrahigh dimensional discriminant analysis. Journal of the American
#' Statistical Association, 110: 630-641. <DOI:10.1080/01621459.2014.920256>

#' @param R a positive integer, number of outcome categories for multinomial (categorical) outcome Y.
#' @param n Number of subjects in the dataset to be simulated.  It will also equal to the number of rows
#' in the  dataset to be simulated, because it is assumed that each row represents a different independent and
#' identically distributed subject.
#' @param p Number of predictor variables (covariates) in the simulated dataset.  These covariates
#' will be the features screened by DC-SIS.
#' @param mu Signal strength;  the larger mu is, the easier the active covariates will be to discover.
#'  #  Specifically, mu is added to the rth predictor for r=1,...,R, so that the probability
#'  that Y equals r will be higher if the rth predictor is higher.  It is assumed that p>>r
#'  so that most predictors will be inactive.  In real data there is no reason why, say, the
#'  first two columns in the matrix should be the important ones, but this is convenient in a
#'  simulation and the choice of permutation of the columns involves no loss of generality.
#' @param heavyTailedCovariates  If TRUE, the covariates will be generated as independent t variates,
#'  plus covariate-specific constants.  If FALSE, they will be generated as independent standard normal variates.
#' @return A list with following components:
#'      X Matrix of predictors to be screened. It will have n rows and p columns.
#'      Y Vector of responses.  It will have length n.
#' @importFrom stats rt
#' @export simulateMVSIS
#' @examples
#' set.seed(12345678)
#' results <- simulateMVSIS()


simulateMVSIS <- function(R = 2,  # number of outcome categories for multinomial (categorical) outcome Y
                          n = 40,  # number of subjects to be generated
                          p = 2000,
                          mu = 3 , # mu = signal strength
                          heavyTailedCovariates = FALSE) {
  # Simulates multinomial Y for demonstrating MVSIS discriminant analysis.
  #  The p columns of X are each generated as independent;  most have
  #  means of zero, but the active covariates have different means depending on the level of Y.
  if ((p<30)|(p>100000)) {stop("Please select a number p of predictors between 30 and 100000.")}

  Y<-sample(1:R,n,replace=T,prob=rep(1/R,R))
  X<-matrix(0,n,p)
  for (r in 1:R) {
    ind<-(Y==r)
    A<-matrix(0,sum(ind),p)
    for(i in 1:sum(ind)){
      # That is, for each subject having Y=r...
      if (heavyTailedCovariates) {
        A[i,] <- c(seq(0,0,length=r-1),mu,seq(0, 0,length=p-r))+rt(p,2)
        # If heavyTailedCovariates = TRUE, then covariates will be generated as a mean
        # plus an independent error following a t distribution with 2 degrees of freedom.
      } else {
        A[i,] <- c(seq(0,0,length=r-1),mu,seq(0, 0,length=p-r))+rnorm(p)
        # If heavyTailedCovariates = TRUE, then covariates will be generated as a mean
        # plus a standard normal distribution, i.e., they will be independently normally
        # distributed with variance 1.
      }
    }
    X[ind,]<-A
  }
  #  The p columns of X are each generated as independent normals with variance one;  most have
  #  means of zero, but the active covariates have different means depending on the level of Y.
  #  Specifically, mu has been added to the rth predictor for r=1,...,R, so that the probability
  #  that Y equals r will be higher if the rth predictor is higher.  It is assumed that p>>r
  #  so that most predictors will be inactive.  In real data there is no reason why, say, the
  #  first two columns in the matrix should be the important ones, but this is convenient in a
  #  simulation and the choice of permutation of the columns involves no loss of generality.
  results <- list()
  results$X <- X
  results$Y <- Y

  return(results)

}
