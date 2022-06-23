#' @title  Perform screening for ultrahigh-dimensional varying coefficient model
#'
#' @description Implements a screening procedure proposed by Liu, Li and
#'  Wu(2014) for varying coefficient models with
#'  ultra-high dimensional predictors.
#'
#'    The function code is adapted from the relevant authors' code. Special thanks are due to
#'    Jingyuan Liu for providing some of the code upon which this function is based.
#'
#' @param X Matrix of predictors to be screened. There should be one row for each
#'          observation.
#' @param Y Vector of responses.  It should have the same length as the number of
#'        rows of X.
#' @param U Covariate, with which coefficient functions vary.
#'
#'
#' @return A list with following components:
#'    CORR_sq A vector of the unconditioned squared correlation with length equal to
#'            the number of columns in the input matrix X. The hgh the unconditioned
#'            squared correlation is, the more desirable it is to retain the corresponding
#'            X covariate in a later predictive model.
#'    rank Vector for the rank of the predictors in terms of the conditional correlation
#'         ( \eqn{\hat{rho}*_j} in the paper). This will have length equal to the number of columns
#'         in the input matrix X, and will consist of a permutation of the integers 1
#'         through that length. A rank of 1 indicates the feature which appears to have
#'         the best marginal predictive performance with largest  \eqn{\hat{rho}*_j}, 2 represents
#'         the second best and so forth.
#'
#' @keywords variable screening
#' @keywords varying-coefficient models
#' @keywords ultra-high dimensional regression
#' @keywords feature selection

#' @references Liu, J., Li, R., & Wu, R. (2014). Feature selection for varying coefficient
#'             models with ultrahigh-dimensional covariates. Journal of the American Statistical
#'             Association, 109: 266-274. <DOI:10.1080/01621459.2013.850086>
#'
#' @export
#' @examples
#' set.seed(12345678)
#' results <- simulateVCM(p=400,
#'                        trueIdx = c(2, 100, 300),
#'                        betaFun = function(U) {
#'                          beta2 <- 2*I(U>0.4)
#'                          beta100 <- 1+U
#'                          beta300 <- (2-3*U)^2
#'                          return(c(beta2,
#'                                   beta100,
#'                                   beta300))
#'                          })
#' screenResults<- screenVCM(X = results$X,
#'                           Y = results$Y,
#'                           U = results$U)
#' rank <- screenResults$rank
#' unlist(rank)
#' trueIdx <- c(2,100,400, 600, 1000)
#' rank[trueIdx]


screenVCM <- function(X, Y, U) {

  if (nrow(X) != length(Y) || nrow(X) != length(U) ) {
    stop("X, Y and U should have same number of rows!")
  }

  Corr<-function(Xj)                                                                              ## function of estimating rho.star.hat using kernel regression
  {                                                                                               ## first choose the common bandwidth h
    resp<-Y*Xj;                                                                                    ## here use the bandwidth for estimating E(Y*xj|U)
    fit.poly<-lm(resp~U+I(U^2)+I(U^3)+I(U^4))                                                      ## use plug-in method (Ruppert, Sheather and Wang, 1995)
    coef<-fit.poly$coeff                                                                           ## To "plug in" the unknown quantities when computing h,
    al2<-coef[3]; al3<-coef[4]; al4<-coef[5]                                                       ## fit a polynomial model Y*xj ~ polynomials of U
    E.second<-mean((2*al2+6*al3*U+12*al4*U^2)^2)                                                   ## plug-in est of second order derivative of E(Y*Xj|U)
    MSE<-(summary(fit.poly)$sigma)^2                                                               ## plug-in est of sigma^2
    n <-length(Y)
    h<-(60*(range(U)[2]-range(U)[1])*MSE/(n*E.second))^0.2                                         ## formula of h with Epanechnikov kernel and uniform U
    K<-function(u)                                                                                 ## weight matrix for a given U
    {
      w<-.75*(1-((U-u)/h)^2)*I(abs(U-u)<=h)                                                         ## Epanechnikov (quadratic) kernel function
      K<-w/(sum(w))                                                                                 ## rescaled kernal function
      return(K)
    }
    wtK<-apply(as.matrix(U),1,K)                                                                   ## kernel matrix
    wtKt<-t(wtK)                                                                                   ## transpose of weight matrix
    EY<-wtKt%*%Y; VarY<-wtKt%*%((Y-EY)^2); EXj<-wtKt%*%Xj; VarXj<-wtKt%*%((Xj-EXj)^2);             ## conditional mean and variance estimations
    YXj<-Y*Xj; EYXj<-wtKt%*%YXj; CovYXj<-EYXj-EY*EXj; CorrYXj<-CovYXj/sqrt(VarY*VarXj);            ## and conditional correlation estimation
    Corr.unc<- mean(CorrYXj^2)                                                                     ## unconditioned squared correlation rho.star.hat
    return(Corr.unc)
  }

  CORR_sq<-apply(X,2,Corr)                                                                        ## compute rho.star.hat for each column of X matrix
  CORR_sq_ord<-CORR_sq[order(CORR_sq, decreasing=TRUE)]
  rank<-match(CORR_sq, CORR_sq_ord)
  results <- list()
  results$rank <- rank

  return(results)

}
