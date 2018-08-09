#' @title  Feature Selection for Ultrahigh-Dimensional Datasets with Independent Subjects,
#'
#' @description Implements one of three screening procedures: Sure Independent Ranking and
#'    Screening (SIRS), Distance Correlation Sure Independence Screening (DC-SIS), or
#'    MV Sure Independence Screening (MV-SIS).  In general they are extensions of the sure
#'    independence screening concept proposed by Fan and Lv (2008), but without a parametric
#'    assumption (e.g., linear or logistic) on the relationship between the predictor variables X
#'    and outcome Y.
#'
#'    Screening methods each rank the predictors based on some measure
#'    of their estimated strength of relationship with Y.  The assumption is that only a few
#'    among the top-ranked variables are likely to be truly significant predictors.
#'
#'    The original version of SIS involved ranking the predictors by their correlation
#'    with Y, implying a linear relationship.  The SIRS method is an extension proposed
#'    by Zhu, Li, Li, & Zhu (2011), which involved ranking the predictors by their
#'    correlation with the rank-ordered Y instead, thereby not assuming a linear
#'    correlation, and potentially outperforming SIS.
#'
#'    DC-SIS was then proposed by Li, Zhong and Zhu (2012) and its relationship measure is
#'    the distance correlation (DC) between a covariate and the outcome, a nonparametric
#'    generalization of the correlation coefficient (Szekely, Rizzo, & Bakirov, 2007).
#'    The function uses the dcor function from the R package energy in order to
#'    calculate this correlation.  Simulations showed that DC-SIS could sometimes provide a further
#'    advantage over SIRS.
#'
#'    The above measures were primarily intended for a numerical Y.  Cui, Li, and Zhong (2015)
#'    proposed MV-SIS, which was developed for categorical Y (including binary Y) as in
#'    discriminant analysis, and which is also robust to heavy-tailed predictor distributions.
#'    The measure used by MV-SIS for the association strength between a particular Xk and Y is a
#'    mean conditional variance measure called MV for short,  namely the expectation in X of
#'    the variance in Y of the conditional cumulative distribution function F(x|Y)=P(X<=x|Y);
#'    note that like the correlation or distance correlation, this is zero if X and Y are
#'    independent because F(x) does not depend on Y in that case.  Cui, Li, and Zhong (2015)
#'    also point out that the MV-SIS can alternatively be used with categorical X variables
#'    and numerical Y, instead of numerical X and categorical Y.  This function supports that
#'    option as "MV-SIS-NY."
#'
#'    Whichever option is chosen, the function returns the ranking of the predictors according
#'    to the appropriate association measure.
#'
#'    The function code is adapted from the relevant authors' code.  Special thanks are due to
#'    Wei Zhong for providing some of the code upon which this function is based.
#'
#' @param X Matrix of predictors to be screened. There should be one row for each
#'          observation.
#' @param Y Vector of responses.  It should have the same length as the number of
#'        rows of X.  The responses should be numerical if SIRS or DC-SIS is used.
#'        The responses should be integers representing response categories if MV-SIS is used.
#'        Binary responses can be used for any method.
#'
#' @param method Screening method.  The options are "SIRS", "DC-SIS", "MV-SIS" and "MV-SIS-NY", as
#'      described above.
#'
#' @importFrom energy dcor
#' @importFrom stats sd
#'
#' @return A list with following components:
#'
#'    \item{measurement:}{A vector of length equal to the number of columns in the input matrix X.
#'                It contains estimated strength of relationship with Y.}
#'
#'    \item{rank:}{The rank of the error measures.  This will have length equal to the number
#'        of columns in the input matrix X, and will consist of a permutation of the integers
#'        1 through that length.  A rank of 1 indicates the feature which appears to have
#'        the best predictive performance, 2 represents the second best and so on.}
#'
#' @keywords Variable screening, variable selection, high-dimensional regression, high-dimensional discriminant analysis
#'
#' @references Cui, H., Li, R., & Zhong, W. (2015). Model-free feature screening for
#'             ultrahigh dimensional discriminant analysis. Journal of the American
#'             Statistical Association, 110: 630-641. <DOI:10.1080/01621459.2014.920256>
#'
#'             Fan, J., & Lv, J. (2008). Sure independence screening for ultrahigh
#'             dimensional feature space. Journal of the Royal Statistical Society, B,
#'             70: 849-911. <DOI:10.1111/j.1467-9868.2008.00674.x>
#'
#'             Li, R., zhong, W., & Zhu, L. (2012). Feature screening via distance
#'             correlation learning. Journal of the American Statistical Association,
#'             107: 1129-1139. <DOI:10.1080/01621459.2012.695654>
#'
#'             Szekely, G. J., Rizzo, M. L., & Bakirov, N. K. (2007). Measuring and
#'             Testing Dependence by Correlation of Distances. Annals of Statistics,
#'             35, 2769-2794. <DOI: 10.1214/009053607000000505>
#'
#'             Zhu, L.-P., Li, L., Li, R., & Zhu, L.-X. (2011) Model-free feature
#'             screening for ultrahigh-dimensional data. Journal of the American
#'             Statistical Association, 106: 1464-1475. <DOI:10.1198/jasa.2011.tm10563>
#'
#'
#' @export screenIID
#' @examples
#' set.seed(12345678)
#' data1 <- simulateDCSIS(n=100,p=1000)
#' answers<- screenIID(X = data1$X, Y = data1$Y, method="DC-SIS")
#' print(which(answers$rank<=10)) # Print the columns of X corresponding to the ten best-ranked
#'                                # predictors. Note that in the simulation, the true predictors
#'                                # are columns 1, 2, 12, and 22, and three of these are included,
#'                                # indicating that the function worked fairly well.




screenIID <- function(X, Y, method = "DC-SIS") {
  if (nrow(X) != length(Y) ) {
    stop("X and Y should have same number of rows!")
  }

  # Make the specification of the method a little more user-friendly:
  method <- toupper(method)
  if (method=="DCSIS") method <- "DC-SIS"
  if (method=="MVSIS") method <- "MV-SIS"
  if (method=="MVSISNY") method <- "MV-SIS-NY"

  # The follwing three functions are defined for use in the MV-SIS and MV-SIS-NY methods.

  Fk<-function(X0,x) {
    # Returns the estimated value of F(X0) = P(X <= X0), specifically the proportion
    # of values in x that are less than a constant X0, where x is a numerical vector
    # representing the observed values of the kth covariate x.  This is the estimated
    # unconditional cumulative distribution function of X.
    Fk=c()
    for (i in 1:length(x))
    { Fk[i]=sum(X0<=x[i])/length(X0) }
    return(Fk)
  }

  Fkr<-function(X0,Y,yr,x) {
    # Returns the estimated value of Fr(X0) = P(X <= X0| Y=yr),
    # for variables (actually vectors) x and Y and for constants X0 and yr.
    # This is the estimated cumulative distribution function of X conditional
    # on the rth value of Y.
    Fkr=c()
    ind_yr=(Y==yr)
    for (i in 1:length(x))
    { Fkr[i]=sum((X0<=x[i])*ind_yr)/sum(ind_yr) }
    return(Fkr)
  }

  MV<-function(X,Y) {
    # Calculates the MV nonparametric association measure between numerical variable X
    # and categorical variable Y, defined as expression (2.2) in Cui, H., Li, R., & Zhong, W. (2015).
    Fk0 <- Fk(X,X)
    # It is not a mistake that X is used twice here.  For each value of this covariate, we want
    # to know what quantile it represents.  The MV  is the Expectation in X of the variance in Y
    # of the conditional cumulative distribution function F(x|Y)=P(X<=x|Y).
    # Note that like the correlation or distance correlation, this is zero if X and Y are
    # independent because F(x) does not depend on Y in that case.
    Yr <- unique(Y)
    if (length(Yr)>15) {stop("A supposedly categorical variable was provided with more than 15 levels.  Did you mean to treat it as numeric?");}
    MVr <- c()
    for (r in 1:length(Yr)) {
      MVr[r] <- (sum(Y==Yr[r])/length(Y))*mean((Fkr(X,Y,Yr[r],X)-Fk0)^2)
    }
    MV <- sum(MVr)
    return(MV)
  }

  results <- list()

  if( method == "DC-SIS" ){
    dcory<-function(Xk){dcor(Xk,Y)}  # The distance correlation measure is from the energy
    # package by Rizzo and Szekely.
    v<-abs(apply(X,2,dcory))
    Dword <- order(v,decreasing = T)   # rank order of estimated association strengths;
    rank<-match(v, v[Dword])

    results$measurement = v
    results$rank <- rank

  } else if(method == "MV-SIS"){

    MVy<-function(Xk){MV(Xk,Y)}
    v<-abs(apply(X,2,MVy))
    MWord <- order(v,decreasing = T)      # rank order of estimated association strengths;
    rank<-match(v, v[MWord])

    results$measurement = v
    results$rank <- rank
  } else if(method == "SIRS"){

    n <- length(Y)
    p <- ncol (X)
    XX <- X
    for (j in 1:p) {
      XX[,j]<-(X[,j]-mean(X[,j]))/sd(X[,j])
    }

    yord <- order(Y)
    SIRSk<-function(Xk) {
      z<-c()
      Xk.ord=Xk[yord]
      for (j in 1:n) {
        z[j]=(sum(Xk.ord[1:j])/n)^2
      }
      return(mean(z))
    }
    v <- abs(apply(XX,2,SIRSk))
    Zword<- order(v,decreasing = T)       # rank order of estimated association strengths;
    rank<-match(v, v[Zword])

    results$measurement = v
    results$rank <- rank
  } else if(method == "MV-SIS-NY"){
    # This is exactly the same as MV-SIS, but with Xk and Y reversed in the call to the MV function.

    MVy<-function(Xk){MV(Y,Xk)}
    v<-abs(apply(X,2,MVy))
    MWord <- order(v,decreasing = T)      # rank order of estimated association strengths;
    rank<-match(v, v[MWord])

    results$measurement = v
    results$rank <- rank
  } else {
    stop("The method requested was not recognized.")
  }
  return(results)


}
