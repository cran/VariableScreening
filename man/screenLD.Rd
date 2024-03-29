% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/screenLD.R
\name{screenLD}
\alias{screenLD}
\title{Perform high-dimensional screening for semiparametric longitudinal regression}
\usage{
screenLD(
  X,
  Y,
  z,
  id,
  subset = 1:ncol(X),
  time,
  degree = 3,
  df = 4,
  corstr = "stat_M_dep",
  M = NULL
)
}
\arguments{
\item{X}{Matrix of features (for example, SNP's).
There should be one row for each observation.}

\item{Y}{Vector of responses.  It should have the same length as the number of
rows of X.}

\item{z}{Optional matrix of covariates to be included in all models. They may
include demographic covariates such as gender or ethnic background, or
some other theoretically important constructs.  It should have the same
number of rows as the number of rows of X.  We suggest a fairly low
dimensional z.  If the model is intended to include an intercept function
(which is recommended), then z should include a column of 1's representing
the constant term.}

\item{id}{Vector of integers identifying the subject to which each observation belongs.
It should have the same length as the number of rows of X.}

\item{subset}{Vector of integers identifying a subset of the features of X to be
screened, the default is 1:ncol(X), i.e., to screen all columns of X.}

\item{time}{Vector of real numbers identifying observation times. It should have the
same length as the number of rows of X.  We suggest using the convention of
scaling time to the interval [0,1].}

\item{degree}{Degree of the piecewise polynomial for the B-spline basis for the
varying coefficient functions; see the documentation for the bs() function
in the splines library.}

\item{df}{Degrees of freedom of the B-spline basis for the varying coefficient
functions; see the documentation for the bs() function in the splines library.}

\item{corstr}{Working correlation structure for the generalized estimation equations
model used to estimate the coefficient functions; see the documentation for
the gee() function in the gee library. Options provided by the gee() function
include "independence", "fixed", "stat_M_dep", "non_stat_M_dep", "exchangeable",
"AR-M" and "unstructured".}

\item{M}{An integer indexing the M value (complexity) of the dependence structure,
if corstr is M-dependent or AR-M; see the documentation for the gee() function
in the gee library.  This will be ignored if the correlation structure does not
require an M parameter. The default value is set to be 1.}
}
\value{
A list with following components:
   error A vector of length equal to the number of columns in the input matrix X.
       It contains sum squared error values for regression models which include the
       time-varying effects of the z covariates (if any) as well as each X covariate by
       itself.  The lower this error is, the more desirable it is to retain the
       corresponding X covariate in a later predictive model.
   rank The rank of the error measures.  This will have length equal to the number
       of columns in the input matrix X, and will consist of a permutation of the integers
       1 through that length.  A rank of 1 indicates the feature which appears to have
       the best predictive performance, 2 represents the second best and so on.
}
\description{
Implements a screening procedure proposed by Chu, Li, and Reimherr (2016)
<DOI:10.1214/16-AOAS912> for varying coefficient longitudinal models with ultra-high
dimensional predictors. The effect of each predictor is allowed to vary over time,
approximated by a low-dimensional B-spline.  Within-subject correlation is handled
using a generalized estimation equation approach with structure specified by the user.
Variance is allowed to change over time, also approximated by a B-spline.
}
\examples{
set.seed(12345678)
results <- simulateLD(p=500)
subset1 <- seq(1,5,2)
subset2 <- seq(100,200,2)
subset3 <- seq(202,400,2)
subset4 <- seq(401,499,2)
set <-c(subset1,subset2,subset3,subset4)
Jmin <- min(table(results$id)) - 1
screenResults <- screenLD(X = results$X,
                            Y = results$Y,
                            z = results$z,
                            id = results$id,
                            subset = set,
                            time = results$time,
                            degree = 3,
                            df = 4,
                            corstr = "stat_M_dep",
                            M = Jmin
                            )


rank <- screenResults$rank
unlist(rank)
trueIdx <- c(5,100,200,400)
rank[which(set \%in\% trueIdx)]
}
\keyword{feature}
\keyword{high-dimensional}
\keyword{regression}
\keyword{screening}
\keyword{selection}
\keyword{variable}
