% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulateMVSIS.R
\name{simulateMVSIS}
\alias{simulateMVSIS}
\title{Simulate a dataset for demonstrating the performance of screenIID with the MV-SIS option with categorical outcome variable}
\usage{
simulateMVSIS(R = 2, n = 40, p = 2000, mu = 3, heavyTailedCovariates = FALSE)
}
\arguments{
\item{R}{a positive integer, number of outcome categories for multinomial (categorical) outcome Y.}

\item{n}{Number of subjects in the dataset to be simulated.  It will also equal to the number of rows
in the  dataset to be simulated, because it is assumed that each row represents a different independent and
identically distributed subject.}

\item{p}{Number of predictor variables (covariates) in the simulated dataset.  These covariates
will be the features screened by DC-SIS.}

\item{mu}{Signal strength;  the larger mu is, the easier the active covariates will be to discover.
#  Specifically, mu is added to the rth predictor for r=1,...,R, so that the probability
that Y equals r will be higher if the rth predictor is higher.  It is assumed that p>>r
so that most predictors will be inactive.  In real data there is no reason why, say, the
first two columns in the matrix should be the important ones, but this is convenient in a
simulation and the choice of permutation of the columns involves no loss of generality.}

\item{heavyTailedCovariates}{If TRUE, the covariates will be generated as independent t variates,
plus covariate-specific constants.  If FALSE, they will be generated as independent standard normal variates.}
}
\value{
A list with following components:
     X Matrix of predictors to be screened. It will have n rows and p columns.
     Y Vector of responses.  It will have length n.
}
\description{
Simulates a dataset that can be used to test screenIID for ultrahigh-dimensional discriminant analysis with the MV-SIS option.
The simulation is based on the balanced scenarios in Example 3.1 of Cui, Li & Zhong (2015).
The simulated dataset has p numerical X-predictors and a categorical Y-response.
 Special thanks are due to Wei Zhong for providing some of the code upon which this function is based.
}
\examples{
set.seed(12345678)
results <- simulateMVSIS()
}
\references{
Cui, H., Li, R., & Zhong, W. (2015). Model-free feature screening
for ultrahigh dimensional discriminant analysis. Journal of the American
Statistical Association, 110: 630-641. <DOI:10.1080/01621459.2014.920256>
}
