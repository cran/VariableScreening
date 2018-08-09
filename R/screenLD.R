#' @title  Perform high-dimensional screening for semiparametric longitudinal regression
#'
#' @description Implements a screening procedure proposed by Chu, Li, and Reimherr (2016)
#' <DOI:10.1214/16-AOAS912> for varying coefficient longitudinal models with ultra-high
#' dimensional predictors. The effect of each predictor is allowed to vary over time,
#' approximated by a low-dimensional B-spline.  Within-subject correlation is handled
#' using a generalized estimation equation approach with structure specified by the user.
#' Variance is allowed to change over time, also approximated by a B-spline.
#'
#'    The function code is adapted from the relevant authors' code.  Special thanks are due to
#'    Wanghuan Chu for providing some of the code upon which this function is based.
#'
#' @param X Matrix of features (for example, SNP's).
#'        There should be one row for each observation.
#' @param Y Vector of responses.  It should have the same length as the number of
#'        rows of X.
#' @param z Optional matrix of covariates to be included in all models. They may
#'        include demographic covariates such as gender or ethnic background, or
#'        some other theoretically important constructs.  It should have the same
#'        number of rows as the number of rows of X.  We suggest a fairly low
#'        dimensional z.  If the model is intended to include an intercept function
#'        (which is recommended), then z should include a column of 1's representing
#'        the constant term.
#' @param id Vector of integers identifying the subject to which each observation belongs.
#'        It should have the same length as the number of rows of X.
#' @param time Vector of real numbers identifying observation times. It should have the
#'        same length as the number of rows of X.  We suggest using the convention of
#'        scaling time to the interval [0,1].
#' @param degree Degree of the piecewise polynomial for the B-spline basis for the
#'        varying coefficient functions; see the documentation for the bs() function
#'        in the splines library.
#' @param df Degrees of freedom of the B-spline basis for the varying coefficient
#'        functions; see the documentation for the bs() function in the splines library.
#' @param corstr Working correlation structure for the generalized estimation equations
#'        model used to estimate the coefficient functions; see the documentation for
#'        the gee() function in the gee library. Options provided by the gee() function
#'        include "independence", "fixed", "stat_M_dep", "non_stat_M_dep", "exchangeable",
#'        "AR-M" and "unstructured".
#' @param M An integer indexing the M value (complexity) of the dependence structure,
#'        if corstr is M-dependent or AR-M; see the documentation for the gee() function
#'        in the gee library.  This will be ignored if the correlation structure does not
#'        require an M parameter. The default value is set to be 1.
#' @importFrom gee gee
#' @importFrom splines bs
#' @importFrom expm sqrtm
#' @importFrom stats Gamma fitted gaussian glm lm
#'
#' @return A list with following components:
#'    \item{error:}{A vector of length equal to the number of columns in the input matrix X.
#'        It contains sum squared error values for regression models which include the
#'        time-varying effects of the z covariates (if any) as well as each X covariate by
#'        itself.  The lower this error is, the more desirable it is to retain the
#'        corresponding X covariate in a later predictive model.}
#'    \item{rank:}{The rank of the error measures.  This will have length equal to the number
#'        of columns in the input matrix X, and will consist of a permutation of the integers
#'        1 through that length.  A rank of 1 indicates the feature which appears to have
#'        the best predictive performance, 2 represents the second best and so on.}
#' @export screenLD
#' @examples
#' set.seed(12345678)
#' data1 <- simulateLD(p=250,trueIdx = c(5, 100, 200))
#' Jmin <- min(table(data1$id)) - 1
#' screenResults <- screenLD(X = data1$X,
#'                             Y = data1$Y,
#'                             z = data1$z,
#'                             id = data1$id,
#'                             time = data1$time,
#'                             degree = 3,
#'                             df = 4,
#'                             corstr = "stat_M_dep",
#'                             M = Jmin
#'                             )
#' rank <- screenResults$rank
#' print(which(rank<=10))   # Show the column indexes of the best 10 predictors;
#' trueIdx <- c(5,100,200)
#' print(rank[trueIdx])  # Show the ranks given to the four truly active predictors;
#' # It is found that they all are given high ranks, indicating the function works well.



screenLD <- function(X,
                     Y,
                     z=NULL,
                     id,
                     time,
                     degree = 3,
                     df = 4,
                     corstr = "stat_M_dep",
                     M = NULL) {
  sort_idx <- order(id)
  id <- id[sort_idx]
  X <- X[sort_idx, ]
  Y <- Y[sort_idx]
  time <- time[sort_idx]
  if (!is.null(z)) {
    z <- z[sort_idx, ,drop=FALSE]
    if (nrow(X) != length(Y) || nrow(X) != nrow(z) || nrow(X) != length(id)) {
      stop("X, Y, z and id should have same number of rows!")
    }
  } else {
    if (nrow(X) != length(Y) || nrow(X) != length(id)) {
      stop("X, Y, and id should have same number of rows!")
    }
  }
  if (corstr == "stat_M_dep" || corstr == "AR-M") {
    if (!is.null(M)) {
      if (M >= min(table(id))) {
        stop(paste(
          "M is too big, should be smaller than the minimum cluster size",
          min(table(id))
        ))
      }
    } else
      M <- min(table(id)) - 1
  }

  N <- length(Y)
  n <- length(unique(id))
  p <- ncol(X)

  ## Function used to generate basis
  Beta1BsGen <- function(u)
    u[-1] * u[1]


  bsInt <- cbind(rep(1, N), bs(time, degree = degree, df = df))
  bsBase <- data.frame(bsInt)

  if (!is.null(z)) {
    for (i in 1:ncol(z)) {
      bsZ <- t(apply(cbind(as.vector(z[, i]), bsInt), 1, Beta1BsGen))
      index1 = ncol(bsBase) + 1
      index2 = index1 + ncol(bsZ) - 1
      bsBase[, index1:index2] <- bsZ
    }
  }


  ## Estimate time-varying variance
  baseModel <- lm(Y ~ 0 + ., data = bsBase)

  res <- baseModel$residuals

  varEst <- fitted(glm(res^2 ~ 0 + ., data = bsBase, family = Gamma(link = "log")))

  ## adjust the variance of Y
  y_va <- Y / (varEst ^ 0.5)

  bsBase_va <- apply(bsBase, 2, function(u)
    u / (varEst ^ 0.5))


  ## Estimate correlation structure

  corModel <-
    gee(
      y_va ~ 0 + .,
      data = as.data.frame(bsBase_va),
      id = id,
      family = gaussian("identity"),
      corstr = corstr,
      Mv = M
    )

  workCorMat <- corModel$working.correlation


  workCovMatList <- list()

  workCovMatList_inv2 <- list()

  y_rw <- y_va

  bsBase_rw <- bsBase

  Jcs <- cumsum(table(id))

  for (i in 1:n) {
    if (i == 1)
      idx <- 1:(Jcs[i])
    else
      idx <- (Jcs[i - 1] + 1):Jcs[i]
    J_i <- length(idx)
    varEst_i <- varEst[idx]
    estCovMat_i <- matrix(0, J_i, J_i)
    for (k in 1:J_i) {
      for (l in 1:J_i)
        estCovMat_i[k, l] <- varEst_i[k] ^ 0.5 *
          varEst_i[l] ^ 0.5 * workCorMat[k, l]
    }
    workCovMatList[[i]] <- estCovMat_i
    if (J_i == 1)
      estCovMat_inv2_i <-
      (1 / estCovMat_i) ^ 0.5
    else
      estCovMat_inv2_i <- sqrtm(solve(estCovMat_i))
    workCovMatList_inv2[[i]] <- estCovMat_inv2_i
    y_rw[idx] <- estCovMat_inv2_i %*% Y[idx]
    bsBase_rw[idx,] <- estCovMat_inv2_i %*% as.matrix(bsBase[idx,])
  }


  ########################
  ### Start Screening ####
  ########################
  rssVec <- rep(0, p)
  for (k in 1:p) {
    bsX_k <- t(apply(cbind(X[, k], bsInt), 1, Beta1BsGen))
    bsX_k_rw <- bsX_k
    for (i in 1:n) {
      if (i == 1)
        idx <- 1:(Jcs[i])
      else
        idx <- (Jcs[i - 1] + 1):Jcs[i]
      if (length(idx) == 1)
        bsX_k_rw[idx,] <-
          workCovMatList_inv2[[i]] * bsX_k[idx,]
      else
        bsX_k_rw[idx,] <-
          workCovMatList_inv2[[i]] %*% as.matrix(bsX_k[idx,])
    }
    Xrw_k <- data.frame(bsBase_rw, bsX_k)
    model_k <- lm(y_rw ~ 0 + ., data = Xrw_k)
    rssVec[k] <- sum(model_k$residuals ^ 2)
  }

  results <- list()

  results$error <- rssVec

  results$rank <- rank(rssVec)

  return(results)

}
