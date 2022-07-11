## Dunnett.R
## Part of the R package: PMCMRplus
##
## Copyright (C) 2022 Thorsten Pohlert
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  A copy of the GNU General Public License is available at
##  http://www.r-project.org/Licenses/
##
##
#' @title Dunnett Distribution
#' @rdname Dunnett
#' @aliases Dunnett
#' @aliases qDunnett
#'
#' @description
#' Distribution function and quantile function
#' for the distribution of Dunnett's many-to-one
#' comparisons test.
#'
#' @param p vector of probabilities.
#' @param n0 sample size for control group.
#' @param n vector of sample sizes for treatment groups.
#'
#' @details
#' Dunnett's distribution is a special case of the
#' multivariate t distribution.
#'
#' Let the total sample size be \eqn{N = n_0 + \sum_i^m n_i}, with \eqn{m} the
#' number of treatment groups, than the quantile \eqn{T_{m v \rho \alpha}}
#' is calculated with \eqn{v = N - k} degree of freedom and
#' the correlation \eqn{\rho}
#'
#' \deqn{
#'  \rho_{ij} = \sqrt{\frac{n_i n_j}
#'              {\left(n_i + n_0\right) \left(n_j+ n_0\right)}} ~~
#'              (i \ne j).
#' }
#'
#' The functions determines \eqn{m} via the length of the input
#' vector \code{n}.
#'
#' Quantiles and p-values are computed with the functions
#' of the package **mvtnorm**.
#'
#' @return
#' \code{pDunnett} gives the distribution function and
#' \code{qDunnett} gives its inverse, the quantile function.
#'
#' @seealso
#' \code{\link[mvtnorm]{qmvt}} \code{\link[mvtnorm]{pmvt}} \code{\link{dunnettTest}}
#'
#' @note
#' The results are seed depending.
#'
#' @examples
#' ## Table gives 2.34 for df = 6, m = 2, one-sided
#' set.seed(112)
#' qval <- qDunnett(p = 0.05, n0 = 3, n = rep(3,2))
#' round(qval, 2)
#' set.seed(112)
#' pDunnett(qval, n0=3, n = rep(3,2), lower.tail = FALSE)
#'
#' ## Table gives 2.65 for df = 20, m = 4, two-sided
#' set.seed(112)
#' qval <- qDunnett(p = 0.05/2, n0 = 5, n = rep(5,4))
#' round(qval, 2)
#' set.seed(112)
#' 2 * pDunnett(qval, n0= 5, n = rep(5,4), lower.tail= FALSE)
#' @export
qDunnett <- function(p, n0, n) {



  N <- sum(n) + n0
  groups <- length(n) + 1
  df <- N - groups

  ## call to internal function
  rho <- dunnettCorrMat(n0, n)

  pp <-  1 - p
  #' @importFrom mvtnorm qmvt
  out <- sapply(pp, function(prop) {
    qmvt(p = prop,
              tail = "lower.tail",
              df = df,
              corr = rho)$quantile
  })

  out
}

#' @rdname Dunnett
#' @aliases pDunnett
#' @param q vector of quantiles.
#' @param lower.tail logical; if TRUE (default),
#' probabilities are \eqn{P[X \leq x]} otherwise, \eqn{P[X > x]}.
#' @export
pDunnett <- function(q, n0, n, lower.tail = TRUE) {

  N <- sum(n) + n0
  groups <- length(n) + 1
  m <- groups - 1
  df <- N - groups

  ## call to internal function
  rho <- dunnettCorrMat(n0, n)

  if (lower.tail) {
    pval <- sapply(q,
                   function(qq) {
                     #' @importFrom mvtnorm pmvt
                     1 - pmvt(
                       lower = rep(qq, m),
                       upper = Inf,
                       df = df,
                       corr = rho
                     )
                   })

  } else {
    pval <- sapply(q,
                   function(qq) {
                     #' @importFrom mvtnorm pmvt
                     1 - pmvt(
                       lower = -Inf,
                       upper = rep(qq, m),
                       df = df,
                       corr = rho
                     )
                   })

  }

  pval

}

## internal function
## builds Dunnett's correlation matrix
## for later use in qmvt or pmvt
dunnettCorrMat <- function(n0, n) {

  groups <- length(n) + 1
  m <- groups - 1

  ## correlation matrix
  rho <- diag(1, m, m)
  for (i in 1:(m-1)) {
    for (j in (i+1):m) {
      rho[i,j] <- ((n[i] * n[j]) /
                     ((n[i] + n0) * (n[j] + n0)))^(1/2)
      rho[j,i] <- rho[i,j]

    }
  }

  rho
}
