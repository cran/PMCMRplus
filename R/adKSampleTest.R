## adKSampleTest.R
##
## Copyright (C) 2017-2020 Thorsten Pohlert
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
#' @name adKSampleTest
#' @title Anderson-Darling k-Sample Test
#'
#' @description
#' Performs Anderson-Darling k-sample test.
#'
#' @details
#' The null hypothesis, H\eqn{_0: F_1 = F_2 = \ldots = F_k}
#' is tested against the alternative,
#' H\eqn{_\mathrm{A}: F_i \ne F_j ~~(i \ne j)}, with at least
#' one unequality beeing strict.
#'
#' This function only evaluates version 1 of the k-sample Anderson-Darling
#' test (i.e. Eq. 6) of Scholz and Stephens (1987).
#' The p-values are estimated with the extended empirical function
#' as implemented in \code{\link[kSamples]{ad.pval}} of
#' the package \pkg{kSamples}.
#'
#' @template class-htest
#' @keywords nonparametric
#'
#' @inherit adAllPairsTest references
#'
#' @seealso
#' \code{\link{adAllPairsTest}}, \code{\link{adManyOneTest}},
#' \code{\link[kSamples]{ad.pval}}.
#' @example examples/kSamples.R
#'
#' @export adKSampleTest
adKSampleTest <- function(x, ...)
  UseMethod("adKSampleTest")

#' @rdname adKSampleTest
#' @method adKSampleTest default
#' @aliases adKSampleTest.default
#' @template one-way-parms
#' @importFrom stats complete.cases
# @importFrom kSamples ad.pval
#' @importFrom kSamples ad.test
#' @export
adKSampleTest.default <-
  function(x, g, ...)
  {
    ## taken from stats::kruskal.test

    if (is.list(x)) {
      if (length(x) < 2L)
        stop("'x' must be a list with at least 2 elements")
      DNAME <- deparse(substitute(x))
      x <- lapply(x, function(u)
        u <- u[complete.cases(u)])
      k <- length(x)
      l <- sapply(x, "length")
      if (any(l == 0))
        stop("all groups must contain data")
      g <- factor(rep(1:k, l))
      x <- unlist(x)
    }
    else {
      if (length(x) != length(g))
        stop("'x' and 'g' must have the same length")
      DNAME <- paste(deparse(substitute(x)), "and",
                     deparse(substitute(g)))
      OK <- complete.cases(x, g)
      x <- x[OK]
      g <- g[OK]
      if (!all(is.finite(g)))
        stop("all group levels must be finite")
      g <- factor(g)
      k <- nlevels(g)
      if (k < 2)
        stop("all observations are in the same group")
    }

    # ix <- order(as.character(g))
    # g <- g[ix]
    # x <- x[ix]
    #
    # ## prepare
    # n <- tapply(x, g, length)
    # N <- sum(n)

    ## Function ADstatV1 with ties
    ## This could be put into a Fortan routine
    ## in order to increase speed.
    # Zstar <- unique(sort(x))
    # L <- length(Zstar)
    # idx <- c(which(!duplicated(g)), N+1)
    # stat <- 0.0
    #
    # print(idx)
    #
    # AsqkN <- .Fortran("adstatv1",
    #          x = as.double(x),
    #          n = as.integer(N),
    #          zstar = as.double(Zstar),
    #          l = as.integer(L),
    #          idx = as.integer(idx),
    #          k = as.integer(k),
    #          stat = as.double(stat))$stat
    #
    # print(AsqkN)
#
#     ADstatV1 <- function(x, g) {
#       lev <- levels(g)
#       Zstar <- sort(x)
#       Zstar <- unique(Zstar)
#       L <- length(Zstar)
#       f <- matrix(0, ncol = L, nrow = k)
#       for (i in 1:k) {
#         tmp <- x[g == lev[i]]
#         f[i, ] <- sapply(Zstar, function(Z) {
#           return(length(tmp[tmp == Z]))
#         })
#       }
#
#       l <- sapply(1:L, function(j)
#         sum(f[, j]))
#
#       print(l)
#
#       tmp <- rep(NA, k)
#       for (i in 1:k) {
#         tm <- 0
#         for (j in 1:(L - 1)) {
#           Bj <- sum(l[1:j])
#           Mij <- sum(f[i, 1:j])
#           tm <- tm + (l[j] / N) *
#             ((N * Mij - n[i] * Bj) ^ 2 /
#                (Bj * (N - Bj)))
#         }
#         tmp[i] <- tm
#       }
#       AsqkN <- sum(1 / n * tmp)
#       return(AsqkN)
#     }
#
#     AsqkN <- ADstatV1(x, g)
#
#     ## Calculate varAsqkN
#     ## Scholz and Stephens, 1987, Eq. 4
#     H <- sum(1 / n)
#
#     h <- 0
#     for (i in seq_len(N - 1)) {
#       h <- h + 1 / i
#     }
#
#     G <- 0
#     for (i in seq_len(N - 2)) {
#       for (j in (i + 1):(N - 1)) {
#         G <- G + 1 / ((N - i) * j)
#
#       }
#
#     }
#
#     ## coefficients
#     a <- (4 * G - 6) * (k - 1) + (10 - 6 * G) * H
#     b <- (2 * G - 4) * k ^ 2 + 8 * h * k +
#       (2 * G - 14 * h - 4) * H - 8 * h + 4 * G - 6
#     c <- (6 * h + 2 * G - 2) *  k ^ 2 +
#       (4 * h - 4 * G + 6) * k +
#       (2 * h - 6) * H + 4 * h
#     d <- (2 * h + 6) * k ^ 2 - 4 * h * k
#
#     varAsqkN <- (a * N ^ 3 + b * N ^ 2 + c * N + d) /
#       ((N - 1) * (N - 2) * (N - 3))
#
      m <- k - 1
#     TkN <- (AsqkN - m) / sqrt(varAsqkN)
#
#     ## p-value from package kSample
#     PVAL <- ad.pval(tx = TkN, m = m, version = 1)

    METHOD <- paste("Anderson-Darling k-Sample Test")

    ##
    out <- ad.test(split(x, g), method = "asymptotic")


    ans <- list(
      method = METHOD,
      data.name = DNAME,
      p.value = out$ad[2,3],
      statistic = c(TkN = out$ad[2,2]),
      parameter = c(m = m)
    #  estimate = c(A2kN = out$ad[2,1])
    )
    class(ans) <- "htest"
    ans
  }

#' @rdname adKSampleTest
#' @method adKSampleTest formula
#' @aliases adKSampleTest.formula
#' @template one-way-formula
#' @export
adKSampleTest.formula <-
  function(formula, data, subset, na.action, ...)
  {
    mf <- match.call(expand.dots = FALSE)
    m <-
      match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(stats::model.frame)

    if (missing(formula) || (length(formula) != 3L))
      stop("'formula' missing or incorrect")
    mf <- eval(mf, parent.frame())
    if (length(mf) > 2L)
      stop("'formula' should be of the form response ~ group")
    DNAME <- paste(names(mf), collapse = " by ")
    names(mf) <- NULL
    y <- do.call("adKSampleTest", as.list(mf))
    y$data.name <- DNAME
    y
  }
