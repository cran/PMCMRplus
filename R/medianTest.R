## medianTest.R
##
## Copyright (C) 2023 Thorsten Pohlert
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
#' @name medianTest
#' @title Brown-Mood Median Test
#'
#' @description
#' Performs Brown-Mood Median Test.
#'
#' @details
#' The null hypothesis, H\eqn{_0: \theta_1 = \theta_2 =
#' \ldots = \theta_k}
#' is tested against the alternative,
#' H\eqn{_\mathrm{A}: \theta_i \ne \theta_j ~~(i \ne j)}, with at least
#' one unequality beeing strict.
#'
#' @return
#' A list with class \sQuote{htest}. For details see
#' \code{\link[stats]{chisq.test}}.
#'
#' @keywords nonparametric
#'
#' @references
#' Brown, G.W., Mood, A.M., 1951,
#' On Median Tests for Linear Hypotheses,
#' in: \emph{Proceedings of the Second Berkeley Symposium on
#' Mathematical Statistics and Probability}.
#' University of California Press, pp. 159–167.
#'
#' @seealso
#' \code{\link[stats]{chisq.test}}.
#'
#' @example examples/kSamples.R
#'
#' @export medianTest
medianTest <- function(x, ...)
  UseMethod("medianTest")

#' @rdname medianTest
#' @method medianTest default
#' @aliases medianTest.default
#' @template one-way-parms
#' @param simulate.p.value a logical indicating whether to compute
#' p-values by Monte-Carlo simulation.
#' @param B an integer specifying the number of replicates used
#' in the Monte-Carlo test.
#'
#' @importFrom stats complete.cases chisq.test
#' @export
medianTest.default <-
  function(x, g, simulate.p.value = FALSE,
           B = 2000, ...)
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

    ## grand median
    med <- median(x)

    mat <- buildMedianFreqTable(x, g, med)

    ## apply Pearson's chisq-test
    out <- chisq.test(mat,
                      simulate.p.value = simulate.p.value,
                      B = B,
                      correct = FALSE)

    if (simulate.p.value) {
      METHOD <- "Brown-Mood Median Test with simulated p-values"
    } else {
      METHOD <- paste("Brown-Mood Median Test")
    }

    out$method <- METHOD
    return(out)
  }

#' @rdname medianTest
#' @method medianTest formula
#' @aliases medianTest.formula
#' @template one-way-formula
#' @export
medianTest.formula <-
  function(formula, data, subset, na.action,
           simulate.p.value = FALSE,
           B = 2000, ...)
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
    y <- do.call("medianTest", c(as.list(mf),
                                 simulate.p.value = simulate.p.value,
                                 B = B))
    y$data.name <- DNAME
    y
  }

## internal
buildMedianFreqTable <- function(x, g, med){
  ## observation per group
  n <- tapply(x, g, length)

  ## count x > median
  fct <- function(x, med) {
    counts <- 0
    for (i in seq_along(x)) {
      if (x[i] > med) {
        counts <- counts + 1
      }
    }
    return(counts)
  }
  gtMed <- tapply(x, g, fct, med)
  leqMed <- n - gtMed

  ## actual frequency table
  mat <- matrix(c(gtMed, leqMed), ncol = 2, byrow = FALSE)
  colnames(mat) <- c("gtMedian", "leqMedian")
  rownames(mat) <- levels(g)

  return(mat)
}
