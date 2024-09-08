## medianAllPairsTest.R
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
#' @name medianAllPairsTest
#' @title Brown-Mood All Pairs Median Test
#'
#' @description
#' Performs Brown-Mood All Pairs Median Test.
#'
#' @details
#' For all-pairs comparisons in an one-factorial layout
#' with non-normally distributed residuals Brown-Mood
#' non-parametric Median test
#' can be performed. A total of \eqn{m = k(k-1)/2}
#' hypotheses can be tested. The null hypothesis
#' H\eqn{_{ij}: \mu_i(x) = \mu_j(x)} is tested in the two-tailed test
#' against the alternative
#' A\eqn{_{ij}: \mu_i(x) \ne \mu_j(x), ~~ i \ne j}.
#'
#' In this procedure the joined median is used for classification,
#' but pairwise Pearson Chisquare-Tests are conducted. Any method
#' as given by \code{\link[stats]{p.adjust.methods}} can be used
#' to account for multiplicity.
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
#' @template class-PMCMR
#' @keywords nonparametric
#'
#' @seealso
#' \code{\link[stats]{chisq.test}}.
#'
#' @example examples/kwAllPairsMC.R
#'
#' @export medianAllPairsTest
medianAllPairsTest <- function(x, ...)
  UseMethod("medianAllPairsTest")

#' @rdname medianAllPairsTest
#' @method medianAllPairsTest default
#' @aliases medianAllPairsTest.default
#' @template one-way-parms
#' @param p.adjust.method method for adjusting p values
#'    (see \code{\link{p.adjust}}).
#'
#' @importFrom stats complete.cases chisq.test
#' @export
medianAllPairsTest.default <-
  function(x, g, p.adjust.method = p.adjust.methods, ...)
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
    } else {
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


    p.adjust.method <- match.arg(p.adjust.method)

    ## grand median
    med <- median(x)

    ## call internal function
    mat <- buildMedianFreqTable(x, g, med)

    ## sequentially apply Pearson's chisq-test
    PSTAT <- matrix(nrow = k-1, ncol = k-1 )
    PVAL <- PSTAT
    for (i in 1:k-1) {
      for (j in (i+1):k) {
        out <- chisq.test(mat[c(i,j), ],
                          simulate.p.value = FALSE,
                          correct = FALSE)

        PSTAT[j-1,i] <- out$statistic
        PVAL[j-1,i] <- out$p.value
      }
    }

    p <- as.vector(PVAL[lower.tri(PVAL, diag = TRUE)])
    p <- p.adjust(p, method = p.adjust.method)
    PVAL[lower.tri(PVAL, diag = TRUE)] <- p

    MODEL <- data.frame(x = x, g = g)
    DIST <- "Chisquare"
    METHOD <- paste("Brown-Mood Median Test")
    colnames(PVAL) <- levels(g)[-k]
    rownames(PVAL) <- levels(g)[-1]
    colnames(PSTAT) <- colnames(PVAL)
    rownames(PSTAT) <- rownames(PVAL)

    ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL,
                statistic = PSTAT, p.adjust.method = p.adjust.method,
                model = MODEL, dist = DIST, alternative = "two.sided")
    class(ans) <- "PMCMR"
    ans
}

#' @rdname medianAllPairsTest
#' @method medianAllPairsTest formula
#' @aliases medianAllPairsTest.formula
#' @template one-way-formula
#' @export
medianAllPairsTest.formula <-
  function(formula, data, subset, na.action,
           p.adjust.method = p.adjust.methods, ...)
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
    p.adjust.method <- match.arg(p.adjust.method)
    y <- do.call("medianAllPairsTest", c(as.list(mf),
                                 p.adjust.method = p.adjust.method))
    y$data.name <- DNAME
    y
  }
