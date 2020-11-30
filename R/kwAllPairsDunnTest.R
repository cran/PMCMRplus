## kwAllPairsDunn.test.R
## Part of the R package: PMCMR
##
## Copyright (C) 2015-2018 Thorsten Pohlert
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

#' @name kwAllPairsDunnTest
#' @title Dunn's All-Pairs Rank Comparison Test
#'
#' @description
#' Performs Dunn's non-parametric all-pairs comparison test
#' for Kruskal-type ranked data.
#'
#' @details
#' For all-pairs comparisons in an one-factorial layout
#' with non-normally distributed residuals Dunn's non-parametric test
#' can be performed. A total of \eqn{m = k(k-1)/2}
#' hypotheses can be tested. The null hypothesis
#' H\eqn{_{ij}: \mu_i(x) = \mu_j(x)} is tested in the two-tailed test
#' against the alternative
#' A\eqn{_{ij}: \mu_i(x) \ne \mu_j(x), ~~ i \ne j}.
#'
#' The p-values are computed from the standard normal distribution using
#' any of the p-adjustment methods as included in \code{\link{p.adjust}}.
#' Originally, Dunn (1964) proposed Bonferroni's p-adjustment method.
#'
#' @references
#' Dunn, O. J. (1964) Multiple comparisons using rank sums,
#' \emph{Technometrics} \emph{6}, 241--252.
#'
#' Siegel, S., Castellan Jr., N. J. (1988) \emph{Nonparametric Statistics
#'  for The Behavioral Sciences}. New York: McGraw-Hill.
#'
#' @template class-PMCMR
#' @keywords nonparametric
#' @concept kruskalranks
#'
#' @seealso
#' \code{\link[stats]{Normal}}, \code{\link[stats]{p.adjust}},
#' \code{\link{kruskalTest}},
#' \code{\link{kwAllPairsConoverTest}}, \code{\link{kwAllPairsNemenyiTest}}
#'
#' @example examples/kwAllPairsMC.R
#' @export
kwAllPairsDunnTest <- function(x, ...) UseMethod("kwAllPairsDunnTest")

#' @rdname kwAllPairsDunnTest
#' @method kwAllPairsDunnTest default
#' @aliases kwAllPairsDunnTest default
#' @template one-way-parms
#' @param p.adjust.method method for adjusting p values
#'    (see \code{\link{p.adjust}}).
#' @importFrom stats pnorm
#' @importFrom stats p.adjust.methods
#' @importFrom stats p.adjust
#' @importFrom stats complete.cases
#' @export
kwAllPairsDunnTest.default <-
function(x, g, p.adjust.method = p.adjust.methods, ...){
    ## taken from stats::kruskal.test

    if (is.list(x)) {
        if (length(x) < 2L)
            stop("'x' must be a list with at least 2 elements")
        DNAME <- deparse(substitute(x))
        x <- lapply(x, function(u) u <- u[complete.cases(u)])
        k <- length(x)
        l <- sapply(x, "length")
        if (any(l == 0))
            stop("all groups must contain data")
        g <- factor(rep(1 : k, l))
        p.adjust.method <- x$p.adjust.method
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

    p.adjust.method <- match.arg(p.adjust.method)
    x.rank <- rank(x)
    R.bar <- tapply(x.rank, g, mean,na.rm=T)
    R.n <- tapply(!is.na(x), g, length)
    g.unique <- unique(g)
    k <- length(g.unique)
    n <- sum(R.n)

    METHOD <- "Dunn's all-pairs test"

    ## get the ties
    C <- gettiesDunn(x.rank)

    if (C != 0) warning("Ties are present. z-quantiles were corrected for ties.")
    compare.stats <- function(i,j) {
        dif <- abs(R.bar[i] - R.bar[j])
        A <- n * (n+1) / 12
        B <- (1 / R.n[i] + 1 / R.n[j])
        zval <- dif / sqrt((A - C) * B)
        return(zval)
    }
    PSTAT <- pairwise.table(compare.stats,levels(g), p.adjust.method="none" )
    compare.levels <- function(i,j) {
        dif <- abs(R.bar[i] - R.bar[j])
        A <- n * (n+1) / 12
        B <- (1 / R.n[i] + 1 / R.n[j])
        zval <- dif / sqrt((A - C) * B)
        pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
        return(pval)
    }
    PVAL <- pairwise.table(compare.levels,levels(g),
                           p.adjust.method=p.adjust.method )
    MODEL <- data.frame(x, g)
    DIST <- "z"
    ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL,
                statistic = PSTAT, p.adjust.method = p.adjust.method,
                model = MODEL, dist = DIST, alternative = "two.sided")
    class(ans) <- "PMCMR"
    ans
}

#' @rdname kwAllPairsDunnTest
#' @method kwAllPairsDunnTest formula
#' @aliases kwAllPairsDunnTest.formula
#' @template one-way-formula
#' @export
kwAllPairsDunnTest.formula <-
function(formula, data, subset, na.action,
         p.adjust.method = p.adjust.methods, ...)
{
    mf <- match.call(expand.dots=FALSE)
    m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(stats::model.frame)

   if(missing(formula) || (length(formula) != 3L))
        stop("'formula' missing or incorrect")
    mf <- eval(mf, parent.frame())
    if(length(mf) > 2L)
       stop("'formula' should be of the form response ~ group")
    DNAME <- paste(names(mf), collapse = " by ")
    p.adjust.method <- match.arg(p.adjust.method)
    names(mf) <- NULL
    y <- do.call("kwAllPairsDunnTest", c(as.list(mf),
                                     p.adjust.method))
    y$data.name <- DNAME
    y
}
