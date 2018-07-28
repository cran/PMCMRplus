## dscfAllPairsTest.R
## Part of the R package: PMCMRplus
##
## Copyright (C) 2017, 2018 Thorsten Pohlert
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

#' @name dscfAllPairsTest
#' @title Multiple Comparisons of Mean Rank Sums
#' @description
#' Performs the all-pairs comparison test for different factor
#' levels according to Dwass, Steel, Critchlow and Fligner.
#' @details
#' For all-pairs comparisons in an one-factorial layout
#' with non-normally distributed residuals the DSCF
#' all-pairs comparison test can be used. A total of \eqn{m = k(k-1)/2}
#' hypotheses can be tested. The null hypothesis
#' H\eqn{_{ij}: F_i(x) = F_j(x)} is tested in the two-tailed test
#' against the alternative
#' A\eqn{_{ij}: F_i(x) \ne F_j(x), ~~ i \ne j}.
#' As opposed to the all-pairs comparison procedures that depend
#' on Kruskal ranks, the DSCF test is basically an extension of
#' the U-test as re-ranking is conducted for each pairwise test.
#'
#' The p-values are estimated from the studentized range distriburtion.
#'
#' @seealso
#' \code{\link{Tukey}}, \code{\link{pairwise.wilcox.test}}
#'
#' @template class-PMCMR
#'
#' @references
#' Douglas, C. E., Fligner, A. M. (1991) On distribution-free multiple
#' comparisons in the one-way analysis of variance, \emph{Communications in
#'  Statistics - Theory and Methods} \bold{20}, 127--139.
#'
#' Dwass, M. (1960) Some k-sample rank-order tests. In \emph{Contributions to
#'   Probability and Statistics}, Edited by: I. Olkin,
#' Stanford: Stanford University Press.
#'
#' Steel, R. G. D. (1960) A rank sum test for comparing all pairs of
#' treatments, \emph{Technometrics} \bold{2}, 197--207
#'
#' @keywords htest nonparametric
#' @concept RankTransformation
#' @concept AllPairsComparisons
#' @export
dscfAllPairsTest <- function(x, ...) UseMethod("dscfAllPairsTest")

#' @rdname dscfAllPairsTest
#' @aliases dscfAllPairsTest.default
#' @method dscfAllPairsTest default
#' @template one-way-parms
#' @importFrom stats ptukey
#' @importFrom stats pairwise.table
#' @importFrom stats complete.cases
#' @export
dscfAllPairsTest.default <-
function(x, g, ...){
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

    k <- nlevels(g)
    n <- tapply(x, g, length)
    glev <- levels(g)

    ## Function to get ties for tie adjustment
    getties <- function(x){
        t <- table(x)
        C <- sum((t^3 - t) / 12)
        C
    }

    ## function for pairwise comparisons
    compare.stats <-function(i, j){
        m <- n[j]
        nn <- n[i]
        xraw <- c(x[g==glev[i]], x[g==glev[j]])
        rankx <- rank(xraw)
        lev <- c(g[g==glev[i]], g[g==glev[j]])
        R <- tapply(rankx, lev, sum)
        U <- c(m*nn + (m * (m + 1) / 2), m * nn + (nn * (nn + 1) / 2)) - R
        Umn <- min(U)
        S <- m + nn
        VAR <- (m * nn / (S * (S - 1))) * ((S^3 - S) / 12 - getties(rankx))
        PSTAT <- sqrt(2) * (Umn - m * nn / 2) / sqrt(VAR)
        PSTAT
    }

    PSTAT <- pairwise.table(compare.stats,levels(g), p.adjust.method="none")
    DIST <- "q"
    PVAL <- ptukey(abs(PSTAT), nmeans = k, df = Inf, lower.tail = FALSE)
    METHOD <- "Dwass-Steele-Critchlow-Fligner all-pairs test"
    MODEL <- data.frame(x, g)
    ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL,
                statistic = PSTAT, p.adjust.method = "single-step",
                model = MODEL, dist=DIST)
    class(ans) <- "PMCMR"
    ans
}

#' @rdname dscfAllPairsTest
#' @aliases dscfAllPairsTest.formula
#' @method dscfAllPairsTest formula
#' @template one-way-formula
#' @export
dscfAllPairsTest.formula <-
    function(formula, data, subset, na.action, ...)
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
    names(mf) <- NULL
    y <- do.call("dscfAllPairsTest", c(as.list(mf)))
    y$data.name <- DNAME
    y
}
