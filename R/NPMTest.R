## NPMTest.R
## Part of the R package: PMCMRplus
##
## Copyright (C) 2017, 2018 Thorsten Pohlert
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## A copy of the GNU General Public License is available at
## http://www.r-project.org/Licenses/

#' @title All-Pairs Comparisons for Simply Ordered Mean Ranksums
#'
#' @description
#' Performs Nashimoto and Wright's all-pairs comparison procedure
#' for simply ordered mean ranksums.
#' Their test denoted as NPM test is basically an
#' extension of Nemenyi's procedure for testing
#' increasingly ordered alternatives.
#'
#' The modified procedure uses the property of a simple order,
#' \eqn{\theta_m' - \theta_m \le \theta_j - \theta_i \le \theta_l' - \theta_l
#' \qquad (l \le i \le m~\mathrm{and}~ m' \le j \le l')}.
#' The null hypothesis H\eqn{_{ij}: \theta_i = \theta_j} is tested against
#' the alternative A\eqn{_{ij}: \theta_i < \theta_j} for any
#' \eqn{1 \le i < j \le k}.
#'
#' The p-values are estimated from the studentized range distribution.
#' If the medians are already increasingly ordered, than the NPM-test simplifies
#' to the ordinary Nemenyi test (see \code{\link{kwAllPairsNemenyiTest}}).
#'
#' @name NPMTest
# @inherit chaAllPairsNashimotoTest
# @references
# Nashimoto, K., Wright, F.T., (2005), Multiple comparison procedures
# for detecting differences in simply ordered means.
# \emph{Comput. Statist. Data Anal.} 48, 291--306.
#'
#' @keywords htest nonparametric
#'
#' @seealso
#' \code{\link{kwAllPairsNemenyiTest}}
#' @template class-PMCMR
#' @examples
#' ## Example from Sachs (1997, p. 402)
#' x <- c(106, 114, 116, 127, 145,
#'        110, 125, 143, 148, 151,
#'        136, 139, 149, 160, 174)
#' g <- gl(3,5)
#' levels(g) <- c("A", "B", "C")
#' NPMTest(x, g)
#'
#' @export
NPMTest <- function(x, ...) UseMethod("NPMTest")

#' @rdname NPMTest
#' @aliases NPMTest.default
#' @method NPMTest default
#' @template one-way-parms
#' @importFrom stats complete.cases
#' @importFrom stats ptukey
#' @export
NPMTest.default <-
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

    ## prepare tukey test
    rij <- rank(x)
    Ri <- tapply(rij, g, mean)
    ni <- tapply(x, g, length)
    k <- nlevels(g)
    n <- length(x)
    df <- Inf

    sigma <- sqrt(n * (n + 1) / 12)

    STAT <- matrix(NA, ncol=k-1, nrow=k-1)
    for (i in 1:(k-1)){
        for(j in (i+1):k){
            u <- j
            m <- i:(u-1)
            tmp <- sapply(m, function(m) (Ri[u] - Ri[m]) /
                                         (sigma / sqrt(2) *
                                          sqrt(1 / ni[m] + 1 /ni[u])))
            STAT[j-1,i] <- max(tmp)
        }
    }

    colnames(STAT) <- levels(g)[1:(k-1)]
    rownames(STAT) <- levels(g)[2:k]

    PVAL <- ptukey(STAT, nmeans = k, df = df, lower.tail=FALSE)

    colnames(PVAL) <- colnames(STAT)
    rownames(PVAL) <- rownames(STAT)

    METHOD <- "Nashimoto-Wright NPM-Test for ordered means \n\t\t of non-normal data"
    MODEL <- data.frame(x, g)
    DIST <- "q"

    ans <- list(method = METHOD,
                data.name = DNAME,
                p.value = PVAL,
                statistic = STAT,
                p.adjust.method = "single-step",
                model = MODEL,
                dist = DIST,
                alternative = "greater")
    class(ans) <- "PMCMR"
    ans
}

#' @rdname NPMTest
#' @aliases NPMTest.formula
#' @method NPMTest formula
#' @template one-way-formula
#' @export
NPMTest.formula <-
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
    y <- do.call("NPMTest", c(as.list(mf)))
    y$data.name <- DNAME
    y
}
