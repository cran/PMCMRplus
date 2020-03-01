# bwsKSampleTest.R
# Part of the R package: PMCMR
#
##  Copyright (C) 2017 Thorsten Pohlert
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

#' @name bwsKSampleTest
#' @title Murakami's k-Sample BWS Test
#' @description
#' Performs Murakami's k-Sample BWS Test.
#' @details
#' The null hypothesis, H\eqn{_0: F_1 = F_2 = \ldots = F_k}
#' is tested against the alternative,
#' H\eqn{_\mathrm{A}: F_i \ne F_j ~~(i \ne j)}, with at least
#' one unequality beeing strict.
#'
#' The p-values are estimated through an assymptotic boot-strap method.
#'
#' @note
#' One may increase the number of permutations to e.g. \code{nperm = 10000}
#' in order to get more precise p-values. However, this will be on
#' the expense of computational time.
#'
#' @inherit bwsAllPairsTest references
# @references
# Murakami, H. (2006) K-sample rank test based on modified
# Baumgartner statistic and its power comparison,
# \emph{J. Jpn. Comp. Statist.} 19, 1--13.
#'
#' @template class-htest
#'
#' @seealso
#' \code{\link{sample}}, \code{\link{bwsAllPairsTest}},
#' \code{\link{bwsManyOneTest}}.
#'
#' @keywords htest nonparametric
#' @example examples/kSamples.R
#'
#' @export
bwsKSampleTest <- function(x, ...) UseMethod("bwsKSampleTest")

#' @rdname bwsKSampleTest
#' @aliases bwsKSampleTest.default
#' @method bwsKSampleTest default
#' @template one-way-parms
#' @param nperm number of permutations for the assymptotic permutation test.
#' Defaults to \code{1000}.
#' @importFrom stats pt
#' @export
bwsKSampleTest.default <-
function(x, g, nperm=1000, ...){
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

    n <- tapply(x, g, length)
    N <- sum(n)
    k <- nlevels(g)

    ## Murakami, 2007, p. 8
    j <- unlist(sapply(n, function(ni) 1:ni, simplify = "array"))
    nn <- unlist(sapply(n, function(ni) rep(ni, ni), simplify = "array"))

    ERij <- (N + 1) / (nn + 1) * j
    VarRij <- (j / (nn + 1)) * (1 - j / (nn + 1)) *
        ((N - nn) * (N + 1) / (nn + 2))
    bkstar.fn <- function(RR){
        tmp <- tapply((RR - ERij)^2 / VarRij, g, sum)
        Bkstar <- 1/k * sum(1/n * tmp)
        return(Bkstar)
    }

    ## Rank the combined sample and sort in increasing order
    ## per group
    Rij <- rank(x)
    o <- order(as.integer(g), Rij)
    STATISTIC <- bkstar.fn(Rij[o])
    names(STATISTIC) <- NULL

    ## asymptotic bootstrap permutation
    mt <- sapply(1:nperm, function(i) {
        ix <- sample(Rij)
        o <- order(as.integer(g), ix)
        bkstar.fn(ix[o])
        })

    ## pvalues
    PVAL <- length(mt[mt >= STATISTIC]) / nperm

    METHOD <- "Murakami's k-Sample BWS Test"
    ans <- list(statistic = c("Bk*" = STATISTIC),
                parameter = c("k" = k),
                method = METHOD,
                p.value = PVAL,
                data.name = DNAME)
    class(ans) <- "htest"
    ans
}

#' @rdname bwsKSampleTest
#' @aliases bwsKSampleTest.formula
#' @method bwsKSampleTest formula
#' @template one-way-formula
#' @export
bwsKSampleTest.formula <-
    function(formula, data, subset, na.action, nperm=1000,...)
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
    y <- do.call("bwsKSampleTest", c(as.list(mf)))
    y$data.name <- DNAME
    y
}
