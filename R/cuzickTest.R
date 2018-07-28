##  cuzickTest.R
##
##  Copyright (C) 2017, 2018 Thorsten Pohlert
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
#' @name cuzickTest
#' @title Testing against Ordered Alternatives (Cuzick's Test)
#' @description
#' Performs Cuzick's test for testing against ordered alternatives.
#'
#' @details
#' The null hypothesis, H\eqn{_0: \theta_1 = \theta_2 = \ldots = \theta_k}
#' is tested against a simple order hypothesis,
#' H\eqn{_\mathrm{A}: \theta_1 \le \theta_2 \le \ldots \le
#' \theta_k,~\theta_1 < \theta_k}.
#'
#' The p-values are estimated from the standard normal distribution.
#'
#' @template class-htest
#' @template trendTests
#' @references
#' Cuzick, J. (1995) A Wilcoxon-type test for trend, \emph{Statistics in Medicine}
#' \bold{4}, 87--90.
#'
#' @export cuzickTest
cuzickTest <- function(x, ...) UseMethod("cuzickTest")

#' @rdname cuzickTest
#' @method cuzickTest default
#' @aliases cuzickTest.default
#' @template one-way-parms
#' @param alternative the alternative hypothesis. Defaults to \code{"two.sided"}.
#' @param scores numeric vector of scores. Defaults to \code{NULL}.
#' @param continuity logical indicator whether a continuity correction
#' shall be performed. Defaults to \code{FALSE}.
#' @importFrom stats pnorm complete.cases
#' @export
cuzickTest.default <-
function(x, g, alternative = c("two.sided", "greater", "less"),
         scores = NULL, continuity = FALSE, ...)
{
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
        if(!is.null(x$alternative)) alternative <- x$alternative
        if(is.null(x$continuity)) {
            continuity <- FALSE
        } else {
            continuity <- TRUE
        }
        if(is.null(x$scores)) {
            scores <- NULL
        } else {
            scores <- x$scores
        }
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

    alternative <- match.arg(alternative)
    n <- length(x)
    if (n < 2)
        stop("not enough observations")
    ni <- tapply(x, g, length)

    rx <- rank(x)

    ## check for ties
    TIES <- FALSE
    TIES <- sum(table(rx) - 1) != 0
    if(TIES){
        warning("Ties are present. 'cuzickTest()' does not correct for ties")
    }
    Ri <- tapply(rx, g, sum)
    if(is.null(scores)) {
        li <- 1:k
    } else {
        if (!is.numeric(scores) | length(scores) != k){
            stop("'scores' must be a numeric vector with 'length == nr of groups'")
        } else {
            li <- scores
        }
    }
    T <- sum(Ri * li)
    L <- sum(ni * li)
    eT <- L * (n + 1) / 2

    varT <- ((n + 1) / 12) * (n * sum(li^2 * ni) - L^2)

    if (!continuity) {
        STAT <- (T - eT) / sqrt(varT)
    } else {
        if (T > eT){
            STAT <- (T - eT - 0.5) / sqrt(varT)
        } else {
            STAT <- (T - eT + 0.5) / sqrt(varT)
        }
    }
    ESTIMATE <- T

    if (alternative == "two.sided"){
        PVAL <- 2 * min(pnorm(abs(STAT), lower.tail = FALSE), 0.5)
    } else if (alternative == "greater"){
        PVAL <- pnorm(STAT, lower.tail = FALSE)
    } else {
        PVAL <- pnorm(STAT)
    }

    names(STAT) <- "z"
    names(ESTIMATE) <- "CU"

    RVAL <- list(statistic = STAT,
                 p.value = PVAL,
                 method = "Cuzick trend test",
                 data.name = DNAME,
                 alternative = alternative,
                 estimate = ESTIMATE)
    class(RVAL) <- "htest"
    return(RVAL)
}

#' @rdname cuzickTest
#' @method cuzickTest formula
#' @aliases cuzickTest.formula
#' @template one-way-formula
#' @export
cuzickTest.formula <-
    function(formula, data, subset, na.action,
             alternative = c("two.sided", "greater", "less"),
             scores = NULL, continuity = FALSE,  ...)
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
    alternative <- match.arg(alternative)
    if(is.null(scores)) scores <- NULL
    names(mf) <- NULL
    y <- do.call("cuzickTest", c(as.list(mf), alternative = alternative,
                                  scores = scores, continuity = continuity))
    y$data.name <- DNAME
    y
}
