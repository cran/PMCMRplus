# chackoTest.R
# Part of the R package: PMCMR
#
# Copyright (C) 2017 Thorsten Pohlert
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#
#' @name chackoTest
#'
#' @title Testing against Ordered Alternatives (Chacko's Test)
#' 
#' @description
#' Performs Chacko's test for testing against ordered alternatives.
#'
#' @details
#' The null hypothesis, H\eqn{_0: \theta_1 = \theta_2 = \ldots = \theta_k}
#' is tested against a simple order hypothesis,
#' H\eqn{_\mathrm{A}: \theta_1 \le \theta_2 \le \ldots \le
#' \theta_k,~\theta_1 < \theta_k}.
#'
#' The p-values are estimated from the chi-square distribution.
#' 
#' @template class-htest
#' 
#' @template trendTests
#' @section Source:
#' The source code for the application of the pool adjacent violators
#' theorem to calculate the isotonic means
#' was taken from the file \code{"pava.f"}, which is included in the
#' package \pkg{Iso}:
#' 
#'  Rolf Turner (2015). Iso: Functions to Perform Isotonic Regression. R
#'  package version 0.0-17. \url{https://CRAN.R-project.org/package=Iso}.
#' 
#' The file \code{"pava.f"} is a Ratfor modification of Algorithm AS 206.1:
#' 
#' Bril, G., Dykstra, R., Pillers, C., Robertson, T. (1984)
#' Statistical Algorithms: Algorithm AS 206: Isotonic
#' Regression in Two Independent Variables, \emph{Appl. Statist.},
#' 34, 352--357.
#' 
#'  The Algorith AS 206 is available from StatLib
#' \url{http://lib.stat.cmu.edu/apstat}. The Royal Statistical Society
#' holds the copyright to these routines,
#' but has given its permission for their distribution provided that
#' no fee is charged.
#'
#' @references
#' Chacko, V. J. (1963), Testing homogenity against ordered alternatives.
#' \emph{Ann. Math. Statist.}, 34, 945--956.
#' 
#' @importFrom stats pchisq complete.cases
#' @useDynLib 'PMCMRplus', .registration = TRUE, .fixes = "F_"
#' @export
chackoTest <- function(x, ...) UseMethod("chackoTest")

#' @rdname chackoTest
#' @aliases chackoTest.default
#' @method chackoTest default
#' @template one-way-parms
#' @export
chackoTest.default <-
    function(x, g, ...)
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

    rij <- rank(x)
    Ri <- tapply(rij, g, mean, na.rm=T)
    ni <- tapply(rij, g, length)
    k <- nlevels(g)
    N <- sum(ni)
    Rmean <- (N + 1) / 2
    
    ## call to own pava
    Riso <- .Fortran("pava",
                     y=as.double(Ri),
                     w=as.double(ni),
                     kt = integer(k),
                     n = as.integer(k))$y
    
    T <- (12 / (N * (N +1))) * sum(ni * (Riso - Rmean)^2)	   
    df <- k - 1
    PVAL <- pchisq(T, df = df, lower.tail = FALSE)

    METHOD <- paste("Chacko's test")
    ans <- list(method = METHOD, data.name = DNAME, 
                p.value = PVAL,
                statistic = c(H = T), 
                parameter = c(df = df),
                alternative = "greater")
    class(ans) <- "htest"
    ans
}

#' @rdname chackoTest
#' @method chackoTest formula
#' @aliases chackoTest.formula
#' @template one-way-formula
#' @export
chackoTest.formula <-
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
    y <- do.call("chackoTest", c(as.list(mf)))
    y$data.name <- DNAME
    y
}
