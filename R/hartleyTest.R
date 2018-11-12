# hartleyTest.R
# Part of the R package: PMCMR
#
# Copyright (C) 2018 Thorsten Pohlert
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
#   Uses pmaxFratio of package SuppDists
#

#' @name hartleyTest
#' @title Hartley's Maximum F-Ratio Test of Homogeneity of
#' Variances
#' @description
#'  Performs Hartley's maximum F-ratio test of the null that
#'  variances in each of the groups (samples) are the same.
#'
#' @details
#' If \code{x} is a list, its elements are taken as the samples
#' to be compared for homogeneity of variances.  In this
#' case, the elements must all be numeric data vectors,
#' \code{g} is ignored, and one can simply use
#' \code{hartleyTest(x)} to perform the test.  If the samples are not
#' yet contained in a list, use \code{hartleyTest(list(x, ...))}.
#'
#' Otherwise, \code{x} must be a numeric data vector, and \code{g} must
#' be a vector or factor object of the same length as \code{x} giving the
#' group for the corresponding elements of \code{x}.
#'
#' Hartley's parametric test requires normality and
#' a nearly balanced design. The p-value of the test
#' is calculated with the function \code{\link[SuppDists]{pmaxFratio}}
#' of the package \pkg{SuppDists}.
#'
#' @references
#' Hartley, H.O. (1950) The maximum F-ratio
#' as a short cut test for heterogeneity of variance,
#' \emph{Biometrika} \bold{37}, 308--312.
#'
#' @seealso
#'  \code{\link{bartlett.test}}, \code{\link[SuppDists]{pmaxFratio}}
#'
#' @template class-htest
#'
#' @keywords htest
#' @examples
#' hartleyTest(count ~ spray, data = InsectSprays)
#'
#' @export
hartleyTest <- function(x, ...) UseMethod("hartleyTest")

#' @rdname hartleyTest
#' @method hartleyTest default
#' @template one-way-parms
#' @importFrom SuppDists pmaxFratio
#' @export
hartleyTest.default <-
    function(x,
             g,
             ...)
{
        ## taken from stats::hartleyTest
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


    ##
    Var <- tapply(x, g, var)
    n <- tapply(x, g, length)
    N <- mean(n)
    k <- nlevels(g)

    if (any(n != N)) {
      warning("Maximum F-ratio test is imprecise for unbalanced designs.")
    }
    mxId <- which.max(Var)
    mnId <- which.min(Var)

    ## Choose n of minimum variance,
    ## if n is equal, than thiswill not have any effect
    df <- n[mnId] - 1

    PSTAT <- Var[mxId] / Var[mnId]
    names(PSTAT) <- "F Max"
    PARMS <- c(df, k)
    names(PARMS) <- c("df", "k")

    PVAL <- pmaxFratio(PSTAT, df = df, k = k, lower.tail = FALSE)

    METHOD <-
      paste("Hartley's maximum F-ratio test of homogeneity of variances")

    ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL,
                statistic = PSTAT, parameter = PARMS)
    class(ans) <- "htest"
    ans
}

#' @rdname hartleyTest
#' @method hartleyTest formula
#' @template one-way-formula
#' @export
hartleyTest.formula <-
    function(formula, data, subset, na.action,
              ...)
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
    y <- do.call("hartleyTest", c(as.list(mf)))
    y$data.name <- DNAME
    y
}
