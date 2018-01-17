##  vanWaerden.test.R
##
##  Copyright (C) 2015, 2016, 2017 Thorsten Pohlert
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
#' @title van der Waerden's Normal Scores Test
#' @description
#'  Performs van der Waerden's normal scores test.
#' @details
#' For one-factorial designs with non-normally distributed
#' residuals van der Waerden's normal scores test can be performed to test
#' the H\eqn{_0: F_1(x) = F_2(x) = \ldots = F_k(x)} against
#' the H\eqn{_\mathrm{A}: F_i (x) \ne F_j(x)~ (i \ne j)} with at least
#' one strict inequality.
#' @name vanWaerdenTest
#' @note
#' A tie correction is not applied in this function.
#'
#' @references
#' W. J. Conover and R. L. Iman (1979)
#' \emph{On multiple-comparisons procedures},
#' Tech. Rep. LA-7677-MS, Los Alamos Scientific Laboratory.
#' 
#' B. L. van der Waerden (1952) Order tests for the two-sample
#' problem and their power, \emph{Indagationes Mathematicae}, 14, 453--458.
#'
#' @seealso
#' \code{\link{kruskalTest}}, \code{\link{normalScoresTest}} 
#' @examples
#' vanWaerdenTest(count ~ spray, data = InsectSprays)
#' 
#' @keywords htest nonparametric
#' @export
vanWaerdenTest <- function(x, ...) UseMethod("vanWaerdenTest")

#' @rdname vanWaerdenTest
#' @method vanWaerdenTest default
#' @aliases vanWeardenTest.default
#' @importFrom stats pchisq
#' @importFrom stats qnorm
#' @template one-way-parms
#' @export
vanWaerdenTest.default <-
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

    n <- length(x)
    if (n < 2)
        stop("not enough observations")
    r <- rank(x)

    # transform to z-scores
    zscores <- qnorm(r / (n+1))
    AJ <- tapply(zscores, g, sum)
    NJ <- tapply(zscores, g, length)
    s2 <- sum(zscores^2) / (n - 1)
    STATISTIC <- (1 / s2) * sum(AJ^2 / NJ)
    PARAMETER <- k - 1

    PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
    names(STATISTIC) <- "Van der Waerden chi-squared"
    names(PARAMETER) <- "df"

    RVAL <- list(statistic = STATISTIC,
                 parameter = PARAMETER,
                 p.value = PVAL,
                 method = "Van der Waerden normal scores test",
                 data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
}

#' @rdname vanWaerdenTest
#' @method vanWaerdenTest formula
#' @aliases vanWaerdenTest.formula
#' @template one-way-formula
#' @export
vanWaerdenTest.formula <-
function(formula, data, subset, na.action, ...)
{
    if(missing(formula) || (length(formula) != 3L))
        stop("'formula' missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if(is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)
    m[[1L]] <- quote(stats::model.frame)
    mf <- eval(m, parent.frame())
    if(length(mf) > 2L)
        stop("'formula' should be of the form response ~ group")
    DNAME <- paste(names(mf), collapse = " by ")
    names(mf) <- NULL
    y <- do.call("vanWaerdenTest", as.list(mf))
    y$data.name <- DNAME
    y
}
