#  Copyright (C) 2020 Thorsten Pohlert
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
#' @name steelsKSampleTest
#' @aliases steelsKSampleTest
#' @title Steel's k-Treatments vs. Control Test
#'
#' @description Performs the non-parametric Steel's test
#' for simultaneously testing k-treatments vs. one control.
#'
#' @details
#' It tests \eqn{H: F(i) = F(0), ~ i \le k}, against
#' \eqn{A: F(i) > F(0)} (greater) with at least one inequality being strict.
#'
#' The function is a wrapper function that calls \code{Steel.test} of
#' the package \pkg{kSamples} with argument \code{method = "asymptotic"}.
#'
#' @template class-htest
#'
#' @references
#' Scholz, F. and Zhu, A. (2019). kSamples: K-Sample Rank Tests and
#' their Combinations. R package version 1.2-9.
#' \url{https://CRAN.R-project.org/package=kSamples}
#'
#' Steel, R. G. D. (1959) A Multiple Comparison Rank Sum Test:
#' Treatments Versus Control, \emph{Biometrics} \bold{15}, 560--572.
#'
#' @seealso
#' \code{\link[kSamples]{Steel.test}}, \code{\link{flignerWolfeTest}}
#'
#' @keywords htest nonparametric
#'
#' @examples
#' ## Example from Sachs (1997, p. 402)
#' x <- c(106, 114, 116, 127, 145,
#' 110, 125, 143, 148, 151,
#' 136, 139, 149, 160, 174)
#' g <- gl(3,5)
#' levels(g) <- c("0", "I", "II")
#'
#' ## Steel's Test
#' steelsKSampleTest(x ~ g, alternative = "greater")
#'
#'
#' ## Example from USEPA (2002):
#' ## Reproduction data from a Ceriodaphnia dubia
#' ## 7-day chronic test to several concentrations
#' ## of effluent. Dose level 50% is excluded.
#' x <- c(20, 26, 26, 23, 24, 27, 26, 23, 27, 24,
#' 13, 15, 14, 13, 23, 26, 0, 25, 26, 27,
#' 18, 22, 13, 13, 23, 22, 20, 22, 23, 22,
#' 14, 22, 20, 23, 20, 23, 25, 24, 25, 21,
#' 9, 0, 9, 7, 6, 10, 12, 14, 9, 13,
#' rep(0,10))
#' g <- gl(6, 10)
#' levels(g) <- c("Control", "3%", "6%", "12%", "25%", "50%")
#'
#' ## NOEC at 3%, LOEC at 6%
#' steelsKSampleTest(x ~ g, subset = g != "50%", alternative = "less")
#'
#' @importFrom kSamples Steel.test
#' @export
steelsKSampleTest <- function(x, ...) UseMethod("steelsKSampleTest")

#' @rdname steelsKSampleTest
#' @method steelsKSampleTest default
#' @aliases steelsKSampleTest.default
#' @template one-way-parms
#' @param alternative the alternative hypothesis. Defaults to \code{two.sided}.
#' @export
steelsKSampleTest.default <-
    function(x, g, alternative = c("two.sided", "greater", "less"), ...)
    {
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

        method <- "asymptotic"
        alternative <- match.arg(alternative)

        ## spelling in Steel.test is "two-sided".
        newAltern <- gsub(pattern = "[.]",
                            replacement = "-",
                            alternative)

        lev <- levels(g)
        y <- lapply(lev, function(i)
            x[g == i])

        ## try with do.call
        tmp <- do.call("Steel.test",
                       list(
                           x = y,
                           g = lev,
                           method = method,
                           alternative = newAltern
                       ))

       S <- tmp$st[1]
       names(S) <- "S"

        ans <- list(
            statistic = S,
            p.value = tmp$st[2],
            method = "Steel's k-Treatments vs. Control Test",
            data.name = DNAME,
            alternative = alternative
        )
        class(ans) <- "htest"
        return(ans)
    }

#' @rdname steelsKSampleTest
#' @aliases steelsKSampleTest.formula
#' @method steelsKSampleTest formula
#' @template one-way-formula
#' @export
steelsKSampleTest.formula <-
    function(formula, data, subset, na.action,
             alternative = c("two.sided", "greater", "less"),
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
    alternative <- match.arg(alternative)
    names(mf) <- NULL
    y <- do.call("steelsKSampleTest", c(as.list(mf),
                               alternative = alternative))
    y$data.name <- DNAME
    y
}
