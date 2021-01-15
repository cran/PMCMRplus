## stepDownTrendTest.R
##
## Copyright (C) 2020 Thorsten Pohlert
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
#' @name stepDownTrendTest
#' @title Step Down Trend Tests
#'
#' @description
#' Performs step-down trend test procedures for monotone responses
#' to detect NOEC (LOEC) according to OECD (2006).
#' @details
#' According to OECD 2006 one can perform a test for trend
#' on responses from all dose groups including the control.
#' If the trend test is significant at the 0.05 level, the
#' high dose group is omitted, and the trend
#' statistic with the remaining dose groups is re-compute
#' The procedure is continued until the trend test is
#' first non-significant at the 0.05 level, then stop.
#'
#' The NOEC is the highest dose
#' remaining at this stage. If this test is significant
#' when only the lowest dose and control remain,
#' then a NOEC cannot be established from the data.
#'
#' @template class-trendPMCMR
#'
#' @references
#' OECD (2006) \emph{Current Approaches in the Statistical
#' Analysis of Ecotoxicity Data: A Guidance to Application},
#' OECD Series on Testing and Assessment \bold{52},
#' Paris: Organisation for Econonomic Co-operation and Development.
#'
#' @examples
#' res <- stepDownTrendTest(Y ~ DOSE, data = trout,
#'                          test = "jonckheereTest",
#'                          alternative = "less")
#' ## print method
#' res
#' ## summary method
#' summary(res)
#' @concept trendtest
#' @keywords htest nonparametric
#' @seealso
#' \code{\link{leTest}}, \code{\link{jonckheereTest}},
#' \code{\link{spearmanTest}}, \code{\link{cuzickTest}},
#' \code{\link{chackoTest}}, \code{\link{johnsonTest}}
#' @export stepDownTrendTest
stepDownTrendTest <- function(x, ...) UseMethod("stepDownTrendTest")

#' @rdname stepDownTrendTest
#' @method stepDownTrendTest default
#' @aliases stepDownTrendTest.default
#' @template one-way-parms
#' @param test the trend test that shall be performed. Defaults to \code{"leTest"}.
#' @param alternative the alternative hypothesis. Defaults to \code{"two.sided"}.
#' @param continuity logical indicator whether a continuity correction
#' shall be performed. Only relevant for \code{"jonckheereTest"}. Defaults to \code{FALSE}.
#' @export
stepDownTrendTest.default <-
    function(x, g, test = c("leTest",
                      "spearmanTest",
                      "jonckheereTest",
                      "cuzickTest",
                      "chackoTest",
                      "johnsonTest"),
             alternative = c("two.sided", "greater", "less"),
             continuity = FALSE,
             ...)
{
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
        if(!is.null(x$alternative)) alternative <- x$alternative
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

    ## check input
    alternative <- match.arg(alternative)
    test <- match.arg(test)


    lev <- levels(g)

    ## sequentially apply the tests
    res <- lapply(2:k, function(i){
        ok <- g %in% lev[1:i]

        inp <- list(alternative = alternative,
                    x = x[ok],
                    g = g[ok],
                    continuity = continuity)

        do.call(test, inp)

    })

    METHOD <- paste0("Step down ", res[[1]]$method)

    pval <- sapply(1:(k-1), function(j) res[[j]]$p.value)
    PVAL <- matrix(pval, nrow = k-1, ncol = 1, byrow = TRUE)
    colnames(PVAL) <- lev[1]
    rownames(PVAL) <- lev[2:k]

    stat <- sapply(1:(k-1), function(j) res[[j]]$statistic)
    STAT <- matrix(stat, nrow = k-1, ncol = 1, byrow = TRUE)
    colnames(STAT) <- lev[1]
    rownames(STAT) <- lev[2:k]

    dist <- names(res[[1]]$statistic)

    ans <- list(method = METHOD,
                data.name = DNAME,
                p.value = PVAL,
                statistic = STAT,
                alternative = alternative,
                p.adjust.method = "none",
                dist = dist)
    class(ans) <- "trendPMCMR"     ### Needs to be a different class
                              ### H0 is imprecise
    ans
}

#' @rdname stepDownTrendTest
#' @method stepDownTrendTest formula
#' @aliases stepDownTrendTest.formula
#' @template one-way-formula
#' @export
stepDownTrendTest.formula <-
function(formula, data, subset, na.action,
         test = c("leTest",
                  "spearmanTest",
                  "jonckheereTest",
                  "cuzickTest",
                  "chackoTest",
                  "johnsonTest"),
         alternative = c("two.sided", "greater", "less"),
         continuity = FALSE,
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
    test <- match.arg(test)
    names(mf) <- NULL
    y <- do.call("stepDownTrendTest", c(as.list(mf),
                                        alternative = alternative,
                                        test = test,
                                        continuity = continuity))
    y$data.name <- DNAME
    y
}
