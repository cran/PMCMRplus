# cochranTest.R
# Part of the R package: PMCMR
#
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
#' @name cochranTest
#' @title Cochran Test
#' @description
#' Performs Cochran's test for testing an outlying (or inlying)
#' variance.
#'
#' @details
#' For normally distributed data the null hypothesis,
#' H\eqn{_0: \sigma_1^2 = \sigma_2^2 = \ldots = \sigma_k^2}
#' is tested against the alternative (greater)
#' H\eqn{_{\mathrm{A}}: \sigma_p > \sigma_i ~~ (i \le k, i \ne p)} with
#' at least one inequality being strict.
#'
#' The p-value is computed with the function \code{\link{pcochran}}.
#'
#' @template class-htest
#' @inherit Cochran references
#'
#' @seealso
#' \code{\link[stats]{bartlett.test}}, \code{\link[stats]{fligner.test}}.
#' @keywords htest
#' @concept outliers
#'
#' @examples
#' data(Pentosan)
#' cochranTest(value ~ lab, data = Pentosan, subset = (material == "A"))
#'
#' @export
cochranTest <- function(x, ...) UseMethod("cochranTest")

#' @rdname cochranTest
#' @method cochranTest default
#' @aliases cochranTest.default
#' @template one-way-parms
#' @param alternative the alternative hypothesis. Defaults to \code{"greater"}
#' @importFrom stats var
#' @export
cochranTest.default <- function(x, g, alternative = c("greater", "less"), ...)
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
        alternative <- x$alternative
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
    Ssq <- tapply(x, g, var)
    n <- length(x)

    if (alternative == "greater"){
        method <- "Cochran test for outlying variance"
        TSsq <- max(Ssq)
        i <- which(TSsq == Ssq)
        names(i) <- NULL
        C <- TSsq / sum(Ssq)
        pval <- pcochran(q = C,
                         k = k,
                         n = n / k,
                         lower.tail = FALSE)

    } else {
        method <- "Cochran test for inlying variance"
        TSsq <- min(Ssq)
        i <- which(TSsq == Ssq)
        names(i) <- NULL
        C <- TSsq / sum(Ssq)
        pval <- pcochran(q = C,
                         k = k,
                         n = n / k,
                         lower.tail = TRUE)
    }

    ans <- list(method = method,
                alternative = alternative,
                statistic = c("C" = C),
                p.value = pval,
                parameter = c("k" = k, "n" = n / k),
                data.name = DNAME,
                estimates = c("group" = i, "var" = TSsq))
    class(ans) <- "htest"
    return(ans)
}

#' @rdname cochranTest
#' @method cochranTest formula
#' @aliases cochranTest.formula
#' @template one-way-formula
#' @export
cochranTest.formula <- function(formula, data, subset, na.action,
         alternative = c("greater", "less"),...)
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
     alternative <- match.arg(alternative)
    y <- do.call("cochranTest", c(as.list(mf), alternative = alternative))
    y$data.name <- DNAME
    y
}
