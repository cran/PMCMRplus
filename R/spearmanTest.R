## spearman.R
##
## Copyright (C) 2017-2020 Thorsten Pohlert
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
#' @name spearmanTest
#' @title Testing against Ordered Alternatives (Spearman Test)
#'
#' @description
#' Performs a Spearman type test for testing against ordered alternatives.
#'
#' @details
#' A one factorial design for dose finding comprises an ordered factor,
#' .e. treatment with increasing treatment levels.
#' The basic idea is to correlate the ranks \eqn{R_{ij}} with the increasing
#' order number \eqn{1 \le i \le k} of the treatment levels (Kloke and McKean 2015).
#' More precisely, \eqn{R_{ij}} is correlated with the expected mid-value ranks
#' under the assumption of strictly increasing median responses.
#' Let the expected mid-value rank of the first group denote \eqn{E_1 = \left(n_1 + 1\right)/2}.
#' The following expected mid-value ranks are
#' \eqn{E_j = n_{j-1} + \left(n_j + 1 \right)/2} for \eqn{2 \le j \le k}.
#' The corresponding number of tied values for the \eqn{i}th group is \eqn{n_i}. #
#' The sum of squared residuals is
#' \eqn{D^2 = \sum_{i=1}^k \sum_{j=1}^{n_i} \left(R_{ij} - E_i \right)^2}.
#' Consequently, Spearman's rank correlation coefficient can be calculated as:
#'
#' \deqn{
#' r_\mathrm{S} = \frac{6  D^2}
#'    {\left(N^3 - N\right)- C},
#' }{%
#'  SEE PDF
#' }
#'
#' with
#' \deqn{
#' C = 1/2 - \sum_{c=1}^r \left(t_c^3 - t_c\right) +
#' 1/2 - \sum_{i=1}^k \left(n_i^3 - n_i \right)
#' }{%
#'  SEE PDF
#' }
#' and \eqn{t_c} the number of ties of the \eqn{c}th group of ties.
#' Spearman's rank correlation coefficient can be tested for
#' significance with a \eqn{t}-test.
#' For a one-tailed test the null hypothesis of \eqn{r_\mathrm{S} \le 0}
#' is rejected and the alternative \eqn{r_\mathrm{S} > 0} is accepted if
#'
#' \deqn{
#' r_\mathrm{S} \sqrt{\frac{\left(n-2\right)}{\left(1 - r_\mathrm{S}\right)}} > t_{v,1-\alpha},
#' }{%
#'  SEE PDF
#' }
#'
#' with \eqn{v = n - 2} degree of freedom.
#'
#' @template class-htest
#' @template trendTests
#' @references
#' Kloke, J., McKean, J. W. (2015) \emph{Nonparametric statistical methods using R}.
#' Boca Raton, FL: Chapman & Hall/CRC.
#'
#' @export spearmanTest
spearmanTest <- function(x, ...) UseMethod("spearmanTest")

#' @rdname spearmanTest
#' @method spearmanTest default
#' @aliases spearmanTest.default
#' @template one-way-parms
#' @param alternative the alternative hypothesis. Defaults to \code{"two.sided"}.
#' @importFrom stats pt
#' @importFrom stats complete.cases
#' @export
spearmanTest.default <-
    function(x, g, alternative = c("two.sided", "greater", "less"),...)
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

    alternative <- match.arg(alternative)
    k <- nlevels(g)
    nk <- tapply(x, g, length)

    n <- length(x)
    gg <- rep(1:k, times=nk)

    ## Function to get ties for tie adjustment
    getties <- function(x){
        t <- table(x)
        C <- sum(t^3 - t)
        C
    }

    ry <- rank(x)
    rx <- rank(gg)

    di <- rx - ry
    Tx <- getties(gg)
    Ty <- getties(x)

    S <- sum(di^2)
    rs <- (n^3 - n - 1/2 * Tx - 1/2 * Ty - 6 * S) /
        sqrt((n^3 - n - Tx) * (n^3 - n - Ty))
    names(rs) <- "rho"
    tval <- rs * sqrt((n - 2) / (1 - rs^2))

    if (alternative == "two.sided"){
        PVAL <- 2 * pt(abs(tval), df=n-2, lower.tail=FALSE)
    } else if (alternative == "greater"){
        PVAL <- pt(tval, df=n-2, lower.tail=FALSE)
    } else {
        PVAL <- pt(tval, df=n-2)
    }

    names(tval) <- "t"
    PARMS <- n-2
    names(PARMS) <- "df"
    ESTIM <- rs
    H0 <- 0
    names(H0) <- "rho"
    METHOD <- paste("Spearman rank correlation test for ordered alternatives")

    ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL,
                statistic = tval, parameter = PARMS, estimate = ESTIM,
                alternative = alternative, null.value = H0)
    class(ans) <- "htest"
    ans
}

#' @rdname spearmanTest
#' @method spearmanTest formula
#' @aliases spearmanTest.formula
#' @template one-way-formula
#' @export
spearmanTest.formula <-
function(formula, data, subset, na.action, alternative = c("two.sided", "greater", "less"),
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
    y <- do.call("spearmanTest", c(as.list(mf), alternative = alternative))
    y$data.name <- DNAME
    y
}
