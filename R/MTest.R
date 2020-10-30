## MTest.R
## Part of the R package: PMCMRplus
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

#' @name MTest
#' @title Extended One-Sided Studentised Range Test
#'
#' @description Performs Nashimoto-Wright's extended
#' one-sided studentised range
#' test against an ordered alternative for normal data
#' with equal variances.
#'
#' @details
#' The procedure uses the property of a simple order,
#' \eqn{\theta_m' - \mu_m \le \mu_j - \mu_i \le \mu_l' - \mu_l
#' \qquad (l \le i \le m~\mathrm{and}~ m' \le j \le l')}.
#' The null hypothesis H\eqn{_{ij}: \mu_i = \mu_j} is tested against
#' the alternative A\eqn{_{ij}: \mu_i < \mu_j} for any
#' \eqn{1 \le i < j \le k}.
#'
#' The all-pairs comparisons test statistics for a balanced design are
#' \deqn{
#'  \hat{h}_{ij} = \max_{i \le m < m' \le j} \frac{\left(\bar{x}_{m'} - \bar{x}_m \right)}
#'  {s_{\mathrm{in}} / \sqrt{n}},
#' }{%
#'  SEE PDF
#' }
#'
#' with \eqn{n = n_i; ~ N = \sum_i^k n_i ~~ (1 \le i \le k)}, \eqn{\bar{x}_i} the arithmetic mean of the \eqn{i}th group,
#' and \eqn{s_{\mathrm{in}}^2} the within ANOVA variance. The null hypothesis is rejected,
#' if \eqn{\hat{h} > h_{k,\alpha,v}}, with \eqn{v = N - k}
#' degree of freedom.
#'
#' For the unbalanced case with moderate imbalance the test statistic is
#' \deqn{
#'  \hat{h}_{ij} = \max_{i \le m < m' \le j} \frac{\left(\bar{x}_{m'} - \bar{x}_m \right)}
#'  {s_{\mathrm{in}} \left(1/n_m + 1/n_{m'}\right)^{1/2}},
#' }{%
#'  SEE PDF
#' }
#'
#' The null hypothesis is rejected, if \eqn{\hat{h}_{ij} > h_{k,\alpha,v} / \sqrt{2}}.
#'
#' The function does not return p-values. Instead the critical h-values
#' as given in the tables of Hayter (1990) for \eqn{\alpha = 0.05} (one-sided)
#' are looked up according to the number of groups (\eqn{k}) and
#' the degree of freedoms (\eqn{v}).
#'
#' @template returnOsrt
#'
#' @note
#' The function will give a warning for the unbalanced case and returns the
#' critical value \eqn{h_{k,\alpha,\infty} / \sqrt{2}}.
#'
#' @references
#' Hayter, A. J.(1990) A One-Sided Studentised Range
#' Test for Testing Against a Simple Ordered Alternative,
#' \emph{Journal of the American Statistical Association}
#' \bold{85}, 778--785.
#'
#' Nashimoto, K., Wright, F.T., (2005) Multiple comparison procedures
#' for detecting differences in simply ordered means.
#' \emph{Comput. Statist. Data Anal.} \bold{48}, 291--306.
#'
#' @seealso
#' \code{\link{osrtTest}}, \code{\link{NPMTest}}
#'
#' @keywords htest
#' @concept parametric
#' @importFrom stats ptukey
#' @example examples/osrtEx.R
#' @export
MTest <- function(x, ...)
    UseMethod("MTest")

#' @rdname MTest
#' @aliases MTest.default
#' @method MTest default
#' @param alternative the alternative hypothesis. Defaults to \code{greater}.
#' @template one-way-parms
#' @importFrom stats var approx
#' @importFrom stats complete.cases
#' @export
MTest.default <-
    function(x, g, alternative = c("greater", "less"), ...) {
        ## taken from stats::kruskal.test

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

        ## check alternative
        alternative <- match.arg(alternative)
        if (alternative == "less") {
            x <- -x
        }

        ## prepare osrt test
        ni <- tapply(x, g, length)
        N <- sum(ni)
        xi <- tapply(x, g, mean)
        s2i <- tapply(x, g, var)
        df <- N - k
        s2in <- 1 / df * sum(s2i * (ni - 1))
        sigma <- sqrt(s2in)

        ## balanced design
        n <- ni[1]
        ## check for all equal
        ok <- sapply(2:k, function(i)
            ni[i] == n)
        is.balanced <- all(ok)

        STAT <- matrix(NA, ncol = k - 1, nrow = k - 1)

        if (is.balanced) {
            for (i in 1:(k - 1)) {
                for (j in (i + 1):k) {
                    u <- j
                    m <- i:(u - 1)
                    tmp <- sapply(m, function(m)
                        (xi[u] - xi[m]) /
                            (sigma / sqrt(n)))
                    STAT[j - 1, i] <- max(tmp)
                }
            }
        } else {
            for (i in 1:(k - 1)) {
                for (j in (i + 1):k) {
                    u <- j
                    m <- i:(u - 1)
                    tmp <- sapply(m, function(m)
                        (xi[u] - xi[m]) /
                            (sigma * sqrt(1 / ni[m] + 1 / ni[u])))
                    STAT[j - 1, i] <- max(tmp)
                }
            }
        }
        colnames(STAT) <- levels(g)[1:(k - 1)]
        rownames(STAT) <- levels(g)[2:k]

        ## interpolate with aux function
        hCrit <- approxHayter(k, df)
        hCrit <- adjust.hCrit(hCrit, is.balanced)

        METHOD <- "Pairwise comparisons using Nashimoto-Wright's M-Test"
        MODEL <- data.frame(x, g)

        parameter = c(k, df)
        names(parameter) <- c("k", "df")

        ans <- list(
            method = METHOD,
            data.name = DNAME,
            crit.value = hCrit,
            statistic = STAT,
            parameter = parameter,
            alternative = alternative,
            dist = "h"
        )
        class(ans) <- "osrt"
        ans
    }

#' @rdname MTest
#' @aliases MTest.formula
#' @method MTest formula
#' @template one-way-formula
#' @export
MTest.formula <-
    function(formula,
             data,
             subset,
             na.action,
             alternative = c("greater", "less"),
             ...)
    {
        mf <- match.call(expand.dots = FALSE)
        m <-
            match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
        mf <- mf[c(1L, m)]
        mf[[1L]] <- quote(stats::model.frame)

        if (missing(formula) || (length(formula) != 3L))
            stop("'formula' missing or incorrect")
        mf <- eval(mf, parent.frame())
        if (length(mf) > 2L)
            stop("'formula' should be of the form response ~ group")
        DNAME <- paste(names(mf), collapse = " by ")
        alternative <- match.arg(alternative)
        names(mf) <- NULL
        y <- do.call("MTest", c(as.list(mf), alternative = alternative))
        y$data.name <- DNAME
        y
    }


#' @rdname MTest
#' @aliases MTest.aov
#' @method MTest aov
## @param x A fitted model object, usually an \link[stats]{aov} fit.
#' @export
MTest.aov <- function(x, alternative = c("greater", "less"),
                      ...) {
    model <- x$model
    DNAME <- paste(names(model), collapse = " by ")
    names(model) <- c("x", "g")
    alternative <- match.arg(alternative)
    y <- do.call("MTest", c(as.list(model),
                 alternative = alternative))
    y$data.name <- DNAME
    y
}
