#  Copyright (C) 2017-2020 Thorsten Pohlert
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
#' @name osrtTest
#' @aliases osrtTest
#' @title One-Sided Studentized Range Test
#'
#' @description Performs Hayter's one-sided studentized range
#' test against an ordered alternative for normal data
#' with equal variances.
#'
#' @details
#' Hayter's one-sided studentized range test (OSRT) can be used
#' for testing several treatment levels with a zero control in a balanced
#' one-factorial design with normally distributed variables that have a
#' common variance. The null hypothesis, H: \eqn{\mu_i = \mu_j ~~ (i < j)}
#' is tested against a simple order alternative,
#' A: \eqn{\mu_i < \mu_j}, with at least one inequality being strict.
#'
#' The test statistic is calculated as,
#' \deqn{
#'  \hat{h} = \max_{1 \le i < j \le k} \frac{ \left(\bar{x}_j - \bar{x}_i \right)}
#'  {s_{\mathrm{in}} / \sqrt{n}},
#' }{%
#' SEE PDF.
#' }
#'
#' with \eqn{k} the number of groups, \eqn{n = n_1, n_2, \ldots, n_k} and
#' \eqn{s_{\mathrm{in}}^2} the within ANOVA variance. The null hypothesis
#' is rejected, if \eqn{\hat{h} > h_{k,\alpha,v}}, with \eqn{v = N - k}
#' degree of freedom.
#'
#' For the unbalanced case with moderate imbalance the test statistic is
#' \deqn{
#'  \hat{h} = \max_{1 \le i < j \le k} \frac{ \left(\bar{x}_j - \bar{x}_i \right)}
#'  {s_{\mathrm{in}} \sqrt{1/n_j + 1/n_i}},
#' }{%
#' SEE PDF.
#' }
#'
#' The function does not return p-values. Instead the critical h-values
#' as given in the tables of Hayter (1990) for \eqn{\alpha = 0.05} (one-sided)
#' are looked up according to the number of groups (\eqn{k}) and
#' the degree of freedoms (\eqn{v}).
#' Non tabulated values are linearly interpolated with the function
#' \code{\link[stats]{approx}}.
#'
#' @note
#' Hayter (1990) has tabulated critical h-values for balanced designs only.
#' For some unbalanced designs some \eqn{k = 3} critical h-values
#' can be found in Hayter et al. 2001. ' The function will give
#' a warning for the unbalanced case and returns the
#' critical value \eqn{h_{k,\alpha,v} / \sqrt{2}}.
#'
#' @return
#' A list with class \code{"osrt"} that contains the following components:
#' @template returnOsrt
#'
#' @references
#' Hayter, A. J.(1990) A One-Sided Studentised Range
#' Test for Testing Against a Simple Ordered Alternative,
#' \emph{Journal of the American Statistical Association}
#' \bold{85}, 778--785.
#'
#' Hayter, A.J., Miwa, T., Liu, W. (2001)
#' Efficient Directional Inference Methodologies for the
#' Comparisons of Three Ordered Treatment Effects.
#' \emph{J Japan Statist Soc} \bold{31}, 153–174.
#'
#' @keywords htest
#' @concept parametric
#'
#' @seealso
#' \code{link{hayterStoneTest}} \code{\link{MTest}}
#'
#' @importFrom stats var approx
#' @example examples/osrtEx.R
#' @export
osrtTest <- function(x, ...)
    UseMethod("osrtTest")

#' @rdname osrtTest
#' @method osrtTest default
#' @aliases osrtTest.default
#' @template one-way-parms
#' @param alternative the alternative hypothesis. Defaults to \code{greater}.
#' @export
osrtTest.default <-
    function(x, g, alternative = c("greater", "less"), ...)
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

        alternative <- match.arg(alternative)
        if (alternative == "less") {
            x <- -x
        }

        xi <- tapply(x, g, mean)
        ni <- tapply(x, g, length)
        k <- nlevels(g)
        N <- length(x)
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


        val <- numeric(length = k * (k - 1) / 2)
        l <- 0
        if (is.balanced) {
            for (i in 1:(k - 1)) {
                for (j in (i + 1):k) {
                    l <- l + 1
                    val[l] <-  (xi[j] - xi[i]) /
                        (sigma / sqrt(n))
                }
            }
        } else {
            for (i in 1:(k - 1)) {
                for (j in (i + 1):k) {
                    l <- l + 1
                    val[l] <-  (xi[j] - xi[i]) /
                        (sigma * sqrt(1 / ni[j] + 1 / ni[i]))
                }
            }
        }
        STAT <- max(val)

        ## aux function
        hCrit <- approxHayter(k, df)
        hCrit <- adjust.hCrit(hCrit, is.balanced)

        METHOD <- "Hayter's One-Sided Studentized Range Test"
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

#' @rdname osrtTest
#' @aliases osrtTest.formula
#' @method osrtTest formula
#' @template one-way-formula
#' @export
osrtTest.formula <-
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
        y <- do.call("osrtTest", c(as.list(mf),
                                   alternative = alternative))
        y$data.name <- DNAME
        y
    }


#' @rdname osrtTest
#' @aliases osrtTest.aov
#' @method osrtTest aov
## @param x A fitted model object, usually an \link[stats]{aov} fit.
#' @export
osrtTest.aov <- function(x, alternative = c("greater", "less"),
                      ...) {
    model <- x$model
    DNAME <- paste(names(model), collapse = " by ")
    names(model) <- c("x", "g")
    alternative <- match.arg(alternative)
    y <- do.call("osrtTest", c(as.list(model),
                            alternative = alternative))
    y$data.name <- DNAME
    y
}
