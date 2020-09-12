# grubbsTest.R
# Part of the R package: PMCMRplus
#
# Copyright (C) 2017, 2018 Thorsten Pohlert
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

#' @name grubbsTest
#' @title Grubbs Outlier Test
#' @description
#' Performs Grubbs single outlier test.
#'
#' @details
#' Let \eqn{X} denote an identically and independently distributed continuous
#' variate with realizations \eqn{x_i ~~ (1 \le i \le k)}.
#' Further, let the increasingly ordered realizations
#' denote \eqn{x_{(1)} \le x_{(2)} \le \ldots \le x_{(n)}}. Then
#' the following model for a single maximum outlier can be proposed:
#'
#' \deqn{
#'   x_{(i)} = \left\{
#'       \begin{array}{lcl}
#'        \mu + \epsilon_{(i)}, & \qquad & i = 1, \ldots, n - 1 \\
#'        \mu + \Delta + \epsilon_{(n)} & & \\
#'       \end{array} \right.}{%
#'   x[(i)] = \mu + \epsilon[(i)] for i = 1, ..., n - 1
#' and x[(i)] = \mu + \Delta + \epsilon[(n)]}
#'
#' with \eqn{\epsilon \approx N(0,\sigma)}. The null hypothesis,
#' H\eqn{_0: \Delta = 0} is tested against the alternative,
#' H\eqn{_{\mathrm{A}}: \Delta > 0}.
#'
#' For testing a single minimum outlier, the model can be proposed
#' as
#'
#' \deqn{
#'   x_{(i)} = \left\{
#'       \begin{array}{lcl}
#'        \mu + \Delta + \epsilon_{(1)} & & \\
#'        \mu + \epsilon_{(i)}, & \qquad & i = 2, \ldots, n \\
#'       \end{array} \right.}{%
#'   x[(i)] = \mu + \Delta + \epsilon[(1)]
#'   and x[(i)] = \mu + \epsilon[(i)] for i = 2, ..., n}
#'
#' The null hypothesis is tested against the alternative,
#' H\eqn{_{\mathrm{A}}: \Delta < 0}.
#'
#' The p-value is computed with the function \code{\link{pgrubbs}}.
#' @param x a numeric vector of data.
#' @param alternative the alternative hypothesis.
#' Defaults to \code{"two.sided"}.
#' @template class-htest
#' @inherit Grubbs references
#' @keywords htest univariate
#' @concept outliers
#' @examples
#' data(Pentosan)
#' dat <- subset(Pentosan, subset = (material == "A"))
#' labMeans <- tapply(dat$value, dat$lab, mean)
#' grubbsTest(x = labMeans, alternative = "two.sided")
#' @importFrom stats na.omit
#' @export
grubbsTest <- function(x, alternative = c("two.sided", "greater", "less")){

    if (!is.numeric(x)){
        stop("'x' must be a numeric vector.")
    }
    alternative <- match.arg(alternative)
    data.name <- deparse(substitute(x))
    x <- na.omit(x)
    Mean <- mean(x)
    S <- sd(x)
    n <- length(x)

    if (alternative == "two.sided"){
        MxAbs <- max(abs(x - Mean))
        G <- MxAbs / S
        pval <- min(1, 2 * pgrubbs(G, n, lower.tail=FALSE))
        i <- which(MxAbs == abs(x-Mean))
        val <- x[i]

    } else if (alternative == "greater"){
        val <- max(x)
        G <- (val - Mean) / S
        pval <- pgrubbs(G, n, lower.tail=FALSE)
        i <- which(val == x)

    } else {
        val <- min(x)
        G <- (Mean - val) / S
        pval <- pgrubbs(G, n, lower.tail = FALSE)
        i <- which(val == x)

    }

    names(val) <- NULL
    names(i) <- NULL

    ans <- list(
        method = "Grubbs single outlier test",
        alternative = alternative,
        statistic = c("G" = G),
        parameter = c("df" = n-2),
        p.value = pval,
        estimate = c(c("i" = i),
                     c("value" = val)),
        data.name = data.name)
    class(ans) <- "htest"
    return(ans)
}
