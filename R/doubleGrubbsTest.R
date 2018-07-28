# doubleGrubbsTest.R
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
#' @name doubleGrubbsTest
#' @title Grubbs Double Outlier Test
#' @description
#' Performs Grubbs double outlier test.
#'
#' @details
#' Let \eqn{X} denote an identically and independently distributed continuous
#' variate with realizations \eqn{x_i ~~ (1 \le i \le k)}.
#' Further, let the increasingly ordered realizations
#' denote \eqn{x_{(1)} \le x_{(2)} \le \ldots \le x_{(n)}}. Then
#' the following model for testing two maximum outliers can be proposed:
#'
#' \deqn{
#'   x_{(i)} = \left\{
#'       \begin{array}{lcl}
#'        \mu + \epsilon_{(i)}, & \qquad & i = 1, \ldots, n - 2 \\
#'        \mu + \Delta + \epsilon_{(j)} & \qquad & j = n-1, n  \\
#'       \end{array} \right.}{%
#'   x[(i)] = \mu + \epsilon[(i)] for i = 1, ..., n - 2
#' and x[(i)] = \mu + \Delta + \epsilon[(j)] for j = n-1, n}
#'
#' with \eqn{\epsilon \approx N(0,\sigma)}. The null hypothesis,
#' H\eqn{_0: \Delta = 0} is tested against the alternative,
#' H\eqn{_{\mathrm{A}}: \Delta > 0}.
#'
#' For testing two minimum outliers, the model can be proposed
#' as
#'
#' \deqn{
#'   x_{(i)} = \left\{
#'       \begin{array}{lcl}
#'        \mu + \Delta + \epsilon_{(j)} & \qquad & j = 1, 2 \\
#'        \mu + \epsilon_{(i)}, & \qquad & i = 3, \ldots, n \\
#'       \end{array} \right.}{%
#'   x[(i)] = \mu + \Delta + \epsilon[(j)] for j = 1, 2
#'   and x[(i)] = \mu + \epsilon[(i)] for i = 3, ..., n}
#'
#' The null hypothesis is tested against the alternative,
#' H\eqn{_{\mathrm{A}}: \Delta < 0}.
#'
#' The p-value is computed with the function \code{\link{pdgrubbs}}.
#' @param x a numeric vector of data.
#' @param alternative the alternative hypothesis.
#' Defaults to \code{"two.sided"}.
#' @param m number of Monte-Carlo replicates.
#' @template class-htest
#' @inherit Grubbs references
#' @keywords htest
#' @concept outliers
#' @examples
#' data(Pentosan)
#' dat <- subset(Pentosan, subset = (material == "A"))
#' labMeans <- tapply(dat$value, dat$lab, mean)
#' doubleGrubbsTest(x = labMeans, alternative = "less")
#' @importFrom stats na.omit
#' @useDynLib 'PMCMRplus', .registration = TRUE, .fixes = "F_"
#' @export
doubleGrubbsTest <- function(x, alternative = c("two.sided", "greater", "less"), m = 1E4)
{
    alternative <- match.arg(alternative)
    DNAME <- deparse(substitute(x))
    x <- na.omit(x)
    n <- length(x)
    d <- as.double(0)
    if (alternative != "two.sided"){
        if (alternative == "less"){
            o <- order(x)
            nn <- c(1,2)
            y <- x[o]
            g <- unlist(
                sapply(nn , function(i) which(y[i] == x))
                )
            estimates <- c(c("i" = g),
                           c("value"= x[g]))
            x <- -x
        } else if (alternative == "greater") {
            o <- order(x)
            nn <- c(n-1,n)
            y <- x[o]
            g <- unlist(
                sapply(nn, function(i) which(y[i] == x))
                )
            estimates <- c(c("i" = g),
                           c("value" = x[g]))
        }

        ## Get statistics
        D <- .Fortran("dstat",
                      x=as.double(x),
                      d=as.double(d),
                      n=as.integer(n))$d
    } else {
        meanx <- mean(x)
        xx <- abs(x - meanx)
        o <- order(xx)
        nn <- c(n-1,n)
        y <- x[o]
        g <- unlist(
            sapply(nn, function(i) which(y[i] == x))
            )
        estimates <- c(c("i" = g),
                       c("value" = x[g]))

        Q <- (x - meanx)^2
        yy <- y[-nn]
        meany <- mean(yy)
        Q2 <- (yy - meany)^2
        D <- sum(Q2) / sum(Q)
    }

    ## Get pvalue
    pval <- pdgrubbs(D, n, m, lower.tail = TRUE)

    ans <- list(statistic = c("D*" = D),
                p.value = pval,
                data.name = DNAME,
                method = "Grubbs test for double outliers",
                alternative = alternative,
                null.value = c("Delta" = 0),
                estimates = estimates)
    class(ans) <- "htest"
    return(ans)
}
