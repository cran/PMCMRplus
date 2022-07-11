## qgrubbs.R
# Part of the R package: PMCMRplus
#
# Copyright (C) 2017 Thorsten Pohlert
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
#' @name Grubbs
#' @title Grubbs distribution
#' @description Distribution function and quantile function
#'     for Grubbs distribution.
#'
#' @aliases qgrubbs pgrubbs
#'
#' @references
#' Grubbs, F. E. (1950) Sample criteria for testing outlying observations.
#' \emph{Ann. Math. Stat.} \bold{21}, 27--58.
#'
#' Wilrich, P.-T. (2011) Critical values of Mandel's h and k,
#' Grubbs and the Cochran test statistic. \emph{Adv. Stat. Anal.}.
#' \doi{10.1007/s10182-011-0185-y}.
#'
#' @param p vector of probabilities.
#' @param n total sample size.
#' @return
#' \code{pgrubbs} gives the distribution function and
#' \code{qgrubbs} gives the quantile function.
#' @seealso
#' \code{\link{TDist}}
#' @keywords distribution
#' @importFrom stats qt
#' @examples
#' qgrubbs(0.05, 7)
#' @export
qgrubbs <- function(p, n)
{
    tq <- qt(p= p/n,
             df = n-2,
             lower.tail=FALSE,
             log.p = FALSE)
    q <-  ((n - 1) * tq) / sqrt(n * (n - 2 + tq^2))
    return(q)
}

#' @rdname Grubbs
#' @param q vector of quantiles.
#' @param lower.tail logical; if TRUE (default),
#' probabilities are \eqn{P[X \leq x]} otherwise, \eqn{P[X > x]}.
#' @importFrom stats pbeta
#' @keywords distribution
#' @export
pgrubbs <- function (q, n, lower.tail = TRUE)
{
    pval <- n * pbeta((1 + abs(q) * sqrt(n)/(n - 1))/2, (n - 2)/2, (n - 2)/2,
                      lower.tail = FALSE,
                      log.p = FALSE)
    pval <- min(1, pval)
    if (!lower.tail){
        return(pval)
    } else {
        return(1 - pval)
    }
}
