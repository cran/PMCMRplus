## Dgrubbs.R
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
#' @name Dgrubbs
#' @title Grubbs D* distribution
#' @description Distribution function for Grubbs D* distribution.
#'
#' @aliases pdgrubbs
#'
#' @references
#' Grubbs, F.E. (1950) Sample criteria for testing outlying observations,
#' \emph{Ann. Math. Stat.} \bold{21}, 27--58.
#'
#' Wilrich, P.-T. (2011) Critical values of Mandel's h and k,
#' Grubbs and the Cochran test statistic, \emph{Adv. Stat. Anal.}.
#' \doi{10.1007/s10182-011-0185-y}.
#'
#' @param q vector of quantiles.
#' @param n total sample size.
#' @param m number of Monte-Carlo replicates. Defaults to \code{10,000}.
#' @param lower.tail  logical; if TRUE (default), probabilities are P[X <= x]
#'          otherwise, P[X > x].
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @return
#' \code{pgrubbs} gives the distribution function
#' @seealso
#' \code{\link{Grubbs}}
#' @keywords distribution
#' @examples
#' pdgrubbs(0.62, 7, 1E4)
#' @useDynLib 'PMCMRplus', .registration = TRUE, .fixes = "F_"
#' @importFrom stats na.omit
#' @export
pdgrubbs <- function(q, n, m = 1E4, lower.tail = TRUE, log.p = FALSE)
{
    q <- na.omit(q)
    n <- na.omit(n)
    nn <- length(q)
    if (length(n) != nn){
        n <- rep(n, nn)
    }
    m <- rep(m, nn)

    p <- sapply(1:nn, function(i)
        .Fortran("pd",
                 q = as.double(q[i]),
                 n = as.integer(n[i]),
                 m = as.integer(m[i]),
                 p = as.double(0))$p)

    if (lower.tail){
        p <- 1 - p
    }
    if (log.p){
        p <- log(p)
    }
    return(p)
}
