## qmandelk.R
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
#' @name Mandel-k
#' @title Mandel's k distribution
#' @description Distribution function and quantile function
#'     for Mandel's k distribution.
#'
#' @seealso \code{\link{mandelkTest}}
#' @aliases qmandelk pmandelk
#' 
#' @references
#' Practice E 691, 2005, \emph{Standard Practice for
#' Conducting an Interlaboratory Study to Determine the
#' Precision of a Test Method}, ASTM International.
#'
#' @importFrom stats qf
#' @param p vector of probabilities.
#' @param k number of groups.
#' @param n number of replicates per group.
#' @param lower.tail logical; if \code{TRUE} (default),
#' probabilities are P[X <= x] otherwise, P[X > x].
#' @param log.p logical; if \code{TRUE}, probabilities
#' are given as log(p).
#' @return
#' \code{pmandelk} gives the distribution function and
#' \code{qmandelk} gives the quantile function.
#' @note The functions are only appropriate for balanced designs.
#' @seealso
#' \code{\link{pmandelh}}, \code{\link{qmandelh}}
#' @keywords distribution
#' @examples
#' qmandelk(0.005, 7, 3, lower.tail=FALSE)
#' @export
qmandelk <- function(p, k, n, lower.tail = TRUE, log.p = FALSE)
{
    fval <- qf(p, df1 = n - 1, df2= ((k-1) * (n-1)),
               lower.tail=lower.tail,
               log.p = log.p)
    kval <- sqrt(k / (1 + (k - 1) / fval))
    return(kval)
}

## taken from package metRology function pmandelk
#' @rdname Mandel-k
#' @param q vector of quantiles.
#' @section Source:
#' The code for \code{pmandelk} was taken from:\cr
#' Stephen L R Ellison. (2017). metRology: Support for Metrological
#' Applications. R package version 0.9-26-2.
#' \url{https://CRAN.R-project.org/package=metRology}
#' @importFrom stats pbeta
#' @keywords distribution
#' @export
pmandelk <- function (q, k, n, lower.tail = TRUE, log.p = FALSE) 
{
    pbeta(q^2/k, (n - 1)/2, (k - 1) * (n - 1)/2, lower.tail = lower.tail, 
        log.p = log.p)
}
