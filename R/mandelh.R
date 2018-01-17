## mandelh.R
## Part of the R package: PMCMRplus
##
## Copyright (C) 2017 Thorsten Pohlert
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
#' @name Mandel-h
#' @title Mandel's h distribution
#' @description Distribution function and quantile function
#'     for Mandel's h distribution.
#'
#' @seealso \code{\link{mandelhTest}}
#' @aliases qmandelh pmandelh
#' 
#' @references
#' Practice E 691, 2005, \emph{Standard Practice for
#' Conducting an Interlaboratory Study to Determine the
#' Precision of a Test Method}, ASTM International.
#'
#' @importFrom stats qt
#' @param p vector of probabilities.
#' @param k number of groups.
#' @param lower.tail logical; if \code{TRUE} (default),
#' probabilities are P[X <= x] otherwise, P[X > x].
#' @param log.p logical; if \code{TRUE}, probabilities
#' are given as log(p).
#' @return
#' \code{pmandelh} gives the distribution function and
#' \code{qmandelh} gives the quantile function.
#'
#' @keywords distribution
#' @examples
#' ## We need a two-sided upper-tail quantile
#' qmandelh(p = 0.005/2, k = 7, lower.tail=FALSE)
#' @export
qmandelh <- function(p, k, lower.tail = TRUE, log.p = FALSE)
{
    tval <- qt(p, df = (k - 2),
               lower.tail=lower.tail,
               log.p = log.p)
    hval <- (k - 1) * tval / sqrt(k * (tval^2 + k - 2 ))
    return(hval)
}

## taken from package metRology function pmandelh
#' @rdname Mandel-h
#' @param q vector of quantiles.
#' @section Source:
#' The code for \code{pmandelh} was taken from:\cr
#' Stephen L R Ellison. (2017). metRology: Support for Metrological
#' Applications. R package version 0.9-26-2.
#' \url{https://CRAN.R-project.org/package=metRology}
#' @importFrom stats pbeta
#' @export
pmandelh <- function (q, k, lower.tail = TRUE, log.p = FALSE) 
{
    pbeta((1 + q * sqrt(k)/(k - 1))/2, (k - 2)/2, (k - 2)/2, 
        lower.tail = lower.tail,
        log.p = log.p)
}
