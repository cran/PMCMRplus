## qcochran.R
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
#' @name Cochran
#' @title Cochran's distribution
#' @description Distribution function and quantile function
#'     for Cochran's distribution.
#'
#' @aliases qcochran
#' 
#' @references
#' Cochran, W.G. (1941) The distribution of the largest of a set of estimated
#' variances as a fraction of their total. \emph{Ann. Eugen.} 11, 47--52.
#'
#' Wilrich, P.-T. (2011) Critical values of Mandel's h and k,
#' Grubbs and the Cochran test statistic. \emph{Adv. Stat. Anal.}.
#' \url{http://dx.doi.org/10.1007/s10182-011-0185-y}.
#' 
#' @param p vector of probabilities.
#' @param k number of groups.
#' @param n (average) sample size of the k groups.
#' @param lower.tail logical; if TRUE (default),
#' probabilities are P[X <= x] otherwise, P[X > x].
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @return
#' \code{pcochran} gives the distribution function and
#' \code{qcochran} gives the quantile function.
#' @seealso
#' \code{\link{FDist}}
#' @keywords distribution
#' @importFrom stats qf
#' @examples
#' qcochran(0.05, 7, 3)
#' @export
qcochran <- function(p, k, n, lower.tail = TRUE, log.p = FALSE)
{
    if (log.p){
        p <- exp(p)
    }
    
    if (!lower.tail){
        pp <- (1 - p) / k
    } else {
        pp <- p / k
    }
 
    C <- 1 / (1 + (k - 1) * qf(p = pp,
                               df1=(k-1) * (n-1),
                               df2 = n-1,
                               lower.tail = TRUE,
                               log.p = FALSE)
    )
    return(C)
}

#' @rdname Cochran
#' @aliases pcochran
#' @param q vector of quantiles.
#' @importFrom stats pf
#' @keywords distribution
#' @export
pcochran <-function (q, k, n, lower.tail = TRUE, log.p = FALSE) 
{
    qF <- (1 / q - 1) / (k - 1)
    pval <- 1 - k * pf(q = qF,
                       df1 = (k - 1) * (n - 1),
                       df2 = n - 1,
                       lower.tail=TRUE)
    
    if (!lower.tail) {
        pval <- 1 - pval
    }
    if (log.p){
        pval <- log(pval)
    }
    return(pval)
}
