## frdAllPairsMillerTest.R
## Part of the R package: PMCMR
##
## Copyright (C)  2017-2020 Thorsten Pohlert
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

#' @name frdAllPairsMillerTest
#' @title Millers's All-Pairs Comparisons Test for Unreplicated Blocked Data
#' @description
#'  Performs Miller's all-pairs comparisons tests of Friedman-type ranked data.
#' @details
#' For all-pairs comparisons in a two factorial unreplicated
#' complete block design
#' with non-normally distributed residuals, Miller's test can be
#' performed on Friedman-type ranked data.
#'
#' A total of \eqn{m = k ( k -1 )/2} hypotheses can be tested.
#' The null hypothesis, H\eqn{_{ij}: \theta_i = \theta_j}, is tested
#' in the two-tailed case against the alternative,
#' A\eqn{_{ij}: \theta_i \ne \theta_j, ~~ i \ne j}.
#'
#' The \eqn{p}-values are computed from the chi-square distribution.
#'
#' @references
#' Bortz J., Lienert, G. A., Boehnke, K. (1990) \emph{Verteilungsfreie
#'    Methoden in der Biostatistik}. Berlin: Springer.
#'
#' Miller Jr., R. G. (1996) \emph{Simultaneous statistical inference}.
#'  New York: McGraw-Hill.
#'
#' Wike, E. L. (2006), \emph{Data Analysis. A Statistical Primer for
#'    Psychology Students}. New Brunswick: Aldine Transaction.
#'
#' @keywords htest nonparametric
#' @concept friedmanranks
#' @seealso
#' \code{\link{friedmanTest}}, \code{\link{friedman.test}},
#' \code{\link{frdAllPairsExactTest}}, \code{\link{frdAllPairsConoverTest}},
#' \code{\link{frdAllPairsNemenyiTest}}, \code{\link{frdAllPairsSiegelTest}}
#' @template class-PMCMR
#' @example examples/frdAllPairs.R
#' @export
frdAllPairsMillerTest <- function(y, ...) UseMethod("frdAllPairsMillerTest")

#' @rdname frdAllPairsMillerTest
#' @aliases frdAllPairsMillerTest.default
#' @method frdAllPairsMillerTest default
#' @template two-way-parms
#' @importFrom stats pchisq
#' @export
frdAllPairsMillerTest.default <-
    function(y, groups, blocks, ...)
{


    ## 2019-10-16
    ## novel external function
    ans <- frdRanks(y, groups, blocks)
    r <- ans$r
    n <- nrow(r)
    k <- ncol(r)
    GRPNAMES <- colnames(r)

    R.mnsum <- colMeans(r)

    compare.stats <- function(i,j) {
        dif <- abs(R.mnsum[i] - R.mnsum[j])
        qval <- dif / sqrt(k * (k + 1) / (6 * n))
        return(qval)
    }

    ## chisquare
    METHOD <- c("Miller, Bortz et al. and Wike all-pairs test",
                " for a two-way",
                " balanced complete block design")
    PSTAT <- pairwise.table(compare.stats,
                            GRPNAMES,
                            p.adjust.method="none")
    PSTAT <- PSTAT^2
    PVAL <- pchisq(PSTAT, df = (k-1), lower.tail = FALSE)
    DIST <- "X"
    p.adjust.method <- "none"

    PARAMETER <- k - 1
    names(PARAMETER) <- "df"
    colnames(PSTAT) <- GRPNAMES[1:(k-1)]
    rownames(PSTAT) <- GRPNAMES[2:k]
    colnames(PVAL) <- GRPNAMES[1:(k-1)]
    rownames(PVAL) <- GRPNAMES[2:k]

    ans <- list(method = METHOD,
                data.name = ans$DNAME,
                p.value = PVAL,
                statistic = PSTAT,
                p.adjust.method = p.adjust.method,
                dist = DIST,
                model = ans$inDF,
                parameter = PARAMETER)
    class(ans) <- "PMCMR"
    ans
}
