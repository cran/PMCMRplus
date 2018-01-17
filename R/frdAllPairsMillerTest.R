## frdAllPairsMillerTest.R
## Part of the R package: PMCMR
##
## Copyright (C)  2017 Thorsten Pohlert
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

#' @rdname frdAllPairsMillerTests
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
#' J. Bortz J, G. A. Lienert, K. Boehnke (1990), \emph{Verteilungsfreie
#'    Methoden in der Biostatistik}. Berlin: Springer.
#' 
#' R. G. Miller Jr. (1996), \emph{Simultaneous statistical inference}.
#'  New York: McGraw-Hill.
#'
#' E. L. Wike (2006), \emph{Data Analysis. A Statistical Primer for
#'    Psychology Students}. New Brunswick: Aldine Transaction.
#' 
#' @keywords htest nonparametric
#' @concept Friedman
#' @concept Rank
#' @concept AllPairs
#' @seealso
#' \code{\link{friedmanTest}}, \code{\link{friedman.test}},
#' \code{\link{frdAllPairsExactTest}}, \code{\link{frdAllPairsConoverTest}},
#' \code{\link{frdAllPairsNemenyiTest}}, \code{\link{frdAllPairsSiegelTest}}
#' @template class-PMCMR
#' @example examples/frdAllPairs.R
#' @export
frdAllPairsMillerTest <- function(y, ...) UseMethod("frdAllPairsMillerTest")

#' @rdname frdAllPairsMillerTests
#' @aliases frdAllPairsMillerTest.default
#' @method frdAllPairsMillerTest default
#' @template two-way-parms
#' @importFrom stats pchisq
#' @export
frdAllPairsMillerTest.default <-
    function(y, groups, blocks, ...)
{
    if ((is.matrix(y)) | (is.data.frame(y))) {
        # corrected 4. Jun 2017
        DNAME <- paste(deparse(substitute(y)))
        GRPNAMES <- colnames(y)
        k <- length(GRPNAMES)
        BLOCKNAMES <- rownames(y)
        n <- length(BLOCKNAMES)
        groups <- factor(rep(GRPNAMES, times = n))
        blocks <- factor(rep(BLOCKNAMES, each = k))
        y <- as.vector(t(y))
    }
    else {
        if (any(is.na(groups)) || any(is.na(blocks))) 
            stop("NA's are not allowed in groups or blocks")
        if (any(diff(c(length(y), length(groups), length(blocks))))) 
            stop("y, groups and blocks must have the same length")
        if (any(table(groups, blocks) != 1)) 
            stop("Not an unreplicated complete block design")
            
        DNAME <- paste(deparse(substitute(y)), ",",
                       deparse(substitute(groups)), "and",
                       deparse(substitute(blocks)) )
        groups <- factor(groups)
        blocks <- factor(blocks)
        k <- nlevels(groups)
        n <- nlevels(blocks)
        GRPNAMES <- as.character(groups[1:k])
    }
    
    mat <- matrix(y, nrow = n, ncol = k, byrow = TRUE)
    r <- t(apply(mat, 1L, rank))
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
                            levels(groups), p.adjust.method="none")
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
    MODEL <- data.frame(x = y, g = groups, b = blocks)
    ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL,
                statistic = PSTAT, p.adjust.method = p.adjust.method,
                dist = DIST, model = MODEL, parameter = PARAMETER)
    class(ans) <- "PMCMR"
    ans
}
