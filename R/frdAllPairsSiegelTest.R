## frdAllPairsSiegelTest.R
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

#' @rdname frdAllPairsSiegelTest
#' @title Siegel and Castellan's All-Pairs Comparisons Test for
#' Unreplicated Blocked Data
#' @description
#'  Performs Siegel and Castellan's all-pairs comparisons tests
#' of Friedman-type ranked data.
#' @details
#' For all-pairs comparisons in a two factorial unreplicated
#' complete block design
#' with non-normally distributed residuals, Siegel and Castellan's test can be
#' performed on Friedman-type ranked data.
#'
#' A total of \eqn{m = k ( k -1 )/2} hypotheses can be tested.
#' The null hypothesis, H\eqn{_{ij}: \theta_i = \theta_j}, is tested
#' in the two-tailed case against the alternative,
#' A\eqn{_{ij}: \theta_i \ne \theta_j, ~~ i \ne j}.
#'
#' The \eqn{p}-values are computed from the standard normal distribution.
#' Any method as implemented in \code{\link{p.adjust}} can be used for
#' \code{p}-value adjustment.
#' 
#' @references
#' S. Siegel, N. J. Castellan Jr. (1988), \emph{Nonparametric
#'    Statistics for the Behavioral Sciences}. 2nd ed.
#'  New York: McGraw-Hill.
#' 
#' @keywords htest nonparametric
#' @concept Friedman
#' @concept Rank
#' @concept AllPairs
#' @seealso
#' \code{\link{friedmanTest}}, \code{\link{friedman.test}},
#' \code{\link{frdAllPairsExactTest}}, \code{\link{frdAllPairsConoverTest}},
#' \code{\link{frdAllPairsNemenyiTest}}, \code{\link{frdAllPairsMillerTest}}
#' @template class-PMCMR
#' @example examples/frdAllPairs.R
#' @export
frdAllPairsSiegelTest <-
    function(y, ...) UseMethod("frdAllPairsSiegelTest")

#' @rdname frdAllPairsSiegelTest
#' @aliases frdAllPairsSiegelTest.default
#' @method frdAllPairsSiegelTest default
#' @template two-way-parms
#' @param p.adjust.method method for adjusting p values
#'    (see \code{\link{p.adjust}}).
#' @importFrom stats pnorm
#' @export
frdAllPairsSiegelTest.default <-
    function(y, groups, blocks, p.adjust.method = p.adjust.methods, ...)
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
       # GRPNAMES <- as.character(levels(groups))
        k <- nlevels(groups)
        n <- nlevels(blocks)
        GRPNAMES <- as.character(groups[1:k])
    }
        
    p.adjust.method = match.arg(p.adjust.method)
    mat <- matrix(y, nrow = n, ncol = k, byrow = TRUE)
    r <- t(apply(mat, 1L, rank))
    R.mnsum <- colMeans(r)
        
    compare.stats <- function(i,j) {
        dif <- abs(R.mnsum[i] - R.mnsum[j])
        zval <- dif / sqrt(k * (k + 1) / (6 * n))
        return(zval)
    }

    ## z value
    METHOD <- c("Siegel-Castellan all-pairs test for a two-way",
               " balanced complete block design")
    PSTAT <- pairwise.table(compare.stats,
                            levels(groups), p.adjust.method="none")
    pstat <- as.numeric(PSTAT)
    pval <- 2 * pnorm(abs(pstat), lower.tail = FALSE)
    pval[pval > 1] <- 1.0   
    padj <- p.adjust(pval, method = p.adjust.method)
    PVAL <- matrix(padj, ncol = (k - 1), nrow = (k -1 ))
    DIST <- "z"
        		
    colnames(PSTAT) <- GRPNAMES[1:(k-1)]
    rownames(PSTAT) <- GRPNAMES[2:k]
    colnames(PVAL) <- GRPNAMES[1:(k-1)]
    rownames(PVAL) <- GRPNAMES[2:k]
    MODEL <- data.frame(x = y, g = groups, b = blocks)
    ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL,
                statistic = PSTAT, p.adjust.method = p.adjust.method,
                dist = DIST, model = MODEL)
    class(ans) <- "PMCMR"
    ans
}
