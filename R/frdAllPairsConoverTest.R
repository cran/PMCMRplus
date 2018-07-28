## posthoc.friedman.conover.test.R
## Part of the R package: PMCMR
##
## Copyright (C) 2015-2018 Thorsten Pohlert
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

#' @name frdAllPairsConoverTest
#' @title Conover's All-Pairs Comparisons Test for Unreplicated Blocked Data
#' @description
#'  Performs Conover's all-pairs comparisons tests of Friedman-type ranked data.
#' @details
#' For all-pairs comparisons in a two factorial unreplicated
#' complete block design
#' with non-normally distributed residuals, Conover's test can be
#' performed on Friedman-type ranked data.
#'
#' A total of \eqn{m = k ( k -1 )/2} hypotheses can be tested.
#' The null hypothesis, H\eqn{_{ij}: \theta_i = \theta_j}, is tested
#' in the two-tailed case against the alternative,
#' A\eqn{_{ij}: \theta_i \ne \theta_j, ~~ i \ne j}.
#'
#' If \code{p.adjust.method == "single-step"} the p-values are computed
#' from the studentized range distribution. Otherwise,
#' the p-values are computed from the t-distribution using
#' any of the p-adjustment methods as included in \code{\link{p.adjust}}.
#'
#' @inherit durbinAllPairsTest references
# @references
# W. J. Conover and R. L. Iman (1979), \emph{On multiple-comparisons
#  procedures}, Tech. Rep. LA-7677-MS, Los Alamos Scientific Laboratory.
#
# W. J. Conover (1999), \emph{Practical nonparametric Statistics},
#  3rd. Edition, Wiley.
#
#' @keywords htest nonparametric
#' @concept Friedman
#' @concept Rank
#' @concept AllPairs
#' @seealso
#' \code{\link{friedmanTest}}, \code{\link{friedman.test}},
#' \code{\link{frdAllPairsExactTest}}, \code{\link{frdAllPairsMillerTest}},
#' \code{\link{frdAllPairsNemenyiTest}}, \code{\link{frdAllPairsSiegelTest}}
#' @template class-PMCMR
#' @example examples/frdAllPairs.R
#' @export
frdAllPairsConoverTest <- function(y, ...)
    UseMethod("frdAllPairsConoverTest")

#' @rdname frdAllPairsConoverTest
#' @aliases frdAllPairsConoverTest.default
#' @method frdAllPairsConoverTest default
#' @template two-way-parms
#' @param p.adjust.method method for adjusting p values
#'    (see \code{\link{p.adjust}}).
#' @importFrom stats pt
#' @importFrom stats ptukey
#' @export
frdAllPairsConoverTest.default <-
    function(y, groups, blocks,
             p.adjust.method = c("single-step", p.adjust.methods), ...)
{
    if ((is.matrix(y)) | (is.data.frame(y))) {
        #groups <- factor(c(col(y)))
        #blocks <- factor(c(row(y)))
        DNAME <- paste(deparse(substitute(y)))
        GRPNAMES <- colnames(y)
        k <- length(GRPNAMES)
        BLOCKNAMES <- rownames(y)
        n <- length(BLOCKNAMES)
        groups <- factor(rep(GRPNAMES, times = n))
        blocks <- factor(rep(BLOCKNAMES, each = k))
        #y <- y[order(groups, blocks)]
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
                       deparse(substitute(blocks)))
        groups <- factor(groups)
        blocks <- factor(blocks)
        k <- nlevels(groups)
        n <- nlevels(blocks)
        GRPNAMES <- as.character(groups[1:k])
    }
    p.adjust.method <- match.arg(p.adjust.method)
    #n <- length(levels(blocks))
    #k <- length(levels(groups))
    #y <- y[order(groups, blocks)]
    #mat <- matrix(y, nrow = n, ncol = k, byrow = FALSE)
    #for (i in 1:length(mat[, 1])) mat[i, ] <- rank(mat[i, ])
    mat <- matrix(y, nrow = n, ncol = k, byrow = TRUE)
    r <- t(apply(mat, 1L, rank))
    R.sum <- colSums(r)
    METHOD <- c("Conover's all-pairs test for a two-way",
                " balanced complete block design")

    ## re coded from Conover and Imam 1978
    m <- 1 # replicates set to 1
    S2 <- m / ( m * k -1 ) * (sum(r^2) - m * k * n *
                              ( m * k + 1)^2 / 4)
    T2 <- 1 / S2 * (sum(R.sum) - n * m * ((m * k + 1) / 2)^2)
    A <-S2 * (2 * n * (m * k - 1)) / ( m * n * k - k - n + 1)
    B <- 1 - T2 / (n * (m * k- 1))

    if (p.adjust.method != "single-step"){

        compare.stats <- function(i, j){
            diff <- R.sum[i] - R.sum[j]
            tval <- diff / (sqrt(A) * sqrt(B))
            tval
        }
        PSTAT <- pairwise.table(compare.stats,levels(groups),
                                p.adjust.method="none")

        compare.levels <- function(i,j) {
            dif <- abs(R.sum[i] - R.sum[j])
            tval <- dif / (sqrt(A) * sqrt(B))
            pval <- 2 * pt(q=abs(tval), df=(m * n * k - k - n + 1),
                           lower.tail=FALSE)
            return(pval)
        }
        PVAL <- pairwise.table(compare.levels, levels(groups),
                               p.adjust.method=p.adjust.method)
        DIST <- "t"
        PARMS <- m * n * k - k - n + 1
        names(PARMS) <- "df"

    } else {
        ## use Tukey distribution for multiple comparisons
        compare.stats <- function(i, j){
            diff <- R.sum[i] - R.sum[j]
            qval <- sqrt(2) * diff / (sqrt(A) * sqrt(B))
            qval
        }
        PSTAT <- pairwise.table(compare.stats,levels(groups),
                                p.adjust.method="none")

        compare.levels <- function(i,j) {
            dif <- abs(R.sum[i] - R.sum[j])
            qval <- sqrt(2) * dif / (sqrt(A) * sqrt(B))
            pval <-  ptukey(q=qval, nmeans = k,
                            df= Inf,
                            lower.tail=FALSE)
            return(pval)
        }

        PVAL <- pairwise.table(compare.levels, levels(groups),
                               p.adjust.method= "none")
        DIST <- "q"
        PARMS <- c(k, Inf)
        names(PARMS) <- c("nmeans", "df")
    }

    colnames(PSTAT) <- GRPNAMES[1:(k-1)]
    rownames(PSTAT) <- GRPNAMES[2:k]
    colnames(PVAL) <- GRPNAMES[1:(k-1)]
    rownames(PVAL) <- GRPNAMES[2:k]
    MODEL <- data.frame(x = y, g = groups, b = blocks)
    ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL,
                statistic = PSTAT, p.adjust.method =p.adjust.method,
                model = MODEL, dist = DIST, parameter = PARMS)
    class(ans) <- "PMCMR"
    ans
}
