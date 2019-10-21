## frdAllPairsExactTest.R
## Part of the R package: PMCMR
##
## Copyright (C) 2017-2019 Thorsten Pohlert
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
##  Literature:
##  Eisinga, Heskes, Pelzer, Te Grotenhuis, 2017,
##  Exact p-values for pairwise comparison of Friedman rank sums, with application to
##  comparing classifiers, BMC Bioinformatics, January 2 2017

#' @name frdAllPairsExactTest
#' @title Exact All-Pairs Comparisons Test for Unreplicated Blocked Data
#' @description
#'  Performs exact all-pairs comparisons tests of Friedman-type ranked data
#' according to Eisinga et al. (2017).
#' @details
#' For all-pairs comparisons in a two factorial unreplicated
#' complete block design
#' with non-normally distributed residuals, an exact test can be
#' performed on Friedman-type ranked data.
#'
#' A total of \eqn{m = k ( k -1 )/2} hypotheses can be tested.
#' The null hypothesis, H\eqn{_{ij}: \theta_i = \theta_j}, is tested
#' in the two-tailed case against the alternative,
#' A\eqn{_{ij}: \theta_i \ne \theta_j, ~~ i \ne j}.
#'
#' The exact \eqn{p}-values
#' are computed using the code of \code{"pexactfrsd.R"}
#' that was a supplement to the publication of Eisinga et al. (2017).
#' Additionally, any of the \eqn{p}-adjustment methods
#' as included in \code{\link{p.adjust}} can be selected, for \eqn{p}-value
#' adjustment.
#'
#' @references
#' Eisinga, R., Heskes, T., Pelzer, B., Te Grotenhuis, M. (2017)
#'  Exact p-values for Pairwise Comparison of Friedman Rank Sums,
#'  with Application to Comparing Classifiers, \emph{BMC Bioinformatics}, 18:68.
#'
#' @keywords htest nonparametric
#' @concept Friedman
#' @concept Rank
#' @concept AllPairs
#'
#' @section Source:
#' The function \code{frdAllPairsExactTest} uses the code
#' of the file \code{pexactfrsd.R} that was a supplement to:\cr
#'
#' R. Eisinga, T. Heskes, B. Pelzer, M. Te Grotenhuis (2017),
#'  Exact p-values for Pairwise Comparison of Friedman Rank Sums,
#'  with Application to Comparing Classifiers, \emph{BMC Bioinformatics}, 18:68.
#' @template class-PMCMR
#' @seealso
#' \code{\link{friedmanTest}}, \code{\link{friedman.test}},
#' \code{\link{frdAllPairsConoverTest}}, \code{\link{frdAllPairsMillerTest}},
#' \code{\link{frdAllPairsNemenyiTest}}, \code{\link{frdAllPairsSiegelTest}}
#' @example examples/frdAllPairs.R
#' @export
frdAllPairsExactTest <-
    function(y, ...) UseMethod("frdAllPairsExactTest")

#' @rdname frdAllPairsExactTest
#' @aliases frdAllPairsExactTest.default
#' @method frdAllPairsExactTest default
#' @template two-way-parms
#' @param p.adjust.method method for adjusting p values
#'    (see \code{\link{p.adjust}}).
#' @importFrom Rmpfr mpfr
#' @import gmp
#' @importFrom stats pairwise.table
#' @importFrom stats p.adjust.methods
#' @importFrom stats p.adjust
#' @export
frdAllPairsExactTest.default <-
    function(y, groups, blocks, p.adjust.method =  p.adjust.methods, ...)
{
    p.adjust.method <- match.arg(p.adjust.method)
    ## 2019-10-16
    ## novel external function
    ans <- frdRanks(y, groups, blocks)
    r <- ans$r
    n <- nrow(r)
    k <- ncol(r)
    GRPNAMES <- colnames(r)

    METHOD <- c("Eisinga, Heskes, Pelzer & Te Grotenhuis all-pairs test",
                " with exact p-values for a two-way",
                " balanced complete block design")
    Ri <- colSums(r)
    compareRi <- function(i, j){
        d <- Ri[i] - Ri[j]
        d
    }
    PSTAT <- pairwise.table(compareRi, GRPNAMES,
                            p.adjust.method="none")
    pstat <- as.numeric(PSTAT)
    pstat <- pstat[!is.na(pstat)]
    pval <- sapply(pstat, function(pstat) {
        pexactfrsd(d = abs(pstat),k = k,n = n, option="pvalue")
    })
    padj <- p.adjust(pval, method = p.adjust.method, n = (k * (k - 1) / 2))
    PVAL <- matrix(NA, ncol = (k - 1), nrow = (k -1 ))
    PVAL[!is.na(PSTAT)] <- padj
    DIST <- "D"
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
                model = ans$inDF)
    class(ans) <- "PMCMR"
    ans
}
