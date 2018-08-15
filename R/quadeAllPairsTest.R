##  quadeAllPairsTest.R
##
##  Copyright (C) 2015-2017 Thorsten Pohlert
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

#' @name quadeAllPairsTest
#' @title All-Pairs Comparisons for
#' Unreplicated Blocked Data (Quade's All-Pairs Test)
#' 
#' @description
#' Performs Quade multiple-comparison test for unreplicated
#' blocked data.
#'
#' @details
#' For all-pairs comparisons of unreplicated blocked data
#' Quade's test can be applied.
#' A total of \eqn{m = k(k-1)/2}
#' hypotheses can be tested. The null hypothesis
#' H\eqn{_{ij}: \theta_i = \theta_j} is tested in the two-tailed test
#' against the alternative
#' A\eqn{_{ij}: \theta_i \ne \theta_j, ~~ i \ne j}.
#'
#' The function has included two methods for approximate p-value estimation:
#' \describe{
#' \item{TDist}{p-values are computed from the t distribution}
#' \item{Normal}{p-values are computed from the standard normal distribution}
#' }
#'
#' If no p-value adjustment is performed (\code{p.adjust.method = "none"}),
#' than a simple protected test is recommended, i.e.
#' all-pairs comparisons should only be applied after a significant
#' \code{\link{quade.test}}. However, any method as implemented in
#' \code{\link{p.adjust.methods}} can be selected by the user.
#'
#' @references
#' W. J. Conover (1999), \emph{Practical nonparametric Statistics},
#' 3rd. Edition, Wiley.
#'
#' N. A. Heckert and J. J. Filliben (2003). NIST Handbook 148:
#' Dataplot Reference Manual, Volume 2: Let Subcommands and Library Functions.
#' National Institute of Standards and Technology Handbook Series, June 2003.
#'
#' D. Quade (1979), Using weighted rankings in the analysis of complete
#' blocks with additive block effects. \emph{Journal of the American
#' Statistical Association}, 74, 680-683.
#'
#' @template class-PMCMR
#'
#' @examples
#' ## Sachs, 1997, p. 675
#' ## Six persons (block) received six different diuretics
#' ## (A to F, treatment).
#' ## The responses are the Na-concentration (mval)
#' ## in the urine measured 2 hours after each treatment.
#' ##
#' y <- matrix(c(
#' 3.88, 5.64, 5.76, 4.25, 5.91, 4.33, 30.58, 30.14, 16.92,
#' 23.19, 26.74, 10.91, 25.24, 33.52, 25.45, 18.85, 20.45,
#' 26.67, 4.44, 7.94, 4.04, 4.4, 4.23, 4.36, 29.41, 30.72,
#' 32.92, 28.23, 23.35, 12, 38.87, 33.12, 39.15, 28.06, 38.23,
#' 26.65),nrow=6, ncol=6,
#' dimnames=list(1:6, LETTERS[1:6]))
#' print(y)
#'
#' ## Global test
#' quade.test(y)
#'
#' ## All-pairs comparisons
#' quadeAllPairsTest(y, dist="TDist", p.adjust.method="holm")
#' 
#' @keywords htest nonparametric
#' @concept AllPairsComparison
#' @concept TwoWayRankAnova
#' @seealso
#' \code{\link{quade.test}}, \code{\link{friedmanTest}}
#' @export
quadeAllPairsTest <- function(y, ...) UseMethod("quadeAllPairsTest")

#' @rdname quadeAllPairsTest
#' @aliases quadeAllPairsTest.default
#' @method quadeAllPairsTest default
#' @template two-way-parms
#' @param dist the test distribution. Defaults to \code{"TDist"}.
#' @param p.adjust.method method for adjusting p values
#' (see \code{\link{p.adjust}}).
#' @importFrom stats pt
#' @importFrom stats pnorm
#' @importFrom stats complete.cases
#' @importFrom stats pairwise.table
#' @importFrom stats p.adjust
#' @importFrom stats p.adjust.methods
#' @export
quadeAllPairsTest.default <-
function(y, groups, blocks, dist=c("TDist", "Normal"),
         p.adjust.method = p.adjust.methods,  ...)
{
    if ((is.matrix(y)) | (is.data.frame(y))) {
        ##

        DNAME <- paste(deparse(substitute(y)))
        GRPNAMES <- colnames(y)
        k <- length(GRPNAMES)
        BLOCKNAMES <- rownames(y)
        b <- length(BLOCKNAMES)
        groups <- factor(rep(GRPNAMES, times = b))
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
                           deparse(substitute(blocks)))
            groups <- factor(groups)
            blocks <- factor(blocks)
            k <- nlevels(groups)
            b <- nlevels(blocks)
            GRPNAMES <- levels(groups)
        }
    mat <- matrix(y, nrow = b, ncol = k, byrow = TRUE)        
    dist <- match.arg(dist)
    p.adjust.method <- match.arg(p.adjust.method)
    r <- t(apply(mat, 1L, rank))
    q <- rank(apply(mat, 1, function(u) max(u) - min(u)))
    s <- q * (r - (k+1)/2)
    w <- q * r
    ## s is a matrix of ranks within blocks (minus the average rank)
    ## multiplied by the ranked ranges of the blocks
    A <- sum(s^2)
    B <- sum(colSums(s)^2) / b
    S <- colSums(s)
    W <- colSums(w)
    if (dist == "TDist") {
        METHOD <- "Quade's test with TDist approximation"
        DIST <- "t"
        denom <- sqrt((2 * b * (A - B))/
                      ((b-1) * (k-1)))
        compare.stats <- function(i,j) {
            dif <- abs(S[i] - S[j]) 
            tval <- dif / denom
            return(tval)
        }
        PSTAT <- pairwise.table(compare.stats,levels(groups),
                                p.adjust.method="none")
        compare.levels <- function(i,j) {
            dif <- abs(S[i] - S[j]) 
            tval <- dif / denom
            pval <- pval <- 2 * pt(abs(tval),
                                   df=(b-1)*(k-1),
                                   lower.tail=FALSE)
            return(pval)
        }
        PVAL <- pairwise.table(compare.levels,levels(groups),
                               p.adjust.method=p.adjust.method)
    } else {
         METHOD <- paste("Quade's test",
                         "with standard-normal approximation",
                         sep="\t")
         DIST <- "z"
         n <- b * k
         denom <- sqrt((k * (k + 1) * (2 * n + 1) * (k-1))/
                           (18 * n * (n + 1)))
         nn <- length(w[,1])
         ff <- 1 / (nn * (nn + 1)/2)
         compare.stats <- function(i,j) {
            dif <- abs(W[i] * ff - W[j] * ff) 
            zval <- dif / denom
            return(zval)
        }
        PSTAT <- pairwise.table(compare.stats,levels(groups),
                                p.adjust.method="none")
        compare.levels <- function(i,j) {
            dif <- abs(W[i] * ff - W[j] * ff) 
            zval <- dif / denom
            pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
            return(pval)
        }
        PVAL <- pairwise.table(compare.levels,levels(groups),
                               p.adjust.method=p.adjust.method)
    }
    colnames(PSTAT) <- GRPNAMES[1:(k-1)]
    rownames(PSTAT) <- GRPNAMES[2:k]
    colnames(PVAL) <- GRPNAMES[1:(k-1)]
    rownames(PVAL) <- GRPNAMES[2:k]
    MODEL <- data.frame(x= y, groups, blocks)
    ans <- list(statistic = PSTAT,
                   p.value = PVAL,
                   method = METHOD,
                   p.adjust.method = p.adjust.method,
                   data.name = DNAME,
                model = MODEL,
                dist = DIST)
    class(ans) <- "PMCMR"
    ans
}
