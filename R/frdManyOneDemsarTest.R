## frdManyOneDemsarTest.R
## Part of the R package: PMCMR
##
## Copyright (C) 2017-2020 Thorsten Pohlert
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

#' @name frdManyOneDemsarTest
#' @title  Demsar's Many-to-One Test
#' for Unreplicated Blocked Data
#'
#' @description
#' Performs Demsar's non-parametric many-to-one comparison test
#' for Friedman-type ranked data.
#'
#' @details
#' For many-to-one comparisons (pairwise comparisons with one control)
#' in a two factorial unreplicated complete block design
#' with non-normally distributed residuals, Demsar's test can be
#' performed on Friedman-type ranked data.
#'
#' Let there be \eqn{k} groups including the control,
#' then the number of treatment levels is \eqn{m = k - 1}.
#' A total of \eqn{m} pairwise comparisons can be performed between
#' the \eqn{i}-th treatment level and the control.
#' H\eqn{_i: \theta_0 = \theta_i} is tested in the two-tailed case against
#' A\eqn{_i: \theta_0 \ne \theta_i, ~~ (1 \le i \le m)}.
#'
#' The \eqn{p}-values are computed from the standard normal distribution.
#' Any of the \eqn{p}-adjustment methods as included in \code{\link{p.adjust}}
#' can be used for the adjustment of \eqn{p}-values.
#'
#' @references
#' Demsar, J. (2006) Statistical comparisons of classifiers over multiple
#'  data sets, \emph{Journal of Machine Learning Research} \bold{7}, 1--30.
#'
#' @concept friedmanranks
#' @keywords htest nonparametric
#' @example examples/frdManyOne.R
#' @seealso
#' \code{\link{friedmanTest}}, \code{\link[stats]{friedman.test}},
#' \code{\link{frdManyOneExactTest}}, \code{\link{frdManyOneNemenyiTest}}.
#'
#' @template class-PMCMR
#' @export
frdManyOneDemsarTest <- function(y, ...) UseMethod("frdManyOneDemsarTest")

#' @rdname frdManyOneDemsarTest
#' @method frdManyOneDemsarTest default
#' @aliases frdManyOneDemsarTest.default
#' @template two-way-parms
#' @param alternative the alternative hypothesis. Defaults to \code{two.sided}.
#' @param p.adjust.method method for adjusting p values
#'    (see \code{\link{p.adjust}}).
#' @importFrom stats pnorm
#' @importFrom stats complete.cases
#' @importFrom stats p.adjust
#' @importFrom stats p.adjust.methods
#' @export
frdManyOneDemsarTest.default <-
    function(y, groups, blocks, alternative = c("two.sided", "greater", "less"),
             p.adjust.method= p.adjust.methods, ...)
{
    ## Check arguments
    p.adjust.method <- match.arg(p.adjust.method)
    alternative <- match.arg(alternative)

    ## 2019-10-16
    ## novel external function
    ans <- frdRanks(y, groups, blocks)
    r <- ans$r
    n <- nrow(r)
    k <- ncol(r)
    GRPNAMES <- colnames(r)

    ## z value
    METHOD <- c("Demsar's many-to-one test for a two-way",
                      " balanced complete block design")
    PSTATv <- rep(NA, k-1)
    R.mnsum <- colMeans(r)
    compare.stats <- function(j) {
        dif <- R.mnsum[j] - R.mnsum[1]
        val <- dif / sqrt(k * (k + 1) / (6 * n))
        return(val)
    }
    for (j in 2:k) {PSTATv[j-1] <- compare.stats(j)}
    ## unadjusted p-values
    if (alternative == "two.sided"){
                PVALv <- 2 * pnorm(abs(PSTATv), lower.tail = FALSE)
    } else if (alternative == "greater"){
        PVALv <- pnorm(PSTATv, lower.tail = FALSE)
    } else {
        PVALv <- pnorm(PSTATv)
            }
    ## adjusted p-values
    PADJv <- p.adjust(PVALv, method = p.adjust.method)

    LNAME <- GRPNAMES[2:k]
    DIST <- "z"

    ## build matrix
    PSTAT <- matrix(data=PSTATv, nrow = (k-1), ncol = 1,
                    dimnames = list(LNAME, GRPNAMES[1]))
    PVAL <- matrix(data=PADJv, nrow = (k-1), ncol = 1,
                   dimnames = list(LNAME, GRPNAMES[1]))
    #MODEL <- data.frame(x = y, g = groups, b = blocks)
    ans <- list(method = METHOD,
                data.name = ans$DNAME,
                p.value = PVAL,
                statistic = PSTAT,
                p.adjust.method = p.adjust.method,
                alternative = alternative,
                dist = DIST,
                model = ans$MODEL)
    class(ans) <- "PMCMR"
    ans
}
