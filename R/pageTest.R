## pageTest.R
## Part of the R package: PMCMRplus
##
## Copyright (C) 2017, 2018 Thorsten Pohlert
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

#' @name pageTest
#' @title Page Rank Sum Test
#'
#' @description
#' Performs Page's ordered aligned rank sum test.
#'
#' @template class-htest
#'
#' @references
#' Page, E. B. (1963) Ordered hypotheses for multiple treatments: A
#' significance test for linear ranks, \emph{Journal of the
#' American Statistical Association} \bold{58}, 216--230.
#'
#' Sachs, L. (1997) \emph{Angewandte Statistik}. Berlin: Springer.
#'
#' @examples
#' ## Sachs (1997), pp. 671 ff.
#' ## 9 reviewers (blocks)
#' ## assigned ranks to 4 objects (groups).
#' data(reviewers)
#' ## See Sachs (1997) p. 677
#' pageTest(reviewers, alternative = "greater")
#'
#' @keywords htest nonparametric
#' @concept OrderedAlternative
#' @concept TwoWayRankAnova
#'
#' @seealso
#' \code{\link{friedmanTest}}
#'
#' @export
pageTest <- function(y, ...) UseMethod("pageTest")

#' @rdname pageTest
#' @aliases pageTest.default
#' @method pageTest default
#' @template two-way-parms
#' @param alternative the alternative hypothesis.
#' Defaults to \code{two.sided}.
#' @importFrom stats pnorm
#' @export
pageTest.default <-
    function(y, groups, blocks, alternative = c("two.sided", "greater",
                                                "less"), ...){
        if ((is.matrix(y)) | (is.data.frame(y))) {
            groups <- factor(c(col(y)))
            blocks <- factor(c(row(y)))
            DNAME <- paste(deparse(substitute(y)))
            GRPNAMES <- colnames(y)
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
            GRPNAMES <- as.character(levels(groups))
        }

        alternative <- match.arg(alternative)
    	n <- length(levels(blocks))
        k <- length(levels(groups))
    	y <- y[order(groups, blocks)]
    	mat <- matrix(y, nrow = n, ncol = k, byrow = FALSE)
        for (i in 1:length(mat[, 1])) mat[i, ] <- rank(mat[i, ])

        # Sachs p. 676
        R.sum <- colSums(mat)
        L <- sum(R.sum * 1:k)
        eL <- n * k * (k + 1)^2 / 4
        varL <- n * k^2 * (k + 1) * (k^2 - 1) / 144

        METHOD <- paste("Page's ordered aligned rank sum test")
        STAT <- (L - eL - 1/2) / sqrt(varL)

        if (alternative == "two.sided"){
            PVAL <- 2 * min(pnorm(abs(STAT), lower.tail = FALSE), 0.5)
        } else if (alternative == "greater") {
            PVAL <- pnorm(STAT, lower.tail = FALSE)
        } else{
            PVAL <- pnorm(STAT)
        }
        names(STAT) <- "z"
        names(L) <- "L"
        ans <- list(statistic=STAT, p.value = PVAL,
                    method = METHOD, estimate = L,
                    data.name=DNAME, alternative = alternative)
        class(ans) <- "htest"
        ans
    }
