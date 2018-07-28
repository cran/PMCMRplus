##  durbinAllPairsTest.R
##  Part of the R package PMCMRplus
##
##  Copyright (C) 2015-2018 Thorsten Pohlert
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

#' @name durbinAllPairsTest
#'
#' @title All-Pairs Comparisons Test for Balanced Incomplete Block Designs
#'
#' @description
#' Performs Conover-Iman all-pairs comparison test for a balanced incomplete
#' block design (BIBD).
#'
#' @details
#' For all-pairs comparisons in a balanced incomplete block design
#' the proposed test of Conover and Imam can be applied.
#' A total of \eqn{m = k(k-1)/2}
#' hypotheses can be tested. The null hypothesis
#' H\eqn{_{ij}: \theta_i = \theta_j} is tested in the two-tailed test
#' against the alternative
#' A\eqn{_{ij}: \theta_i \ne \theta_j, ~~ i \ne j}.
#'
#' The p-values are computed from the t distribution. If no p-value adjustment
#' is performed (\code{p.adjust.method = "none"}),
#' than a simple protected test is recommended, i.e.
#' the all-pairs comparisons should only be applied after a significant
#' \code{\link{durbinTest}}. However, any method as implemented in
#' \code{\link{p.adjust.methods}} can be selected by the user.
#'
#' @references
#' Conover, W. J., Iman, R. L. (1979) \emph{On multiple-comparisons
#'  procedures}, Tech. Rep. LA-7677-MS, Los Alamos Scientific Laboratory.
#'
#' Conover, W. J. (1999) \emph{Practical nonparametric Statistics},
#' 3rd. Edition, Wiley.
#'
#' @examples
#' ## Example for an incomplete block design:
#' ## Data from Conover (1999, p. 391).
#' y <- matrix(c(2,NA,NA,NA,3, NA,  3,  3,  3, NA, NA, NA,  3, NA, NA,
#'   1,  2, NA, NA, NA,  1,  1, NA,  1,  1,
#' NA, NA, NA, NA,  2, NA,  2,  1, NA, NA, NA, NA,
#'  3, NA,  2,  1, NA, NA, NA, NA,  3, NA,  2,  2),
#' ncol=7, nrow=7, byrow=FALSE, dimnames=list(1:7, LETTERS[1:7]))
#' durbinAllPairsTest(y)
#' @keywords htest nonparametric
#' @template class-PMCMR
#' @seealso
#' \code{\link{durbinTest}}
#' @export
durbinAllPairsTest <- function(y, ...) UseMethod("durbinAllPairsTest")

#' @rdname durbinAllPairsTest
#' @aliases durbinAllPairsTest.default
#' @method durbinAllPairsTest default
#' @template two-way-parms
#' @param p.adjust.method method for adjusting p values
#' (see \code{\link{p.adjust}})
#' @importFrom stats pt
#' @importFrom stats complete.cases
#' @importFrom stats pairwise.table
#' @importFrom stats p.adjust
#' @importFrom stats p.adjust.methods
#' @export
durbinAllPairsTest.default <- function(y, groups, blocks,
                                       p.adjust.method = p.adjust.methods, ...)
{
    DNAME <- deparse(substitute(y))
    if (is.matrix(y)) {
        GRPNAME <- colnames(y)
        g <- factor(c(col(y)))
        b <- factor(c(row(y)))
        yy <- as.vector(y)
        datf1 <- data.frame(yy, b, g)
        datf2 <- datf1[!is.na(datf1[,1]),]
        blocks <- factor(datf2[,2])
        groups <- factor(datf2[,3])
        y <- datf2[,1]
        ## Clean up
        rm(b, g, yy, datf1, datf2)
    }
    else {
        if (anyNA(groups) || anyNA(blocks))
            stop("NA's are not allowed in 'groups' or 'blocks'")
        if (any(diff(c(length(y), length(groups), length(blocks))) != 0L))
            stop("'y', 'groups' and 'blocks' must have the same length")
        DNAME <- paste(DNAME, ", ", deparse(substitute(groups)),
                       " and ", deparse(substitute(blocks)), sep = "")
        groups <- factor(groups)
        blocks <- factor(blocks)
        GRPNAME <- levels(groups)
    }

    ## Need to ensure consistent order.
    o <- order(blocks, groups)
    y <- y[o]
    groups <- groups[o]
    blocks <- blocks[o]

    p.adjust.method = match.arg(p.adjust.method)
    t <- nlevels(groups)
    b <- nlevels(blocks)
    r <- unique(table(groups))
    k <- unique(table(blocks))
    rij <- unlist(tapply(y, blocks, rank))
    Rj <- tapply(rij, groups, sum)
    ## Taken from NIST
    A <- sum(rij^2)
    C <- (b * k * (k + 1)^2) / 4
    D <- sum(Rj^2) - r * C
    T1 <- (t - 1) / (A - C) * D

    denom <- sqrt(((A - C) * 2 * r) / (b * k - b - t + 1) *
                      (1 - T1 / (b * (k -1))))
    df <- b * k - b - t + 1
   # Pairwise comparisons
    compare.stats <- function(i,j) {
        dif <- abs(Rj[i] - Rj[j])
        tval <- dif / denom
        return(tval)
    }
    PSTAT <- pairwise.table(compare.stats,levels(groups),
                            p.adjust.method="none" )

    compare.levels <- function(i,j) {
        dif <- abs(Rj[i] - Rj[j])
        tval <- dif / denom
        pval <- 2 * pt(q=abs(tval), df=df, lower.tail=FALSE)
        return(pval)
    }
    PVAL <- pairwise.table(compare.levels,levels(groups),
                           p.adjust.method=p.adjust.method)

    METHOD <- c("Durbin's all-pairs test for a two-way",
                " balanced incomplete block design")
    colnames(PSTAT) <- GRPNAME[1:(t-1)]
    rownames(PSTAT) <- GRPNAME[2:t]
    colnames(PVAL) <- colnames(PSTAT)
    rownames(PVAL) <- rownames(PSTAT)
    MODEL <- data.frame(x = y, groups, blocks)
    DIST <- "t"
    PARAMETER <- df
    names(PARAMETER) <- "df"
    ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL,
                statistic = PSTAT, p.adjust.method = p.adjust.method,
                model = MODEL, dist = DIST, parameter = PARAMETER)
    class(ans) <- "PMCMR"
    ans
}
