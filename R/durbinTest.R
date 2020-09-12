##  durbin.test.R
##  Part of the R package PMCMRplus
##
##  Copyright (C) 2015-2020 Thorsten Pohlert
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

#' @name durbinTest
#' @title Durbin Test
#'
#' @description
#' Performs Durbin's tests whether k groups
#' (or treatments) in a two-way balanced incomplete block design (BIBD)
#' have identical effects.
#'
#' @details
#' For testing a two factorial layout of a balanced incomplete
#' block design whether the \eqn{k} groups have identical effects,
#' the Durbin test can be performed. The null hypothesis,
#' H\eqn{_0: \theta_i = \theta_j ~ (1 \le i < j \le k)},
#' is tested against the alternative that at least
#' one \eqn{\theta_i \ne \theta_j}.
#'
#' The p-values are computed from the chi-square distribution.
#'
#' @note
#' The function does not test, whether it is a true BIBD.
#' This function does not test for ties.
#'
#' @references
#' Conover,W. J. (1999) \emph{Practical nonparametric Statistics},
#' 3rd. Edition, Wiley.
#'
#' Heckert, N. A., Filliben, J. J. (2003) \emph{NIST Handbook 148:
#' Dataplot Reference Manual}, Volume 2:
#' Let Subcommands and Library Functions.
#' National Institute of Standards and Technology Handbook Series, June 2003.
#'
#' @keywords htest
#' @keywords nonparametric
#' @concept friedmanranks
#' @examples
#' ## Example for an incomplete block design:
#' ## Data from Conover (1999, p. 391).
#' y <- matrix(c(
#' 2,NA,NA,NA,3, NA,  3,  3,  3, NA, NA, NA,  3, NA, NA,
#'   1,  2, NA, NA, NA,  1,  1, NA,  1,  1,
#' NA, NA, NA, NA,  2, NA,  2,  1, NA, NA, NA, NA,
#'  3, NA,  2,  1, NA, NA, NA, NA,  3, NA,  2,  2
#' ), ncol=7, nrow=7, byrow=FALSE,
#' dimnames=list(1:7, LETTERS[1:7]))
#' durbinTest(y)
#' @template class-htest
#' @export
durbinTest <- function(y, ...) UseMethod("durbinTest")

#' @rdname durbinTest
#' @aliases durbinTest.default
#' @method durbinTest default
#' @template two-way-parms
#' @importFrom stats pchisq
#' @export
durbinTest.default <- function(y, groups, blocks, ...)
{
    DNAME <- deparse(substitute(y))

    if (is.matrix(y)) {
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
    }

    ## Need to ensure consistent order.
    o <- order(blocks, groups)
    y <- y[o]
    groups <- groups[o]
    blocks <- blocks[o]

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

    STATISTIC <- c("Durbin chi-squared" = T1)
    PARAMETER <- c(df = t - 1)
    PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)

    ans <- list(statistic = STATISTIC,
                parameter = PARAMETER,
                p.value = PVAL,
                method = c("Durbin's rank sum test for a two-way",
                           " balanced incomplete block design"),
                data.name = DNAME)
    class(ans) = "htest"
    return(ans)
}
