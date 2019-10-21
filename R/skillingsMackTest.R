## skillingsMackTest.R
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
##
## PBIBD or PRCB
#' @title Skillings-Mack Test
#' @description
#' Performs Skillings-Mack rank sum test for partially balanced
#' incomplete block designs or partially balanced random block designs.
#' The null hypothesis
#' H\eqn{_0: \theta_i = \theta_j~~(i \ne j)} is tested against the
#' alternative H\eqn{_{\mathrm{A}}: \theta_i \ne \theta_j}, with at least
#' one inequality beeing strict.
#' @name skillingsMackTest
#' @template class-htest
#' @details
#' The function has implemented the test of Skillings and Mack (1981).
#' The test statistic is assymptotically chi-squared distributed with
#' df = k - 1 degrees of freedom.
#'
#' @references
#' Skillings, J. H., Mack, G.A. (1981) On the use of a Friedman-type
#' statistic in balanced and unbalanced block designs,
#' \emph{Technometrics} \bold{23}, 171--177.
#'
#' @note
#' The input vector/matrix \code{'y'} must contain \code{NA}.
#'
#' @examples
#' ## Example from Hollander and Wolfe 1999,
#' ## originally appeared in Brady 1969.
#' x <- cbind(c(3,1,5,2,0,0,0,0),
#'            c(5,3,4,NA,2,2,3,2),
#'            c(15,18,21,6,17,10,8,13))
#' colnames(x) <- c("R", "A", "B")
#' rownames(x) <- 1:8
#' skillingsMackTest(x)
#'
#' ## Compare with Friedman Test for CRB
#' ## Sachs, 1997, p. 675
#' ## Six persons (block) received six different diuretics
#' ## (A to F, treatment).
#' ## The responses are the Na-concentration (mval)
#' ## in the urine measured 2 hours after each treatment.
#'  y <- matrix(c(
#' 3.88, 5.64, 5.76, 4.25, 5.91, 4.33, 30.58, 30.14, 16.92,
#' 23.19, 26.74, 10.91, 25.24, 33.52, 25.45, 18.85, 20.45,
#' 26.67, 4.44, 7.94, 4.04, 4.4, 4.23, 4.36, 29.41, 30.72,
#' 32.92, 28.23, 23.35, 12, 38.87, 33.12, 39.15, 28.06, 38.23,
#' 26.65),nrow=6, ncol=6,
#' dimnames=list(1:6, LETTERS[1:6]))
#' print(y)
#' friedmanTest(y)
#' skillingsMackTest(y)
#' @keywords htest nonparametric
#' @seealso
#' \code{\link{friedmanTest}}, \code{\link{durbinTest}}
#' @export
skillingsMackTest <- function(y, ...) UseMethod("skillingsMackTest")


#' @rdname skillingsMackTest
#' @aliases skillingsMackTest.default
#' @method skillingsMackTest default
#' @template two-way-parms
#' @importFrom MASS ginv
#' @importFrom stats pchisq
#' @export
skillingsMackTest.default <- function(y, groups, blocks, ...)
{

    ##
    if ((is.matrix(y)) | (is.data.frame(y))) {
        groups <- factor(c(col(y)))
        blocks <- factor(c(row(y)))
        DNAME <- paste(deparse(substitute(y)))
        GRPNAMES <- colnames(y)
        y <- as.vector(y)
    }
    else {
      ##  if (any(is.na(groups)) || any(is.na(blocks)))
      ##      stop("NA's are not allowed in groups or blocks")
        if (any(diff(c(length(y), length(groups), length(blocks)))))
            stop("y, groups and blocks must have the same length")
      ##  if (any(table(groups, blocks) != 1))
      ##      stop("Not an unreplicated complete block design")

        DNAME <- paste(deparse(substitute(y)), ",",
                       deparse(substitute(groups)), "and",
                       deparse(substitute(blocks)))
        groups <- factor(groups)
        blocks <- factor(blocks)
        GRPNAMES <- as.character(levels(groups))
    }

    ny <- length(y)
    ng <- length(groups)
    nb <- length(blocks)
    ok <- all.equal(ny, ng, nb)
    if(!ok){
        stop("vectors must all have the same length.")
    }

    ## Nr of treatments in the j-th block
    ok <- complete.cases(y, blocks)
    k <- tapply(y[ok], blocks[ok], length)

    nb <- nlevels(blocks)
    ng <- nlevels(groups)

    ## Friedman type ranking
    y <- y[order(groups, blocks)]
    Y <- matrix(y, nrow = nb, ncol = ng, byrow = FALSE)
    R <- Y
    for (j in 1:length(R[, 1])) R[j, ] <- rank(Y[j, ], na.last="keep")

    ## insert values for NA
    for (i in 1:ng){
        for (j in 1:nb){
            if (is.na(R[j, i])){
                R[j, i] <- (k[j] + 1) / 2
            }
        }
    }

    ## adjusted treatment sums
    A <- sapply(1:ng, function(i) {
        sum(sqrt(12 / (k + 1)) *
            (R[,i] - (k + 1) / 2))
        })


    ## create cov-matrix
    lambda <- matrix(NA, ng, ng)
    for (t in (1: (ng-1))){
        for (q in ((t+1):ng)){
            ok <- complete.cases(Y[,q], Y[,t])
            lambda[q,t] <- length(ok[ok == TRUE])
            lambda[t,q] <- lambda[q,t]
        }
    }
    tmp <- rowSums(lambda, na.rm=TRUE)
    COV <- -lambda
    for (i in (1: ng)){
        COV[i,i] <- tmp[i]
    }

    ## Statistic , requires ginv of package MASS
    T <- t(A) %*% ginv(COV) %*% A

    ## assymptotic chi-square
    df <- ng - 1
    PVAL <- pchisq(T, df = df, lower.tail=FALSE)

    METHOD <- paste("Skillings-Mack test")

    ans <- list(p.value = PVAL,
                statistic = c("Skillings-Mack chi-squared" = T),
                parameter = c("df" = df),
                method = METHOD,
                data.name = DNAME
                )
    class(ans) <- "htest"
    ans
}
