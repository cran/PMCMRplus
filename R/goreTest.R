## goreTest.R
## Part of the R package: PMCMR
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
##
##
#' @title Gore Test
#' @description
#' Performs Gore's test. The null hypothesis
#' H\eqn{_0: \theta_i = \theta_j~~(i \ne j)} is tested against the
#' alternative H\eqn{_{\mathrm{A}}: \theta_i \ne \theta_j}, with at least
#' one inequality beeing strict.
#' @name goreTest
#' @param y a numeric vector of data values.
#' @param groups a vector or factor object giving the group for the
#'          corresponding elements of \code{"y"}.
#' @param blocks a vector or factor object giving the group for the
#'          corresponding elements of \code{"y"}.
#' @template class-htest
#' @details
#' The function has implemented Gore's test for testing
#' main effects in unbalanced CRB designs,
#' i.e. there are one ore more observations per cell.
#' The statistic is assymptotically chi-squared distributed.
#' @references
#' Gore, A. P. (1975) Some nonparametric tests and selection
#' procedures for main effects in two-way layouts.
#' \emph{Ann. Inst. Stat. Math.} \bold{27}, 487--500.
#'
#' @keywords htest nonparametric
#' @seealso
#' \code{\link{friedmanTest}}, \code{\link{skillingsMackTest}},
#' \code{\link{durbinTest}}
#' @examples
#' ## Crop Yield of 3 varieties on two
#' ## soil classes
#' X <-c("130,A,Light
#' 115,A,Light
#' 123,A,Light
#' 142,A,Light
#' 117,A,Heavy
#' 125,A,Heavy
#' 139,A,Heavy
#' 108,B,Light
#' 114,B,Light
#' 124,B,Light
#' 106,B,Light
#' 91,B,Heavy
#' 111,B,Heavy
#' 110,B,Heavy
#' 155,C,Light
#' 146,C,Light
#' 151,C,Light
#' 165,C,Light
#' 97,C,Heavy
#' 108,C,Heavy")
#' con <- textConnection(X)
#' x <- read.table(con, header=FALSE, sep=",")
#' close(con)
#' colnames(x) <- c("Yield", "Variety", "SoilType")
#' goreTest(y = x$Yield, groups = x$Variety, blocks = x$SoilType)
#' @importFrom stats pchisq
#' @export
goreTest <- function(y, groups, blocks)
{
    DNAME <- paste(deparse(substitute(y)), ",",
                   deparse(substitute(groups)), "and",
                   deparse(substitute(blocks)))
    if (any(is.na(groups)) || any(is.na(blocks)))
        stop("NA's are not allowed in groups or blocks")
    if (any(diff(c(length(y), length(groups), length(blocks)))))
        stop("y, groups and blocks must have the same length")
    if (any(is.na(y)))
        stop("NA's are not allowed in y")

    groups <- factor(groups)
    blocks <- factor(blocks)

    ny <- length(y)
    ng <- length(groups)
    nb <- length(blocks)
    ok <- all.equal(ny, ng, nb)
    if(!ok){
        stop("vectors must all have the same length.")
    }

    ## Nr of blocks
    b <- nlevels(blocks)

    ## Nr. of treatment levels
    v <- nlevels(groups)

    ## Total number of observations
    N <- length(y)

    n <- matrix(NA, ncol=b, nrow=v)
    for(i in (1:v)){
        for (j in (1:b)){
            n[i,j] <- length(y[groups == levels(groups)[i]
                               & blocks == levels(blocks)[j]])
        }
    }

    p <- n / N
    q <- 1 / p
    qi <- rowSums(q)
    qstar <- sum(1 / qi)

    ## Create U statistic
    ## see Gore, 1975, p. 488, Eq. 2.1
    phi.fn <- function(x){
        (sign(x) + 1) / 2
    }

    ##
    U <- rep(0, v)
    c <- 0
    for (i in  1:v){
        for (ii in (1:v)){
            if (ii != i){
                for (j in (1:b)){
                    c <- c + 1
                    tmp <- 0
                    nij <- n[i,j]
                    Xijk <- y[groups == levels(groups)[i]
                              & blocks == levels(blocks)[j]]
                    niij <- n[ii,j]
                    Xiijl <-  y[groups == levels(groups)[ii]
                                & blocks == levels(blocks)[j]]
                    for (k in (1:nij)){
                        for (l in (1:niij)){
                            tmp <- tmp + phi.fn(Xijk[k] - Xiijl[l])
                        }
                    }
                    U[i] <- U[i] + tmp / (nij * niij)
                }
            }
        }
    }

    S <- (12 * N / v^2) *
        (sum((U - (v - 1) * b / 2)^2 / qi) -
         sum((U - (v - 1) * b / 2) / qi)^2 / qstar)
    df <- v - 1

    PVAL <- pchisq(S, df = df, lower.tail=FALSE)
    METHOD <- paste(c("Gore test with multiple observations per cell"))

    ans <- list(p.value = PVAL,
                statistic = c("Gore chi-squared" = S),
                parameter = c("df" = df),
                method = METHOD,
                data.name = DNAME
                )
    class(ans) <- "htest"
    ans
}
