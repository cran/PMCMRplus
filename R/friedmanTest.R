## friedmanTest.R
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
#' @title Friedman Rank Sum Test
#' @description
#' Performs a Friedman rank sum test. The null hypothesis
#' H\eqn{_0: \theta_i = \theta_j~~(i \ne j)} is tested against the
#' alternative H\eqn{_{\mathrm{A}}: \theta_i \ne \theta_j}, with at least
#' one inequality beeing strict.
#' @name friedmanTest
#' @template class-htest
#' @details
#' The function has implemented Friedman's test as well as
#' the extension of Conover anf Iman (1981). Friedman's
#' test statistic is assymptotically chi-squared distributed.
#' Consequently, the default test distribution is \code{dist = "Chisquare"}.
#'
#' If \code{dist = "FDist"} is selected, than the approach of
#' Conover and Imam (1981) is performed.
#' The Friedman Test using the \eqn{F}-distribution leads to
#' the same results as doing an two-way Analysis of Variance without
#' interaction on rank transformed data.
#'
#' @references
#'  Conover, W. J., Iman, R. L. (1981) Rank transformations as a bridge
#'  between parametric and nonparametric statistics, \emph{The American
#'    Statistician} \bold{35}, 124--129.
#'
#'  Sachs, L. (1997) \emph{Angewandte Statistik}. Berlin: Springer.
#'
#' @keywords htest nonparametric
#' @seealso
#' \code{\link{friedman.test}}
#' @examples
#' ## Hollander & Wolfe (1973), p. 140ff.
#' ## Comparison of three methods ("round out", "narrow angle", and
#' ##  "wide angle") for rounding first base.  For each of 18 players
#' ##  and the three method, the average time of two runs from a point on
#' ##  the first base line 35ft from home plate to a point 15ft short of
#' ##  second base is recorded.
#' RoundingTimes <-
#' matrix(c(5.40, 5.50, 5.55,
#'         5.85, 5.70, 5.75,
#'         5.20, 5.60, 5.50,
#'         5.55, 5.50, 5.40,
#'         5.90, 5.85, 5.70,
#'         5.45, 5.55, 5.60,
#'         5.40, 5.40, 5.35,
#'         5.45, 5.50, 5.35,
#'         5.25, 5.15, 5.00,
#'         5.85, 5.80, 5.70,
#'         5.25, 5.20, 5.10,
#'         5.65, 5.55, 5.45,
#'         5.60, 5.35, 5.45,
#'         5.05, 5.00, 4.95,
#'         5.50, 5.50, 5.40,
#'         5.45, 5.55, 5.50,
#'         5.55, 5.55, 5.35,
#'         5.45, 5.50, 5.55,
#'         5.50, 5.45, 5.25,
#'         5.65, 5.60, 5.40,
#'         5.70, 5.65, 5.55,
#'         6.30, 6.30, 6.25),
#'       nrow = 22,
#'       byrow = TRUE,
#'       dimnames = list(1 : 22,
#'                       c("Round Out", "Narrow Angle", "Wide Angle")))
#'
#' ## Chisquare distribution
#' friedmanTest(RoundingTimes)
#'
#' ## check with friedman.test from R stats
#' friedman.test(RoundingTimes)
#'
#' ## F-distribution
#' friedmanTest(RoundingTimes, dist = "FDist")
#'
#' ## Check with One-way repeated measure ANOVA
#' rmat <- RoundingTimes
#' for (i in 1:length(RoundingTimes[,1])) rmat[i,] <- rank(rmat[i,])
#' dataf <- data.frame(
#'     y = y <- as.vector(rmat),
#'     g = g <- factor(c(col(RoundingTimes))),
#'     b = b <- factor(c(row(RoundingTimes))))
#' summary(aov(y ~ g + Error(b), data = dataf))
#'
#' @export
friedmanTest <- function(y, ...)
    UseMethod("friedmanTest")

#' @rdname friedmanTest
#' @aliases friedmanTest.default
#' @method friedmanTest default
#' @template two-way-parms
#' @param dist the test distribution. Defaults to \code{Chisquare}.
#' @importFrom stats pf
#' @importFrom stats pchisq
#' @export
friedmanTest.default <-
    function(y, groups, blocks, dist = c("Chisquare", "FDist"), ...){

    	dist <- match.arg(dist)
    	## 2019-10-16
    	## novel external function
    	ans <- frdRanks(y, groups, blocks)
    	mat <- ans$r
    	n <- nrow(mat)
    	k <- ncol(mat)
    	GRPNAMES <- colnames(mat)

        R.sum <- colSums(mat)
        METHOD <- paste("Friedman rank sum test")

        ## Friedman's T1 value
        A1 <- 0
        for (i in 1:n){
            for (j in 1:k){
                A1 <- A1 + mat[i,j]^2
            }
        }
        C1 <- (n * k * (k + 1)^2) / 4
        TT <- 0
        for (j in 1:k) {
            TT <- TT + (R.sum[j] - ((n * (k + 1))/2))^2
        }
        T1 <- ((k - 1) * TT) / (A1 - C1)

        if (dist == "Chisquare"){
            PARMS <- k - 1
            PVAL <- pchisq(T1, PARMS,lower.tail = FALSE)
            names(PARMS) <- "df"
            names(T1) <- "Friedman chi-squared"
            STAT <- T1
        } else {
        ## F-value
            T2 <- ((n - 1) * T1) / (n * (k - 1) - T1)
            df1 <- k - 1
            df2 <- (n-1) * (k-1)
            PARMS <- c(df1, df2)
            names(PARMS) <- c("num df", "denom df")
            names(T2) <- "Conover's F"
            PVAL <- pf(T2, df1, df2, lower.tail = FALSE)
            STAT <- T2
        }
        ans <- list(statistic=STAT, p.value = PVAL,
                    parameter= PARMS, method = METHOD,
                    data.name=ans$DNAME)
        class(ans) <- "htest"
        ans
    }
