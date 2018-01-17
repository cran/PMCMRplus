##  jonckheere.test.R
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
##
## Note: 
## jonckheere.test(x, g, alternative = "two.sided", continuity = TRUE)
##
## is equivalent to:
## cor.test(x, g, method = "kendall",
## alternative = "two.sided", continuity = TRUE)
##
#' @name jonckheereTest
#' @title Testing against Ordered Alternatives (Jonckheere-Terpstrata Test)
#' 
#' @description
#' Performs the Jonckheere-Terpstrata test for testing against ordered alternatives.
#' @details
#' The null hypothesis, H\eqn{_0: \theta_1 = \theta_2 = \ldots = \theta_k}
#' is tested against a simple order hypothesis,
#' H\eqn{_\mathrm{A}: \theta_1 \le \theta_2 \le \ldots \le
#' \theta_k,~\theta_1 < \theta_k}.
#'
#' The p-values are estimated from the standard normal distribution.
#' 
#' @note
#' \code{jonckheere.test(x, g, alternative = "two.sided", continuity = TRUE)} is
#' equivalent to
#' 
#' \code{cor.test(x, as.numeric(g), method = "kendall", alternative = "two.sided", continuity = TRUE)}
#' 
#' @section Source:
#' The code for the computation of the standard deviation
#' for the Jonckheere-Terpstrata test in the presence of ties was taken from:\cr
#' 
#' Kloke, J., McKean, J. (2016).
#' npsm: Package for Nonparametric Statistical Methods using R.
#' R package version 0.5. \url{https://CRAN.R-project.org/package=npsm}
#'
#' @references
#' Jonckheere, A. R. (1954). A distribution-free k-sample test
#' against ordered alternatives. \emph{Biometrica}, 41, 133–145.
#'
#' Kloke, J., McKean, J. W. (2015).
#' \emph{Nonparametric statistical methods using R}.
#' Boca Raton, FL: Chapman & Hall/CRC.
#' @template class-htest
#' @template trendTests
#' @export jonckheereTest
jonckheereTest <- function(x, ...) UseMethod("jonckheereTest")

#' @rdname jonckheereTest
#' @method jonckheereTest default
#' @aliases jonckheereTest.default
#' @template one-way-parms
#' @param alternative the alternative hypothesis. Defaults to \code{"two.sided"}.
#' @param continuity logical indicator whether a continuity correction
#' shall be performed. Defaults to \code{FALSE}.
#' @importFrom stats pnorm complete.cases
#' @export
jonckheereTest.default <-
    function(x, g, alternative = c("two.sided", "greater", "less"),
             continuity = FALSE, ...)
{
    if (is.list(x)) {
        if (length(x) < 2L)
            stop("'x' must be a list with at least 2 elements")
        DNAME <- deparse(substitute(x))
        x <- lapply(x, function(u) u <- u[complete.cases(u)])
        k <- length(x)
        l <- sapply(x, "length")
        if (any(l == 0))
            stop("all groups must contain data")
        g <- factor(rep(1 : k, l))
        ## check incoming from formula
        if(is.null(x$alternative)){
            alternative <- "two.sided"
        } else { 
            alternative <- x$alternative
        }
        if(is.null(x$continuity)) {
            continuity <- FALSE
        } else {
            continuity <- TRUE
        }
        x <- unlist(x)
		
    }
    else {
        if (length(x) != length(g))
            stop("'x' and 'g' must have the same length")
        DNAME <- paste(deparse(substitute(x)), "and",
                       deparse(substitute(g)))
        OK <- complete.cases(x, g)
        x <- x[OK]
        g <- g[OK]
        if (!all(is.finite(g)))
            stop("all group levels must be finite")
        g <- factor(g)
        k <- nlevels(g)
        if (k < 2)
            stop("all observations are in the same group")
    }
    if (!is.logical(continuity))
        stop("'continuity' must be 'FALSE' or 'TRUE'")
    alternative <- match.arg(alternative)
    n <- length(x)
    if (n < 2)
        stop("not enough observations")
    nij <- tapply(x, g, length)
    X <- matrix(NA, ncol= k, nrow = max(nij))
    j <- 0
    for (i in 1:k) {
        for (l in 1:nij[i]) {
            j = j + 1
            X[l,i] <- x[j]
        }
    }
    psi.f <- function(u) {
        psi <- (sign(u) + 1) / 2
        psi
    }
    Uij <- function(i,j,X){
        ni <- nij[i]
        nj <- nij[j]
        sumUij <- 0
        for (s in (1:ni)) {
            for (t in (1:nj)) {
                sumUij <- sumUij + psi.f(X[t,j] - X[s,i])
            }
        }
        sumUij
    }
    J <- 0
    for (i in (1:(k-1))) {
        for (j in ((i+1):k)) {
            J = J + Uij(i,j,X)
        }
    }
    mu <- (n^2 - sum(nij^2)) / 4
    st <- 0
    for (i in (1:k)) {
        st <- st + nij[i]^2 * (2 * nij[i] + 3)
        st
    }
	### check for ties
    TIES <- FALSE
    TIES <- (sum(table(rank(x)) - 1) > 0)
	
    if(!TIES){
        s <- sqrt((n^2 * (2 * n + 3) - st) / 72)
    
        S <- J - mu
 
    } else {
        ## if ties are present, no continuity correction will be done
        warning("Ties are present. Jonckheere z was corrected for ties.")
        S <- J - mu
        ## n : total sample size, nij group sizes
        ##
        ## taken from Kloke and McKean
        ## function jonckheere of package npsm 
        ## see citation("npsm")
        ##
        nt <- as.vector(table(x))
        s <- sqrt((n * (n - 1) * (2 * n + 5) - sum(nij * (nij - 1) * (2 * nij + 5)) -
        sum(nt * (nt - 1) * (2 * nt + 5))) /
        72 + (sum(nij * (nij - 1) * (nij - 2)) * sum(nt * (nt - 1) * (nt - 2)))/
        (36 * n * (n - 1) * (n - 2)) + (sum(nij * (nij - 1)) *
        sum(nt * (nt - 1)))/(8 * n * (n - 1)))
    }
	
    ## Check for continuity correction
    ## like in Kendall's tau
    if (continuity){
        S <- sign(S) * (abs(S) - 0.5)
    }
    STATISTIC <- S / s
    if (alternative == "two.sided") {
        PVAL <- 2 * min(pnorm(abs(STATISTIC), lower.tail = FALSE), 0.5)		
    } else if (alternative == "greater") {
        PVAL <- pnorm(STATISTIC, lower.tail = FALSE)
    } else {
        PVAL <- pnorm(STATISTIC)
    }
    ESTIMATES <- J 
    names(ESTIMATES) <- "JT"
    names(STATISTIC) <- "Jonckheere z-value"
    RVAL <- list(statistic = STATISTIC,
                 p.value = PVAL,
                 method = "Jonckheere-Terpstrata test",
                 data.name = DNAME,
                 alternative = alternative,
                 estimates = ESTIMATES)
    class(RVAL) <- "htest"
    return(RVAL)
}

#' @rdname jonckheereTest
#' @method jonckheereTest formula
#' @aliases jonckheereTest.formula
#' @template one-way-formula
#' @export
jonckheereTest.formula <-
function(formula, data, subset, na.action, alternative = c("two.sided", "greater", "less"), 
         continuity = FALSE, ...)
{
    mf <- match.call(expand.dots=FALSE)
    m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(stats::model.frame)
                 
    if(missing(formula) || (length(formula) != 3L))
        stop("'formula' missing or incorrect")
    mf <- eval(mf, parent.frame())  
    if(length(mf) > 2L)
        stop("'formula' should be of the form response ~ group")
    DNAME <- paste(names(mf), collapse = " by ")
    alternative <- match.arg(alternative)
    names(mf) <- NULL
    y <- do.call("jonckheereTest",
                 c(as.list(mf), alternative = alternative, continuity = continuity))
    y$data.name <- DNAME
    y
}
