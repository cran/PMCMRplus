## MTest.R
## Part of the R package: PMCMRplus
##
## Copyright (C) 2017 Thorsten Pohlert
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

#' @name MTest
#' @title Extended One-Sided Studentised Range Test
#'
#' @description Performs Nashimoto-Wright's extended one-sided studentised range
#' test against an ordered alternative for normal data
#' with equal variances.
#'
#' This test is an extension of Hayter's OSRT
#' (see \code{\link{osrtTest}} by
#' applying a simple order restriction of
#' \eqn{\mu_{m'} - \mu_m \le \mu_j - \mu_i \le \
#' \mu_{l'} - \mu_{l}} for any \eqn{l \le i \le m}
#' and \eqn{m' \le j \le l'}. It tests all-pairs
#' \eqn{\mathrm{H}_{ij}: \mu_i \ge \mu_j} against
#' \eqn{\mathrm{A}_{ij}: \mu_i < \mu_j$ for any $1 \le i < j \le k}.
#'
#' @template class-PMCMR
#' 
#' @references
#' Nashimoto, K., Wright, F.T. (2005),
#' Multiple comparison procedures for detecting differences
#' in simply ordered means. \emph{Comput. Statist. Data Anal.} 48, 291--306.
#' @keywords htest
#' @concept AllPairsComparison
#' @importFrom stats ptukey
#' @examples
#' MTest(weight ~ group, data = PlantGrowth)
#' @export
MTest <- function(x, ...) UseMethod("MTest")

#' @rdname MTest
#' @aliases MTest.default
#' @method MTest default
#' @template one-way-parms
#' @importFrom stats ptukey
#' @importFrom stats var
#' @importFrom stats complete.cases
#' @export
MTest.default <-
function(x, g, ...){
        ## taken from stats::kruskal.test
        
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

    ## prepare tukey test
    ni <- tapply(x, g, length)
    n <- sum(ni)
    xi <- tapply(x, g, mean)
    s2i <- tapply(x, g, var)
    df <- n - k
    s2in <- 1 / df * sum(s2i * (ni - 1))

    sigma <- sqrt(s2in)
  
    STAT <- matrix(NA, ncol=k-1, nrow=k-1)
    for (i in 1:(k-1)){
        for(j in (i+1):k){
            u <- j
            m <- i:(u-1)
            tmp <- sapply(m, function(m) (xi[u] - xi[m]) / 
                                         (sigma / sqrt(2) *
                                          sqrt(1 / ni[m] + 1 /ni[u])))
            STAT[j-1,i] <- max(tmp)
        }
    }

    colnames(STAT) <- levels(g)[1:(k-1)]
    rownames(STAT) <- levels(g)[2:k]
    
    PVAL <- ptukey(STAT, nmeans = k, df = df, lower.tail=FALSE)

    colnames(PVAL) <- colnames(STAT)
    rownames(PVAL) <- rownames(STAT)
        
    METHOD <- "Nashimoto-Wright M-Test for ordered means \n\t\t of normal data with equal variance"
    MODEL <- data.frame(x, g)
    DIST <- "q"

    ans <- list(method = METHOD,
                data.name = DNAME,
                p.value = PVAL,
                statistic = STAT,
                p.adjust.method = "single-step",
                model = MODEL,
                dist = DIST,
                alternative = "greater")
    class(ans) <- "PMCMR"
    ans
}

#' @rdname MTest
#' @aliases MTest.formula
#' @method MTest formula
#' @template one-way-formula
#' @export
MTest.formula <-
function(formula, data, subset, na.action, ...)
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
    names(mf) <- NULL
    y <- do.call("MTest", c(as.list(mf)))
    y$data.name <- DNAME
    y
}
