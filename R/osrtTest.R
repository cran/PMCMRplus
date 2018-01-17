#  Copyright (C) 2017 Thorsten Pohlert
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#
#' @title One-Sided Studentised Range Test
#'
#' @description Performs Hayter's one-sided studentised range
#' test against an ordered alternative for normal data
#' with equal variances.
#'
#' @name osrtTest
#' @aliases osrtTest
#'
#' @template class-htest
#' 
#' @references
#' A. J. Hayter (1990) A One-Sided Studentised Range
#' Test for Testing Against a Simple Ordered Alternative.
#' \emph{Journal of the American Statistical Association}
#' 85, 778--785.
#'
#' @keywords htest
#' @importFrom stats ptukey
#' @examples
#' osrtTest(weight ~ group, data = PlantGrowth)
#' @export
osrtTest <- function(x, ...) UseMethod("osrtTest")

#' @rdname osrtTest
#' @method osrtTest default
#' @aliases osrtTest.default
#' @template one-way-parms
#' @export
osrtTest.default <-
function(x, g, ...)
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
    
    xi <- tapply(x, g, mean)
    ni <- tapply(x, g, length)
    k <- nlevels(g)
    n <- length(x)
    df <- n - k
    
    sigma2 <- 0
    c <- 0
    for (i in 1:k){
        for (j in 1:ni[i]){
            c <- c + 1
            sigma2 <- sigma2 + (x[c] - xi[i])^2 / df
	}
    }
    sigma <- sqrt(sigma2)
    
    compare.stats <- function(j,i) {
        dif <- xi[j] - xi[i] 
        A <- sigma / sqrt(2) * sqrt(1 / ni[i] + 1 / ni[j])
        qval <- abs(dif) / A
        return(qval)
    }
    
    val <- pairwise.table(compare.stats,levels(g), p.adjust.method="none" )
    STAT <- max(val, na.rm=TRUE)
    PVAL <- ptukey(STAT, nmeans = k, df = df, lower.tail=FALSE)
    METHOD <- "Hayter's One-Sided Studentised Range Test"
    ans <- list(method = METHOD,
                statistic = c(q = STAT),
                p.value = PVAL,
                data.name = DNAME,
                alternative = "increasing")
    class(ans) <- "htest"
    ans
}

#' @rdname osrtTest
#' @aliases osrtTest.formula
#' @method osrtTest formula
#' @template one-way-formula
#' @export
osrtTest.formula <-
    function(formula, data, subset, na.action,  ...)
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
    y <- do.call("osrtTest", c(as.list(mf)))
    y$data.name <- DNAME
    y
}
