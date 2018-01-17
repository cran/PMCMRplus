# normalScoresTest.R
# Part of the R package: PMCMR
#
# Copyright (C) 2017 Thorsten Pohlert
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
#   Uses pNormScore and normOrder from package SuppDists
#
#   source: 
#   package SuppDists

#' @title Lu-Smith Normal Scores Test
#' @description Performs the Lu-Smith normal score test
#' @details
#' For one-factorial designs with non-normally distributed
#' residuals the Lu-Smith normal score test can be performed to test
#' the H\eqn{_0: F_1(x) = F_2(x) = \ldots = F_k(x)} against
#' the H\eqn{_\mathrm{A}: F_i (x) \ne F_j(x) ~ (i \ne j)} with at least
#' one strict inequality. This function is basically a wrapper function to
#' \code{\link[SuppDists]{pNormScore}} of the package \pkg{SuppDists}.
#' @name normalScoresTest
#' @template class-htest
#' @keywords htest nonparametric
#' @importFrom SuppDists normOrder
#' @importFrom SuppDists pNormScore
#' @seealso
#' \code{\link{vanWaerdenTest}}, \code{\link{kruskalTest}},
#' \code{\link[SuppDists]{pNormScore}}
#' @references
#' Lu, H., Smith, P. (1979). Distribution of normal scores statistic
#' for nonparametric one-way analysis of variance.
#' \emph{Journal of the American Statistical Association}, 74, 715--722.
#' @examples
#' normalScoresTest(count ~ spray, data = InsectSprays)
#' @export
normalScoresTest <- function(x, ...) UseMethod("normalScoresTest")

#' @rdname normalScoresTest
#' @aliases normalScoresTest.default
#' @method normalScoresTest default
#' @template one-way-parms
#' @export
normalScoresTest.default <-
    function(x, g, ...)
{
        ## taken from stats::normalScoresTest
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

    n <- length(x)
    ni <- tapply(x, g, length)
    rij <- rank(x, ties.method = "random")
    eN <- normOrder(n)
    eij <- eN[rij]
    Si <- tapply(eij, g, sum) 

    ## Statistic
    STAT <- (n - 1) * 1 / sum(eN^2) * sum(Si^2 / ni)
    U <- sum(1 / ni)

    ## pvalue
    PVAL <- pNormScore(STAT, k, n, U, lower.tail=FALSE)    
    PARMS <- c(k = k, n= n,U =U)
    METHOD <- paste("Normal Scores test")

    ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL,
                parameter = PARMS, statistic = c(Score = STAT))
    class(ans) <- "htest"
    ans
}

#' @rdname normalScoresTest
#' @aliases normalScoresTest.formula
#' @method normalScoresTest formula
#' @template one-way-formula
#' @export
normalScoresTest.formula <-
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
    y <- do.call("normalScoresTest", c(as.list(mf)))
    y$data.name <- DNAME
    y
}
