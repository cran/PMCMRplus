## leTest.R
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
#' @name leTest
#' @title Testing against Ordered Alternatives (Le's Test)
#'
#' @description
#' Performs Le's test for testing against ordered alternatives.
#' @details
#' The null hypothesis, H\eqn{_0: \theta_1 = \theta_2 = \ldots = \theta_k}
#' is tested against a simple order hypothesis,
#' H\eqn{_\mathrm{A}: \theta_1 \le \theta_2 \le \ldots \le
#' \theta_k,~\theta_1 < \theta_k}.
#'
#' The p-values are estimated from the standard normal distribution.
#'
#' @template class-htest
#' @template trendTests
#' @references
#' Le, C. T. (1988) A new rank test against ordered alternatives
#' in k-sample problems, \emph{Biometrical Journal} \bold{30}, 87--92.
#' @export leTest
leTest <- function(x, ...) UseMethod("leTest")

#' @rdname leTest
#' @method leTest default
#' @aliases leTest.default
#' @template one-way-parms
#' @param alternative the alternative hypothesis. Defaults to \code{"two.sided"}.
#' @importFrom stats pnorm complete.cases
#' @export
leTest.default <-
    function(x, g, alternative = c("two.sided", "greater", "less"),...)
{
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
        if(!is.null(x$alternative)) alternative <- x$alternative
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

    ## test
    alternative <- match.arg(alternative)

    ## prepare
    k <- nlevels(g)
    n <- tapply(x, g, length)
    r <- rank(x)
    R <- tapply(r, g, mean)
    N <- sum(n)

    ## local variables
    lambda <- numeric(k)
    L <- numeric(k)
    M <- numeric(k)

    for (i in 2:k)
    {
        L[i] <- sum(n[1:(i-1)])
    }

    for (i in 1:(k-1))
    {
        M[i] <- sum(n[(i+1):k])
    }

    lambda <- n * (L - M)

    ## statistic
    W <- sum(lambda * R)

    sigma <-  sqrt((N * (N +1) / 12) * sum( n * (L - M)^2))

    ## z value
    z <- W / sigma

    if (alternative == "two.sided"){
        PVAL <- 2 * (min(0.5, pnorm(abs(z), lower.tail=FALSE)))
    } else if (alternative == "greater"){
        PVAL <- pnorm(z, lower.tail=FALSE)
    } else {
        PVAL <- pnorm(z)
    }


    METHOD <- paste("Le's rank test for ordered alternatives")

    ans <- list(method = METHOD,
                data.name = DNAME,
                p.value = PVAL,
                statistic = c(z = z),
                alternative = alternative,
                estimate = c(W = W))
    class(ans) <- "htest"
    ans
}

#' @rdname leTest
#' @method leTest formula
#' @aliases leTest.formula
#' @template one-way-formula
#' @export
leTest.formula <-
function(formula, data, subset, na.action, alternative = c("two.sided", "greater", "less"),
         ...)
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
    y <- do.call("leTest", c(as.list(mf), alternative = alternative))
    y$data.name <- DNAME
    y
}
