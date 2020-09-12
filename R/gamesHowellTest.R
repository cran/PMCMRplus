## gamesHowellTest.R
## Part of the R package: PMCMR
##
## Copyright (C) 2017-2020 Thorsten Pohlert
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
## Source:
## L. Sachs (1997) Angewandte Statistik, New York: Springer.
#' @name gamesHowellTest
#' @title Games-Howell Test
#' @description
#' Performs Games-Howell all-pairs comparison test for normally distributed
#' data with unequal group variances.
#'
#' @details
#' For all-pairs comparisons in an one-factorial layout
#' with normally distributed residuals but unequal between-groups variances
#' the Games-Howell Test can be performed. A total of \eqn{m = k(k-1)/2}
#' hypotheses can be tested. The null hypothesis
#' H\eqn{_{ij}: \mu_i(x) = \mu_j(x)} is tested in the two-tailed test
#' against the alternative
#' A\eqn{_{ij}: \mu_i(x) \ne \mu_j(x), ~~ i \ne j}.
#'
#' The p-values are computed from the studentized range distribution.
#' @seealso
#' \code{\link{ptukey}}
#' @keywords htest
#' @concept parametric
#'
#' @template class-PMCMR
#'
#' @examples
#' fit <- aov(weight ~ feed, chickwts)
#' shapiro.test(residuals(fit))
#' bartlett.test(weight ~ feed, chickwts) # var1 = varN
#' anova(fit)
#'
#' ## also works with fitted objects of class aov
#' res <- gamesHowellTest(fit)
#' summary(res)
#' summaryGroup(res)
#' @importFrom stats ptukey
#' @importFrom stats complete.cases
#' @importFrom stats var
#' @export
gamesHowellTest <- function(x, ...) UseMethod("gamesHowellTest")

#' @rdname gamesHowellTest
#' @method gamesHowellTest default
#' @aliases gamesHowellTest.default
#' @template one-way-parms-aov
#' @export
gamesHowellTest.default <-
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

    ## prepare gamesHowell test
    ni <- tapply(x, g, length)
    n <- sum(ni)
    xi <- tapply(x, g, mean)
    s2i <- tapply(x, g, var)

    s2in <- 1 / (n - k) * sum(s2i * (ni - 1))

    compare.stats <- function(i,j) {
        dif <- xi[i] - xi[j]
        A <- (s2i[i] / ni[i] + s2i[j] / ni[j])
        qval <- dif / sqrt(A) * sqrt(2)
        return(qval)
    }

    PSTAT <- pairwise.table(compare.stats,levels(g), p.adjust.method="none" )

    compare.levels <- function(i,j) {
        dif <- xi[i] - xi[j]
        A <- (s2i[i] / ni[i] + s2i[j] / ni[j])
        qval <- dif / sqrt(A) * sqrt(2)
        df <- A^2 / (s2i[i]^2 / (ni[i]^2 * (ni[i] - 1)) +
                    s2i[j]^2 / (ni[j]^2 * (ni[j] - 1)))
        pval <- ptukey(abs(qval), nmeans = k, df = df, lower.tail = FALSE)
        return(pval)
    }

    PVAL <-  pairwise.table(compare.levels, levels(g), p.adjust.method="none")

    MODEL <- data.frame(x, g)
    DIST <- "q"
    METHOD <- "Games-Howell test"
    ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL,
                statistic = PSTAT, p.adjust.method = "none",
                model = MODEL, dist = DIST, alternative = "two.sided")
    class(ans) <- "PMCMR"
    ans
}

#' @rdname gamesHowellTest
#' @method gamesHowellTest formula
#' @aliases gamesHowellTest.formula
#' @template one-way-formula
#' @export
gamesHowellTest.formula <-
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
    y <- do.call("gamesHowellTest", c(as.list(mf)))
    y$data.name <- DNAME
    y
}

#' @rdname gamesHowellTest
#' @aliases gamesHowellTest.aov
#' @method gamesHowellTest aov
# @param obj A fitted model object, usually an \link[stats]{aov} fit.
#' @export
gamesHowellTest.aov <- function(x, ...) {
    model <- x$model
    DNAME <- paste(names(model), collapse = " by ")
    names(model) <- c("x", "g")
    y <- do.call("gamesHowellTest", as.list(model))
    y$data.name <- DNAME
    y
}
