# kruskal_R
# Part of the R package: PMCMR
#
# Copyright (C) 2017-2020 Thorsten Pohlert
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
#   Uses pKruskalWallis of package SuppDists
#
#   source:
#   W. J. Conover, R. L. Iman (1981) Rank transformations
#   as a bridge between parametric and nonparametric statistics,
#   The American Statistician 35 (3), 124--129.

#' @title Kruskal-Wallis Rank Sum Test
#' @description
#'  Performs a Kruskal-Wallis rank sum test.
#' @details
#' For one-factorial designs with non-normally distributed
#' residuals the Kruskal-Wallis rank sum test can be performed to test
#' the H\eqn{_0: F_1(x) = F_2(x) = \ldots = F_k(x)} against
#' the H\eqn{_\mathrm{A}: F_i (x) \ne F_j(x)~ (i \ne j)} with at least
#' one strict inequality.
#'
#' As the Kruskal-Wallis H-statistic is assymptotically
#' chi-squared distributed with \eqn{v = k - 1} degree
#' of freedom, the default test distribution is consequently
#' \code{dist = "Chisquare"}. If \code{dist = "KruskalWallis"} is selected,
#' an incomplete beta approximation is used for the calculation
#' of p-values as implemented in the function
#' \code{\link[SuppDists]{pKruskalWallis}} of the package
#' \pkg{SuppDists}. For \code{dist = "FDist"}
#' the proposed method of Conover and Imam (1981) is used, which is
#' equivalent to a one-way ANOVA F-test using rank transformed data
#' (see examples).
#'
#' @inherit friedmanTest references
#'
#' @seealso
#'  \code{\link{kruskal.test}}, \code{\link[SuppDists]{pKruskalWallis}},
#'  \code{\link{Chisquare}}, \code{\link{FDist}}
#'
#' @template class-htest
#'
#' @example examples/kSamples.R
#' @export
kruskalTest <- function(x, ...) UseMethod("kruskalTest")

#' @rdname kruskalTest
#' @method kruskalTest default
#' @template one-way-parms
#' @param dist the test distribution. Defaults's to \code{"Chisquare"}.
#' @importFrom stats pchisq pf
#' @importFrom SuppDists pKruskalWallis
#' @export
kruskalTest.default <-
    function(x, g, dist=c("Chisquare", "KruskalWallis", "FDist"), ...)
{
        ## taken from stats::kruskalTest
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
        if(is.null(x$dist)){
            dist <- "Chisquare"
        } else {
            dist <- x$dist
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

    dist <- match.arg(dist)
    x.rank <- rank(x)
    R.bar <- tapply(x.rank, g, mean,na.rm=T)
    R.n <- tapply(!is.na(x), g, length)
    k <- nlevels(g)
    n <- sum(R.n)

    getties <- function(x) {
        n <- length(x)
        t <- table(x)
        C <- 1 - sum(t^3 - t) / (n^3 - n)
        C <- min(1, C)
        return(C)
    }


    C <- getties(x.rank)
    if (C != 1) warning("Ties are present. Quantiles were corrected for ties.")
    ## Kruskal-Wallis statistic
    H <- (12 / (n * (n + 1))) *
        sum(tapply(x.rank, g, "sum")^2 / R.n) - 3 * (n + 1)
    PSTAT <- H / C

    if (dist == "Chisquare"){
        PARMS <- k - 1
        PVAL <- pchisq(PSTAT, df = PARMS, lower.tail = FALSE)
        names(PSTAT) <- "chi-squared"
        names(PARMS) <- "df"
    } else if (dist == "KruskalWallis"){

        ## pKruskalWallis from package SuppDists
        U <- sum(1 / R.n)
        c <- k
        N <- n
        PARMS <- c(c, U, N)
        PVAL <- pKruskalWallis(PSTAT, c = c,
                               N = N, U = U,
                               lower.tail = FALSE)
        names(PSTAT) <- "H"
        names(PARMS) <- c("k", "U", "N")

    } else {
        ## F distribution
        N <- n
        H <- PSTAT
        df1 <- k - 1
        df2 <- N - k
        PSTAT <- (H / ( k - 1)) / (( N - 1 - H) / (N - k))
        PVAL <- pf(PSTAT, df1 = df1, df2=df2, lower.tail = FALSE)
        PARMS <- c(df1, df2)
        names(PARMS) <- c("num df", "denom df")
        names(PSTAT) <- "Conover's F"
    }

    METHOD <- paste("Kruskal-Wallis test")

    ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL,
                statistic = PSTAT, parameter = PARMS)
    class(ans) <- "htest"
    ans
}

#' @rdname kruskalTest
#' @method kruskalTest formula
#' @template one-way-formula
#' @export
kruskalTest.formula <-
    function(formula, data, subset, na.action,
             dist=c("Chisquare", "KruskalWallis", "FDist"), ...)
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
    dist <- match.arg(dist)
    names(mf) <- NULL
    y <- do.call("kruskalTest", c(as.list(mf), dist = dist))
    y$data.name <- DNAME
    y
}
