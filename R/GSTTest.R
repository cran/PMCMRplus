# GSTTest.R
# Part of the R package: PMCMR
#
# Copyright (C) 2018 Thorsten Pohlert
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

#' @name GSTTest
#' @title Generalized Siegel-Tukey Test of Homogeneity of
#' Scales
#' @description
#'  Performs a Siegel-Tukey k-sample rank dispersion test.
#' @details
#' Meyer-Bahlburg (1970) has proposed a generalized Siegel-Tukey
#' rank dispersion test for the \eqn{k}-sample case.
#' Likewise to the \code{\link{fligner.test}}, this test
#' is a nonparametric test for testing the homogegeneity of
#' scales in several groups.
#' Let \eqn{\theta_i}{theta_i}, and \eqn{\lambda_i}{lambda_i} denote
#'  location and scale parameter of the \eqn{i}th group,
#'  then for the two-tailed case, the null hypothesis
#'  H: \eqn{\lambda_i / \lambda_j = 1 | \theta_i = \theta_j, ~ i \ne j}{%
#'  lambda_i / lambda_j = 1 | theta_i = theta_j, i != j} is
#'  tested against the alternative,
#'  A: \eqn{\lambda_i / \lambda_j \ne 1}{lambda_i / lambda_j != 1}
#'  with at least one inequality beeing strict.
#'
#'  The data are combinedly ranked according to Siegel-Tukey.
#'  The ranking is done by alternate extremes (rank 1 is lowest,
#'  2 and 3 are the two highest, 4 and 5 are the two next lowest, etc.).
#'
#' Meyer-Bahlburg (1970) showed, that the Kruskal-Wallis H-test
#' can be employed on the Siegel-Tukey ranks.
#' The H-statistic is assymptotically
#' chi-squared distributed with \eqn{v = k - 1} degree
#' of freedom, the default test distribution is consequently
#' \code{dist = "Chisquare"}. If \code{dist = "KruskalWallis"} is selected,
#' an incomplete beta approximation is used for the calculation
#' of p-values as implemented in the function
#' \code{\link[SuppDists]{pKruskalWallis}} of the package
#' \pkg{SuppDists}.
#'
#' @note
#' If ties are present, a tie correction is performed and
#' a warning message is given. The GSTTest is sensitive to
#' median differences, likewise to the Siegel-Tukey test.
#' It is thus appropriate to apply this test on the residuals
#' of a one-way ANOVA, rather than on the original data
#' (see example).
#'
#' @references
#' H.F.L. Meyer-Bahlburg (1970), A nonparametric test for relative
#' spread in k unpaired samples, \emph{Metrika} \bold{15}, 23--29.
#'
#' @seealso
#'  \code{\link{fligner.test}}, \code{\link[SuppDists]{pKruskalWallis}},
#'  \code{\link{Chisquare}}, \code{\link{fligner.test}}
#'
#' @template class-htest
#'
#' @examples
#' GSTTest(count ~ spray, data = InsectSprays)
#'
#' ## as means/medians differ, apply the test to residuals
#' ## of one-way ANOVA
#' ans <- aov(count ~ spray, data = InsectSprays)
#' GSTTest( residuals( ans) ~ spray, data =InsectSprays)
#'
#' @keywords htest
#' @keywords nonparametric
#' @export
GSTTest <- function(x, ...) UseMethod("GSTTest")

#' @rdname GSTTest
#' @method GSTTest default
#' @template one-way-parms
#' @param dist the test distribution. Defaults's to \code{"Chisquare"}.
#' @importFrom stats pchisq
#' @importFrom SuppDists pKruskalWallis
#' @export
GSTTest.default <-
    function(x,
             g,
             dist=c("Chisquare", "KruskalWallis"),
             ...)
{
        ## taken from stats::GSTTest
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

    ## Siegel-Tukey ranking
    ## rank sequence
    ## based on code from Daniel Malter
    ord <- order(x)
    gord <- g[ord]
    n <- tapply(x[ord], gord, length)
    N <- sum(n)


    a <- rep(seq(ceiling(N / 4)), each=2)
    b <- rep(c(0, 1), ceiling(N)/4)
    suppressWarnings(
      rk.up <- c(1, (a * 4 + b))[1:ceiling(N / 2)]
    )
    suppressWarnings(
      rk.down <- rev(c(a * 4 + b - 2)[1:floor(N / 2)])
    )

    r <- c(rk.up, rk.down)
    T <- tapply(r, gord, sum)

    if(sum(T) != N * (N + 1) /2) {
      warning("Does not sum up to check sum!")
    }

    getties <- function(x) {
        n <- length(x)
        t <- table(x)
        C <- 1 - sum(t^3 - t) / (n^3 - n)
        C <- min(1, C)
        return(C)
    }


    C <- getties(x)
    if (C != 1) warning("Ties are present. Quantiles were corrected for ties.")
    ## Kruskal-Wallis statistic
    H <- (12 / (N * (N + 1))) *
        sum(T * T / n) - 3 * (N + 1)
    PSTAT <- H / C

    if (dist == "Chisquare"){
        PARMS <- k - 1
        PVAL <- pchisq(PSTAT, df = PARMS, lower.tail = FALSE)
        names(PSTAT) <- "chi-squared"
        names(PARMS) <- "df"
    } else if (dist == "KruskalWallis"){

        ## pKruskalWallis from package SuppDists
        U <- sum(1 / n)
        c <- k
        PARMS <- c(c, U, N)
        PVAL <- pKruskalWallis(PSTAT, c = c,
                               N = N, U = U,
                               lower.tail = FALSE)
        names(PSTAT) <- "H"
        names(PARMS) <- c("k", "U", "N")

    }

    METHOD <- paste("Generalized Siegel-Tukey test of homogeneity of scales")

    ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL,
                statistic = PSTAT, parameter = PARMS)
    class(ans) <- "htest"
    ans
}

#' @rdname GSTTest
#' @method GSTTest formula
#' @template one-way-formula
#' @export
GSTTest.formula <-
    function(formula, data, subset, na.action,
             dist=c("Chisquare", "KruskalWallis"),
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
    dist <- match.arg(dist)
    names(mf) <- NULL
    y <- do.call("GSTTest", c(as.list(mf), dist = dist))
    y$data.name <- DNAME
    y
}
