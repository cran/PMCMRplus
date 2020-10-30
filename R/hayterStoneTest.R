#  Copyright (C) 2020 Thorsten Pohlert
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
#' @name hayterStoneTest
#' @aliases hayterStoneTest
#' @title Hayter-Stone Test
#'
#' @description Performs the non-parametric Hayter-Stone procedure
#' to test against an monotonically increasing alternative.
#'
#' @details
#' Let \eqn{X} be an identically and idepentendly distributed variable
#' that was \eqn{n} times observed at \eqn{k} increasing treatment levels.
#' Hayter and Stone (1991) proposed a non-parametric procedure
#' to test the null hypothesis, H: \eqn{\theta_i = \theta_j ~~ (i < j \le k)}
#' against a simple order alternative, A: \eqn{\theta_i < \theta_j}, with at least
#' one inequality being strict.
#'
#' The statistic for a global test is calculated as,
#' \deqn{
#'  h = \max_{1 \le i < j \le k} \frac{2 \sqrt{6} \left(U_{ij} - n_i n_j / 2 \right)}
#'  {\sqrt{n_i n_j \left(n_i + n_j + 1 \right)}},
#' }{%
#' SEE PDF.
#' }
#'
#' with the Mann-Whittney counts:
#' \deqn{
#' U_{ij} =  \sum_{a=1}^{n_i} \sum_{b=1}^{n_j} I\left\{x_{ia} < x_{ja}\right\}.
#' }{%
#' SEE PDF
#' }
#'
#' Under the large sample approximation, the test statistic \eqn{h} is distributed
#' as \eqn{h_{k,\alpha,v}}. Thus, the null hypothesis is rejected, if \eqn{h > h_{k,\alpha,v}}, with \eqn{v = \infty}
#' degree of freedom.
#'
#' If \code{method = "look-up"} the function will not return
#' p-values. Instead the critical h-values
#' as given in the tables of Hayter (1990) for
#' \eqn{\alpha = 0.05} (one-sided)
#' are looked up according to the number of groups (\eqn{k}) and
#' the degree of freedoms (\eqn{v = \infty}).
#'
#' If \code{method = "boot"} an asymetric permutation test
#' is conducted and a \eqn{p}-value is returned.
#'
#' @return
#' Either a list of class \code{htest} or a
#' list with class \code{"osrt"} that contains the following
#' components:
#' @template returnOsrt
#'
#' @references
#' Hayter, A. J.(1990) A One-Sided Studentised Range
#' Test for Testing Against a Simple Ordered Alternative,
#' \emph{J Amer Stat Assoc}
#' \bold{85}, 778--785.
#'
#' Hayter, A.J., Stone, G. (1991)
#' Distribution free multiple comparisons for monotonically ordered treatment effects.
#' \emph{Austral J Statist} \bold{33}, 335--346.
#'
#' @seealso
#' \code{\link{osrtTest}}, \code{\link{hsAllPairsTest}},
#' \code{\link{sample}}
#'
#' @keywords htest nonparametric
#'
#' @example examples/shirleyEx.R
#' @export
hayterStoneTest <- function(x, ...) UseMethod("hayterStoneTest")

#' @rdname hayterStoneTest
#' @method hayterStoneTest default
#' @aliases hayterStoneTest.default
#' @template one-way-parms
#' @param alternative the alternative hypothesis. Defaults to \code{greater}.
#' @param method a character string specifying the test statistic to use.
#' Defaults to \code{"look-up"} that uses published Table values.
#' @param nperm number of permutations for the asymptotic permutation test.
#' Defaults to \code{1000}. Ignored, if \code{method = "look-up"}.
#' @export
hayterStoneTest.default <-
function(x, g, alternative = c("greater", "less"),
         method = c("look-up", "boot"),
         nperm = 1E4,...)
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

    method <- match.arg(method)
    alternative <- match.arg(alternative)
    if (alternative == "less") {
        x <- -x
    }

    ## Mann-Whittney counts
    U <- function(xi, xj) {
        S <- 0
        for (i in seq_along(xi)) {
            S <- S + sum(ifelse(xi[i] < xj, 1, 0))
        }
        S
    }



    n <- tapply(x, g, length)

    hsStat <- function(x, ix, g, n) {
        x <- x[ix]
        k <- nlevels(g)
        ## Statistic
        Sij <- numeric(k * (k - 1) / 2)
        l <- 0
        SQRT6 <- sqrt(6)
        for (i in 1:(k - 1)) {
            G <- levels(g)[i]
            ok <- g == G
            xi <- x[ok]
            for (j in (i + 1):k) {
                l <- l + 1
                ## select xj
                G <- levels(g)[j]
                ok <- g == G
                xj <- x[ok]
                Uij <- U(xi, xj)

                Sij[l] <- 2 * SQRT6 * (Uij - n[i] * n[j] / 2) /
                    (sqrt(n[i] * n[j] * (n[i] + n[j] + 1)))

            }
        }

        STAT <- max(Sij)
        return(STAT)
    }

    k <- nlevels(g)
    N <- sum(n)
    l <- 1:N
    STAT <- hsStat(x, l, g, n)

    ## for output description
    METHOD <- "Hayter-Stone Test"
    parameter = c(k, Inf)
    names(parameter) <- c("k", "df")

    if (method == "boot") {
        ## permutation
        mt <- matrix(NA, ncol = 1, nrow = nperm)
        for (i in 1:nperm) {
            ix <- sample(l)
            mt[i, ] <- hsStat(x, ix, g, n)
        }

        ## pvalues
        PVAL <- sapply(1:1, function(j) {
            p <- (sum(mt[, j] >= STAT[j])) / nperm
            p
        })

        ans <- list(
            statistic = c(h = STAT),
            p.value = PVAL,
            method = METHOD,
            data.name = DNAME,
            alternative = alternative,
            parameter = parameter
        )
        class(ans) <- "htest"
        return(ans)

    } else {

    ## look for k and v = inf
    nrows <- nrow(TabCrit$hayter.h005)
    kk <- as.numeric(colnames(TabCrit$hayter.h005))
    ## check for kk
    if (k > max(kk) | k < min(kk)) stop("No critical values for k = ", k)

    hCrit <- unlist(TabCrit$hayter.h005[nrows, paste0(k)])

    ans <- list(
        method = METHOD,
        data.name = DNAME,
        crit.value = hCrit,
        statistic = STAT,
        parameter = parameter,
        alternative = alternative,
        dist = "h"
    )
    class(ans) <- "osrt"
    return(ans)
    }
}

#' @rdname hayterStoneTest
#' @aliases hayterStoneTest.formula
#' @method hayterStoneTest formula
#' @template one-way-formula
#' @export
hayterStoneTest.formula <-
    function(formula, data, subset, na.action,
             alternative = c("greater", "less"),
             method = c("look-up", "boot"),
             nperm = 1E4,
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
    method <- match.arg(method)
    names(mf) <- NULL
    y <- do.call("hayterStoneTest", c(as.list(mf),
                               alternative = alternative,
                               method = method,
                               nperm = nperm))
    y$data.name <- DNAME
    y
}
