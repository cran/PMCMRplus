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
#' @name hsAllPairsTest
#' @aliases hsAllPairsTest
#' @title Hayter-Stone All-Pairs Comparison Test
#'
#' @description Performs the non-parametric Hayter-Stone all-pairs procedure
#' to test against monotonically increasing alternatives.
#'
#' @details
#' Let \eqn{X} be an identically and idepentendly distributed variable
#' that was \eqn{n} times observed at \eqn{k} increasing treatment levels.
#' Hayter and Stone (1991) proposed a non-parametric procedure
#' to test the null hypothesis, H: \eqn{\theta_i = \theta_j ~~ (i < j \le k)}
#' against a simple order alternative, A: \eqn{\theta_i < \theta_j}.
#'
#' The statistic for all-pairs comparisons is calculated as,
#' \deqn{
#'  S_{ij} = \frac{2 \sqrt{6} \left(U_{ij} - n_i n_j / 2 \right)}
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
#' Under the large sample approximation, the test statistic \eqn{S_{ij}} is distributed
#' as \eqn{h_{k,\alpha,v}}. Thus, the null hypothesis is rejected,
#' if \eqn{S_{ij} > h_{k,\alpha,v}}, with \eqn{v = \infty} degree of freedom.
#'
#' If \code{method = "look-up"} the function will not return
#' p-values. Instead the critical h-values
#' as given in the tables of Hayter (1990) for
#' \eqn{\alpha = 0.05} (one-sided)
#' are looked up according to the number of groups (\eqn{k}) and
#' the degree of freedoms (\eqn{v = \infty}).
#'
#' If \code{method = "boot"} an asymetric permutation test
#' is conducted and \eqn{p}-values are returned.
#'
#'  If \code{method = "asympt"} is selected the asymptotic
#' \eqn{p}-value is estimated as implemented in the
#' function \code{pHayStonLSA} of the package \pkg{NSM3}.
#'
#' @source
#' If \code{method = "asympt"} is selected, this function calls
#' an internal probability function \code{pHS}. The GPL-2 code for
#' this function was taken from \code{pHayStonLSA} of the
#' the package \pkg{NSM3}:
#'
#' Grant Schneider, Eric Chicken and Rachel Becvarik (2020) NSM3:
#' Functions and Datasets to Accompany Hollander, Wolfe, and
#' Chicken - Nonparametric Statistical Methods, Third Edition. R
#' package version 1.15. \url{https://CRAN.R-project.org/package=NSM3}
#'
#' @return
#' Either a list of class \code{"PMCMR"} or a
#' list with class \code{"osrt"} that contains the following
#' components:
#' @template returnOsrt
#' @template class-PMCMR
#'
#' @references
#' Hayter, A. J.(1990) A One-Sided Studentised Range
#' Test for Testing Against a Simple Ordered Alternative,
#' \emph{Journal of the American Statistical Association}
#' \bold{85}, 778--785.
#'
#' Hayter, A.J., Stone, G. (1991)
#' Distribution free multiple comparisons for monotonically ordered treatment effects.
#' \emph{Austral J Statist} \bold{33}, 335--346.
#'
#' @seealso
#' \code{\link{hayterStoneTest}} \code{\link{sample}}
#'
#' @keywords htest nonparametric
#'
#' @example examples/shirleyEx.R
#'
#' @export
hsAllPairsTest <- function(x, ...) UseMethod("hsAllPairsTest")

#' @rdname hsAllPairsTest
#' @method hsAllPairsTest default
#' @aliases hsAllPairsTest.default
#' @template one-way-parms
#' @param alternative the alternative hypothesis. Defaults to \code{greater}.
#' @param method a character string specifying the test statistic to use.
#' Defaults to \code{"look-up"} that uses published Table values of Williams (1972).
#' @param nperm number of permutations for the asymptotic permutation test.
#' Defaults to \code{1000}. Ignored, if \code{method = "look-up"}.
#' @export
hsAllPairsTest.default <-
function(x, g, alternative = c("greater", "less"),
         method = c("look-up", "boot", "asympt"),
         nperm = 1E4,
         ...)
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


    hStat <- function(x0, ix, g) {
        x <- x0[ix]
        k <- nlevels(g)
        n <- tapply(x, g, length)
        ## Statistic
        Sij <- matrix(NA, nrow = k - 1, ncol = k - 1)
        # Sij <- numeric(k * (k - 1) / 2)
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

                Sij[j - 1, i] <- 2 * SQRT6 * (Uij - n[i] * n[j] / 2) /
                    (sqrt(n[i] * n[j] * (n[i] + n[j] + 1)))

            }
        }
        return(Sij)
    }

    ## first call
    l <- seq_along(x)
    STAT <- hStat(x, l, g)

    k <- nlevels(g)
    #n <- tapply(x, g, length)
    #STAT <- Sij
    colnames(STAT) <- levels(g)[1:(k-1)]
    rownames(STAT) <- levels(g)[2:k]

    if (method == "boot" | method == "asympt") {


        if (method == "boot") {
            ## permutation
            hValue <- as.numeric(STAT)
            m <- length(hValue)
            mt <- matrix(NA, ncol = m, nrow = nperm)
            for (i in 1:nperm) {
                ix <- sample(l)
                tmp <- hStat(x, ix, g)
                mt[i, ] <- as.numeric(tmp)
            }

            ## pvalues
            PVAL <- sapply(1:m, function(j) {
                p <- sum(mt[, j] >= hValue[j]) / nperm
                p
            })

        } else {
            hValue <- as.numeric(STAT)
            PVAL <- sapply(hValue, function(hh)
                pHS(h = hh, k = k))
        }

        ## to matrix
        P <- matrix(PVAL,
                    nrow = k - 1,
                    ncol = k - 1,
                    byrow = FALSE)
        colnames(P) <- colnames(STAT)
        row.names(P) <- row.names(STAT)

        #DAT <- data.frame(x, g)
        METH <- c("Hayter-Stone's all-pairs test")
        ans <- list(
            statistic = STAT,
            p.value = P,
            data = NULL,
            method = METH,
            data.name = DNAME,
            alternative = alternative,
            dist = "h",
            p.adjust.method = method
        )
        class(ans) <- "PMCMR"
        return(ans)

    } else {
    ## look for k and v = inf
    nrows <- nrow(TabCrit$hayter.h005)
    kk <- as.numeric(colnames(TabCrit$hayter.h005))
    ## check for kk
    if (k > max(kk) | k < min(kk)) stop("No critical values for k = ", k)

    hCrit <- unlist(TabCrit$hayter.h005[nrows, paste0(k)])
    METHOD <- "Pairwise comparisons using Hayter-Stone's all-pairs test"
    parameter = c(k, Inf)
    names(parameter) <- c("k", "df")

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

#' @rdname hsAllPairsTest
#' @aliases hsAllPairsTest.formula
#' @method hsAllPairsTest formula
#' @template one-way-formula
#' @export
hsAllPairsTest.formula <-
    function(formula, data, subset, na.action,
             alternative = c("greater", "less"),
             method = c("look-up", "boot", "asympt"),
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
    y <- do.call("hsAllPairsTest", c(as.list(mf),
                               alternative = alternative,
                               method = method,
                               nperm = nperm))
    y$data.name <- DNAME
    y
}
