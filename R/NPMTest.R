## NPMTest.R
## Part of the R package: PMCMRplus
##
## Copyright (C) 2017-2020 Thorsten Pohlert
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## A copy of the GNU General Public License is available at
## http://www.r-project.org/Licenses/

#' @name NPMTest
#' @title All-Pairs Comparisons for Simply Ordered Mean Ranksums
#' @description
#' Performs Nashimoto and Wright's all-pairs comparison procedure
#' for simply ordered mean ranksums.
#'
#' @template returnOsrt
#'
#' @details
#' The procedure uses the property of a simple order,
#' \eqn{\theta_m' - \theta_m \le \theta_j - \theta_i \le \theta_l' - \theta_l
#' \qquad (l \le i \le m~\mathrm{and}~ m' \le j \le l')}.
#' The null hypothesis H\eqn{_{ij}: \theta_i = \theta_j} is tested against
#' the alternative A\eqn{_{ij}: \theta_i < \theta_j} for any
#' \eqn{1 \le i < j \le k}.
#'
#' The all-pairs comparisons test statistics for a balanced design are
#' \deqn{
#'  \hat{h}_{ij} = \max_{i \le m < m' \le j} \frac{\left(\bar{R}_{m'} - \bar{R}_m \right)}{\sigma_a / \sqrt{n}},
#' }{%
#'  SEE PDF
#' }
#'
#' with \eqn{n = n_i; ~ N = \sum_i^k n_i ~~ (1 \le i \le k)}, \eqn{\bar{R}_i} the mean rank for the \eqn{i}th group,
#' and \eqn{\sigma_a = \sqrt{N \left(N + 1 \right) / 12}}. The null hypothesis is rejected,
#' if \eqn{h_{ij} > h_{k,\alpha,\infty}}.
#'
#' For the unbalanced case with moderate imbalance the test statistic is
#' \deqn{
#'  \hat{h}_{ij} = \max_{i \le m < m' \le j} \frac{\left(\bar{R}_{m'} - \bar{R}_m \right)}
#'  {\sigma_a \left(1/n_m + 1/n_{m'}\right)^{1/2}},
#' }{%
#'  SEE PDF
#' }
#'
#' The null hypothesis is rejected, if \eqn{\hat{h}_{ij} > h_{k,\alpha,\infty} / \sqrt{2}}.
#'
#' If \code{method = "look-up"} the function will not return
#' p-values. Instead the critical h-values
#' as given in the tables of Hayter (1990) for
#' \eqn{\alpha = 0.05} (one-sided)
#' are looked up according to the number of groups (\eqn{k}) and
#' the degree of freedoms (\eqn{v = \infty}).
#'
#' If \code{method = "boot"} an asymetric permutation test
#' is conducted and \eqn{p}-values is returned.
#'
#' @note
#' The function will give a warning for the unbalanced case and returns the
#' critical value \eqn{h_{k,\alpha,\infty} / \sqrt{2}}.
#'
#' @references
#' Hayter, A. J.(1990) A One-Sided Studentised Range
#' Test for Testing Against a Simple Ordered Alternative,
#' \emph{Journal of the American Statistical Association}
#' \bold{85}, 778--785.
#'
#' Nashimoto, K., Wright, F.T. (2007)
#'  Nonparametric Multiple-Comparison Methods for Simply
#'  Ordered Medians.
#'  \emph{Comput Stat Data Anal} \bold{51}, 5068–5076.
#'
#' @keywords htest nonparametric
#'
#' @return
#' Either a list of class \code{"PMCMR"} or a
#' list with class \code{"osrt"} that contains the following
#' components:
#' @template returnOsrt
#' @template class-PMCMR
#'
#' @seealso
#' \code{\link{MTest}}
#'
#' @example examples/shirleyEx.R
#'
#' @export
NPMTest <- function(x, ...)
    UseMethod("NPMTest")

#' @rdname NPMTest
#' @aliases NPMTest.default
#' @method NPMTest default
#' @template one-way-parms
#' @importFrom stats complete.cases
#' @param alternative the alternative hypothesis. Defaults to \code{greater}.
#' @param method a character string specifying the test statistic to use.
#' Defaults to \code{"look-up"} that uses published Table values of Williams (1972).
#' @param nperm number of permutations for the asymptotic permutation test.
#' Defaults to \code{1000}. Ignored, if \code{method = "look-up"}.
#' @export
NPMTest.default <-
    function(x,
             g,
             alternative = c("greater", "less"),
             method = c("look-up", "boot"),
             nperm = 1E4,
             ...) {
        ## taken from stats::kruskal.test

        if (is.list(x)) {
            if (length(x) < 2L)
                stop("'x' must be a list with at least 2 elements")
            DNAME <- deparse(substitute(x))
            x <- lapply(x, function(u)
                u <- u[complete.cases(u)])
            k <- length(x)
            l <- sapply(x, "length")
            if (any(l == 0))
                stop("all groups must contain data")
            g <- factor(rep(1:k, l))
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

        ni <- tapply(x, g, length)
        ## balanced design
        n <- ni[1]
        ## check for all equal
        ok <- sapply(2:k, function(i) ni[i] == n)
        is.balanced <- all(ok)

        hStat <- function(x0, ix, g, is.balanced) {
            ## prepare tukey test
            x <- x0[ix]
            rij <- rank(x)
            ni <- tapply(x, g, length)
            k <- nlevels(g)
            df <- Inf
            Ri <- tapply(rij, g, mean)
            N <- length(x)
            sigma <- sqrt(N * (N + 1) / 12)



            STAT <- matrix(NA, ncol = k - 1, nrow = k - 1)
            if (is.balanced) {
                for (i in 1:(k - 1)) {
                    for (j in (i + 1):k) {
                        u <- j
                        m <- i:(u - 1)
                        tmp <- sapply(m, function(m)
                            (Ri[u] - Ri[m]) /
                                (sigma / sqrt(n)))
                        STAT[j - 1, i] <- max(tmp)
                    }
                }

            } else {
                for (i in 1:(k - 1)) {
                    for (j in (i + 1):k) {
                        u <- j
                        m <- i:(u - 1)
                        tmp <- sapply(m, function(m)
                            (Ri[u] - Ri[m]) /
                                (sigma * sqrt(1 / ni[u] + 1 / ni[m])))
                        STAT[j - 1, i] <- max(tmp)
                    }
                }
            }

            return(STAT)
        }

        ## first cal
        l <- seq_along(x)
        STAT <- hStat(x, l, g, is.balanced)

        k <- nlevels(g)
        colnames(STAT) <- levels(g)[1:(k - 1)]
        rownames(STAT) <- levels(g)[2:k]

        if (method == "boot") {
            ## permutation
            hValue <- as.numeric(STAT)
            m <- length(hValue)
            mt <- matrix(NA, ncol = m, nrow = nperm)
            for (i in 1:nperm) {
                ix <- sample(l)
                tmp <- hStat(x, ix, g, is.balanced )
                mt[i,] <- as.numeric(tmp)
            }

            ## pvalues

            PVAL <- sapply(1:m, function(j) {
                p <- sum(mt[, j] >= hValue[j]) / nperm
                p
            })

            ## to matrix
            P <- matrix(PVAL,
                        nrow = k - 1,
                        ncol = k - 1,
                        byrow = FALSE)
            colnames(P) <- colnames(STAT)
            row.names(P) <- row.names(STAT)

            #DAT <- data.frame(x, g)
            METH <- c("Nashimoto-Wright's NPM-Test")
            ans <- list(
                statistic = STAT,
                p.value = P,
                data = NULL,
                method = METH,
                data.name = DNAME,
                alternative = alternative,
                dist = "h",
                p.adjust.method = "boot"
            )
            class(ans) <- "PMCMR"
            return(ans)

        } else {
            ## get critical h-value with k = k and v = Inf
            nrows <- nrow(TabCrit$hayter.h005)
            kk <- as.numeric(colnames(TabCrit$hayter.h005))
            ## check for kk
            if (k > max(kk) |
                k < min(kk))
                stop("No critical values for k = ", k)

            hCrit <- unlist(TabCrit$hayter.h005[nrows, paste0(k)])
            hCrit <- adjust.hCrit(hCrit, is.balanced)

            METHOD <-
                "Pairwise comparisons using Nashimoto-Wright's NPM-Test"
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

#' @rdname NPMTest
#' @aliases NPMTest.formula
#' @method NPMTest formula
#' @template one-way-formula
#' @export
NPMTest.formula <-
    function(formula,
             data,
             subset,
             na.action,
             alternative = c("greater", "less"),
             method = c("look-up", "boot"),
             nperm = 1E4,
             ...)
    {
        mf <- match.call(expand.dots = FALSE)
        m <-
            match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
        mf <- mf[c(1L, m)]
        mf[[1L]] <- quote(stats::model.frame)

        if (missing(formula) || (length(formula) != 3L))
            stop("'formula' missing or incorrect")
        mf <- eval(mf, parent.frame())
        if (length(mf) > 2L)
            stop("'formula' should be of the form response ~ group")
        DNAME <- paste(names(mf), collapse = " by ")
        names(mf) <- NULL
        alternative <- match.arg(alternative)
        method <- match.arg(method)
        y <- do.call("NPMTest",
                     c(
                         as.list(mf),
                         alternative = alternative,
                         method = method,
                         nperm = nperm
                     ))
        y$data.name <- DNAME
        y
    }
