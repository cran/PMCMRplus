##  shanTest.R
##
##  Copyright (C) 2020 Thorsten Pohlert
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
##
#' @name shanTest
#' @title Testing against Ordered Alternatives (Shan-Young-Kang Test)
#'
#' @description
#' Performs the Shan-Young-Kang test for testing against ordered alternatives.
#' @details
#' The null hypothesis, H\eqn{_0: \theta_1 = \theta_2 = \ldots = \theta_k}
#' is tested against a simple order hypothesis,
#' H\eqn{_\mathrm{A}: \theta_1 \le \theta_2 \le \ldots \le
#' \theta_k,~\theta_1 < \theta_k}.
#'
#' Let \eqn{R_{ij}} be the rank of \eqn{X_{ij}},
#' where \eqn{X_{ij}} is jointly ranked
#' from \eqn{\left\{1, 2, \ldots, N \right\}, ~~ N = \sum_{i=1}^k n_i},
#' the the test statistic is
#'
#' \deqn{
#'  S = \sum_{i = 1}^{k-1} \sum_{j = i + 1}^k D_{ij},
#' }{%
#'  SEE PDF
#' }
#'
#' with
#' \deqn{
#'  D_{ij} = \sum_{l = 1}^{n_i} \sum_{m=1}^{n_j}  \left(R_{jm} - R_{il} \right)~ \mathrm{I}\left(X_{jm} > X_{il} \right),
#' }{%
#'  SEE PDF
#' }
#'
#' where
#'
#' \deqn{
#'  \mathrm{I}(u) = \left\{
#'  \begin{array}{c}
#'      1, \qquad \forall~ u > 0 \\
#'      0, \qquad \forall~ u \le 0
#'   \end{array}
#'   \right.
#' .}{%
#'  SEE PDF.
#' }
#'
#' The test statistic is asymptotically normal distributed:
#' \deqn{
#'  z = \frac{S - \mu_{\mathrm{S}}}{\sqrt{s^2_{\mathrm{S}}}}
#' }{
#'  SEE PDF
#' }
#'
#' The p-values are estimated from the standard normal distribution.
#'
#' @note
#' The variance estimation (see Theorem 2.1, Shan et al. 2014)
#' can become negative for certain combinations of \eqn{N,~n_i,~k
#' \qquad (1 \le i \le k)}. In these cases the function will return
#' a warning and the returned p-value will be \code{NaN}.
#'
#' @references
#' Shan, G., Young, D., Kang, L. (2014) A New Powerful Nonparametric
#' Rank Test for Ordered Alternative Problem. PLOS ONE 9, e112924.
#' https://doi.org/10.1371/journal.pone.0112924
#'
#' @template class-htest
#' @template trendTests
#' @export shanTest
shanTest <- function(x, ...)
    UseMethod("shanTest")

#' @rdname shanTest
#' @method shanTest default
#' @aliases shanTest.default
#' @template one-way-parms
#' @param alternative the alternative hypothesis.
#' Defaults to \code{"greater"}.
#' @importFrom stats pnorm complete.cases
#' @export
shanTest.default <-
    function(x,
             g,
             alternative = c("greater", "less"),
             ...)
    {
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
            ## check incoming from formula
            if (is.null(x$alternative)) {
                alternative <- "greater"
            } else {
                alternative <- x$alternative
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

        if (k < 3) {
            stop("a minimum of k = 3 groups is required.")
        }
        alternative <- match.arg(alternative)
        N <- length(x)
        if (N < 2)
            stop("not enough observations")
        nij <- tapply(x, g, length)


        if (alternative == "less") {
            x <- -x
        }

        ## Indicator function
        Igt <- function(x, y) {
            ifelse(x > y, 1, 0)
        }

        ### START HERE FOR PERMUTATION TEST
        ## rank data
        R <- rank(x)

        ## create and populate matrix
        X <- matrix(NA, ncol = k, nrow = max(nij))
        j <- 0
        for (i in 1:k) {
            for (l in 1:nij[i]) {
                j = j + 1
                X[l, i] <- R[j]
            }
        }

        ## Intermediate Dij
        Dij <- function(i, j, X) {
            ni <- nij[i]
            nj <- nij[j]
            Ril <- R[g == levels(g)[i]]
            Rjm <- R[g == levels(g)[j]]
            sumDij <- 0
            for (l in (1:ni)) {
                for (m in (1:nj)) {
                    sumDij <- sumDij +
                        (Rjm[m] - Ril[l]) * Igt(X[m, j], X[l, i])
                }
            }
            sumDij
        }

        ## Test statistic
        S <- 0
        for (i in (1:(k - 1))) {
            for (j in ((i + 1):k)) {
                S = S + Dij(i, j, X)
            }
        }

        ## mean
        Sum1 <- 0
        for (i in 1:(k - 1)) {
            for (j in (i + 1):k) {
                Sum1 <- Sum1 + nij[i] * nij[j]
            }
        }
        ES <- Sum1 * (N + 1) / 6

        ## variance
        CovA <- (2 * N ^ 2 + N - 1) / 90
        CovB <- (-7 * N ^ 2 - 11 * N - 4) / 360

        term1 <- ((N ^ 2 + N) / 12 - (N + 1) ^ 2 / 36) * Sum1

        Sum3 <- 0
        for (i in 1:(k - 2)) {
            for (j in (i + 1):(k - 1)) {
                for (l in (j + 1):k) {
                    Sum3 <- Sum3 + nij[i] * nij[j] * nij[l]
                }
            }
        }
        term3 <- 2 * Sum3 * CovB

        Sum1 <- 0
        for (i in 1:(k - 1)) {
            for (j in (i + 1):k) {
                Sum1 <- Sum1 + nij[i] * choose(nij[j], 2)
            }
        }

        Sum2 <- 0
        for (i in 2:k) {
            for (j in 1:(i - 1)) {
                Sum2 <- Sum2 + nij[i] * choose(nij[j], 2)
            }
        }

        term2 <- 2 * (Sum1 + Sum2) * CovA

        VarS <- term1 + term2 + term3

        ## z-statistic
        STATISTIC <- (S - ES) / sqrt(VarS)
        ### END HERE FOR PERMUTATION TEST ###

        if (alternative == "less") {
            STATISTIC <- -STATISTIC
            PVAL <- pnorm(STATISTIC, lower.tail = TRUE)
        } else {
            PVAL <- pnorm(STATISTIC, lower.tail = FALSE)
        }

        ## check for ties
        TIES <- (sum(table(x) - 1) > 0)
        if (TIES) {
            warning("Ties are present. No correction for ties.")
        }

        ESTIMATES <- S
        names(ESTIMATES) <- "S"
        names(STATISTIC) <- "z"
        RVAL <- list(
            statistic = STATISTIC,
            p.value = PVAL,
            method = "Shan-Young-Kang test",
            data.name = DNAME,
            alternative = alternative,
            estimates = ESTIMATES
        )
        class(RVAL) <- "htest"
        return(RVAL)
    }

#' @rdname shanTest
#' @method shanTest formula
#' @aliases shanTest.formula
#' @template one-way-formula
#' @export
shanTest.formula <-
    function(formula,
             data,
             subset,
             na.action,
             alternative = c("greater", "less"),
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
        alternative <- match.arg(alternative)
        names(mf) <- NULL
        y <- do.call("shanTest",
                     c(as.list(mf), alternative = alternative))
        y$data.name <- DNAME
        y
    }
