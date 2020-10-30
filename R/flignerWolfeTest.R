# flignerWolfeTest.R
# Part of the R package: PMCMR
#
# Copyright (C) 2020 Thorsten Pohlert
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
#' @name flignerWolfeTest
#' @title Testing Several Treatments With One Control
#' @description
#' Performs Fligner-Wolfe non-parametric test for
#' simultaneous testing of several locations of treatment groups
#' against the location of the control group.
#'
#' @details
#' For a one-factorial layout with non-normally distributed residuals
#' the Fligner-Wolfe test can be used.
#'
#' Let there be \eqn{k-1}-treatment groups and one control group, then
#' the null hypothesis, H\eqn{_0: \theta_i - \theta_c = 0 ~ (1 \le i \le k-1)}
#' is tested against the alternative (greater),
#' A\eqn{_1: \theta_i - \theta_c > 0 ~ (1 \le i \le k-1)},
#' with at least one inequality being strict.
#'
#' Let \eqn{n_c} denote the sample size of the control group,
#' \eqn{N^t = \sum_{i=1}^{k-1} n_i} the sum of all treatment
#' sample sizes and \eqn{N = N^t + n_c}. The test statistic without taken
#' ties into account is
#'
#' \deqn{
#'  W = \sum_{j=1}^{k-1} \sum_{i=1}^{n_i} r_{ij} -
#'  \frac{N^t \left(N^t + 1 \right) }{2}
#' }{%
#'  SEE PDF
#' }
#'
#' with \eqn{r_{ij}} the rank of variable \eqn{x_{ij}}.
#' The null hypothesis is rejected,
#' if \eqn{W > W_{\alpha,m,n}} with
#' \eqn{m = N^t} and \eqn{n = n_c}.
#'
#' In the presence of ties, the statistic is
#'
#'  \deqn{
#'     \hat{z} = \frac{W - n_c N^t / 2}{s_W},
#'  }{%
#'   SEE PDF
#'  }
#'
#' where
#' \deqn{
#'   \frac{n_c N^t}{12 N \left(N - 1 \right)}
#'   \sum_{j=1}^g t_j \left(t_j^2 - 1\right),
#' }{%
#'  SEE PDF
#' }
#'
#' with \eqn{g} the number of tied groups and \eqn{t_j}
#' the number of tied values in the \eqn{j}th group. The null hypothesis
#' is rejected, if \eqn{\hat{z} > z_\alpha} (as cited in EPA 2006).
#'
#' If \code{dist = Wilcoxon}, then the \eqn{p}-values are estimated from the  \code{\link[stats]{Wilcoxon}}
#' distribution, else the \code{\link[stats]{Normal}} distribution is used. The latter can be used,
#' if ties are present.
#'
#' @template class-htest
#' @template trendTests
#'
#' @references
#' EPA (2006) \emph{Data Quality Assessment:
#' Statistical Methods for Practitioners}
#' (Guideline No. EPA QA/G-9S), US-EPA.
#'
#' Fligner, M.A., Wolfe, D.A. (1982)
#'  Distribution-free tests for comparing several
#'  treatments with a control. \emph{Stat Neerl} \bold{36},
#'  119--127.
#'
#' @concept wilcoxonranks
#' @importFrom stats pwilcox complete.cases
#' @export
flignerWolfeTest <- function(x, ...) UseMethod("flignerWolfeTest")

#' @rdname flignerWolfeTest
#' @aliases flignerWolfeTest.default
#' @method flignerWolfeTest default
#' @template one-way-parms
#' @param alternative the alternative hypothesis. Defaults to \code{"greater"}.
#' @param dist the test distribution. Defaults to \code{"Wilcoxon"}.
#' @export
flignerWolfeTest.default <-
    function(x, g, alternative = c("greater", "less"),
             dist = c("Wilcoxon", "Normal"), ...)
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

    alternative <- match.arg(alternative)
    dist <- match.arg(dist)

    ## rank
    rij <- rank(x)
    ni <- tapply(rij, g, length)
    rSj <- tapply(rij, g, sum)
    Ns <- sum(ni[-1])
    nc <- ni[1]

    ## Wilcoxon type statistic
    ## leave out control group
    FW0 <- sum(rSj[-1]) - (Ns * (Ns + 1)) / 2

    ## check for ties

    if (dist == "Wilcoxon") {
#        TIES <- FALSE
        TIES <- sum(table(rij) - 1) != 0
        if (TIES) {
            warning("Cannot compute exact p-value with ties. You may use the normal distribution to account for ties.")
        }

        PVAL <- switch(
            alternative,
            greater = pwilcox(
                q = FW0,
                n = nc,
                m = Ns,
                lower.tail = FALSE
            ),
            less = pwilcox(FW0, nc, Ns)
        )

        statistic <- c("W" = FW0)
        parameter <- c(nc, Ns)
        names(parameter) <- c("n", "m")


    } else {


        ## are ties present
        getties <- function(x){
            n <- length(x)
            t <- table(x)
            C <- sum(t * (t^2 - 1))
            return(C)
        }

        C <- getties(rij)
        if (C != 0) warning("Ties are present. z-quantiles were corrected for ties.")
        ## compute variance
        N <- nc + Ns
        Var0 <- (nc * Ns * (N + 1)) / 12 -
                (nc * Ns / (12 * N * (N - 1))) * C

        z <- (FW0 - nc * Ns / 2) / sqrt(Var0)

        PVAL <- switch(alternative,
                       greater = pnorm(z, lower.tail = FALSE),
                       less = pnorm(z))
        statistic <- z
        names(statistic) <- "z"
        parameter <- NULL
    }

    ## prepare output
    METHOD <- paste("Fligner-Wolfe test")

    ans <- list(method = METHOD, data.name = DNAME,
                p.value = PVAL,
                statistic = statistic,
                parameter = parameter,
                alternative = alternative,
                null.value = c("shift of treatment location(s) vs control location" = 0))

    class(ans) <- "htest"
    ans
}

#' @rdname flignerWolfeTest
#' @method flignerWolfeTest formula
#' @aliases flignerWolfeTest.formula
#' @template one-way-formula
#' @export
flignerWolfeTest.formula <-
    function(formula, data, subset, na.action,
             alternative = c("greater", "less"),
             dist = c("Wilcoxon", "Normal"), ...)
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
    dist <- match.arg(dist)
    names(mf) <- NULL
    y <- do.call("flignerWolfeTest", c(as.list(mf),
                             alternative = alternative,
                             dist = dist))
    y$data.name <- DNAME
    y
}
