## tamhaneT2Test.R
## Part of the R package: PMCMRplus
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
#' @name tamhaneT2Test
#' @title Tamhane's T2 Test
#' @description
#' Performs Tamhane's T2 (or T2') all-pairs comparison test for normally distributed
#' data with unequal variances.
#'
#' @template class-PMCMR
#'
#' @details
#' For all-pairs comparisons in an one-factorial layout
#' with normally distributed residuals but unequal groups variances
#' the T2 test (or T2' test) of Tamhane can be performed.
#' Let \eqn{X_{ij}} denote a continuous random variable
#' with the \eqn{j}-the realization (\eqn{1 \le j \le n_i})
#' in the \eqn{i}-th group (\eqn{1 \le i \le k}). Furthermore, the total
#' sample size is \eqn{N = \sum_{i=1}^k n_i}. A total of \eqn{m = k(k-1)/2}
#' hypotheses can be tested: The null hypothesis is
#' H\eqn{_{ij}: \mu_i = \mu_j ~~ (i \ne j)} is tested against the alternative
#' A\eqn{_{ij}: \mu_i \ne \mu_j} (two-tailed). Tamhane T2 all-pairs
#' test statistics are given by
#'
#' \deqn{
#'  t_{ij} \frac{\bar{X}_i - \bar{X_j}}
#'  {\left( s^2_j / n_j + s^2_i / n_i \right)^{1/2}}, ~~
#'  (i \ne j)
#' }{%
#'  SEE PDF
#' }
#'
#' with \eqn{s^2_i} the variance of the \eqn{i}-th group.
#' The null hypothesis is rejected (two-tailed) if
#'
#' \deqn{
#'  \mathrm{Pr} \left\{ |t_{ij}| \ge t_{v_{ij}\alpha'/2} | \mathrm{H} \right\}_{ij} =
#'  \alpha.
#' }{%
#'  SEE PDF
#' }
#'
#' T2 test uses Welch's approximate solution for
#' calculating the degree of freedom.
#'
#' \deqn{
#'  v_{ij} = \frac{\left( s^2_i / n_i + s^2_j / n_j \right)^2}
#'  {s^4_i / n^2_i \left(n_i - 1\right) + s^4_j / n^2_j \left(n_j - 1\right)}.
#' }{%
#'  SEE PDF
#' }
#'
#' T2' test applies the following approximation for the degree of freedom
#' \deqn{
#'  v_{ij} = n_i + n_j - 2
#' }{%
#'  SEE PDF
#' }
#'
#' The p-values are computed from the \code{\link[stats]{TDist}}-distribution
#' and adjusted according to Dunn-Sidak.
#' \deqn{
#'  p'_{ij} = \min \left\{1, ~ (1 - (1 - p_{ij})^m)\right\}
#' }{%
#'  SEE PDF
#' }
#'
#' @note
#' T2 test is basically an all-pairs pairwise-t-test. Similar results
#' can be obtained with \code{pairwise.t.test(..., var.equal=FALSE, p.adjust.mehod = FALSE)}.
#'
#' A warning message appears
#' in the modified T2' test, if none of in Tamhane (1979) given conditions
#' for nearly balanced
#' sample sizes and nearly balanced standard errors is true.
#'
#' Thanks to Sirio Bolaños for his kind suggestion for adding T2' test
#' into this function.
#'
#' @references
#'  Tamhane, A. C. (1979) A Comparison of Procedures for Multiple Comparisons
#'  of Means with Unequal Variances, \emph{Journal of the American
#'  Statistical Association} \bold{74}, 471--480.
#'
#' @examples
#' fit <- aov(weight ~ feed, chickwts)
#' shapiro.test(residuals(fit))
#' bartlett.test(weight ~ feed, chickwts) # var1 = varN
#' anova(fit)
#'
#' ## also works with fitted objects of class aov
#' res <- tamhaneT2Test(fit)
#' summary(res)
#' summaryGroup(res)
#' res
#'
#' ## compare with pairwise.t.test
#' WT <- pairwise.t.test(chickwts$weight,
#'                       chickwts$feed,
#'                       pool.sd = FALSE,
#'                       p.adjust.method = "none")
#' p.adj.sidak <- function(p, m) sapply(p, function(p) min(1, 1 - (1 - p)^m))
#' p.raw <- as.vector(WT$p.value)
#' m <- length(p.raw[!is.na(p.raw)])
#' PADJ <- matrix(ans <- p.adj.sidak(p.raw, m),
#'                nrow = 5, ncol = 5)
#' colnames(PADJ) <- colnames(WT$p.value)
#' rownames(PADJ) <- rownames(WT$p.value)
#' PADJ
#'
#' ## same without Welch's approximate solution
#' summary(T2b <- tamhaneT2Test(fit, welch = FALSE))
#'
#' @seealso
#' \code{\link{dunnettT3Test}} \code{\link{uryWigginsHochbergTest}}
#' @importFrom stats pt
#' @importFrom stats complete.cases
#' @importFrom stats var
#' @importFrom stats pairwise.table
#' @export
tamhaneT2Test <- function(x, ...) UseMethod("tamhaneT2Test")

#' @rdname tamhaneT2Test
#' @method tamhaneT2Test default
#' @aliases tamhaneT2Test.default
#' @template one-way-parms-aov
#' @param welch indicates, whether Welch's approximate solution for
#' calculating the degree of freedom shall be used or, as usually,
#' \eqn{df = N - 2}. Defaults to \code{TRUE}.
#' @export
tamhaneT2Test.default <-
function(x, g, welch = TRUE, ...){
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

    ## prepare tamhaneT2 test
    ni <- tapply(x, g, length)
    n <- sum(ni)
    xi <- tapply(x, g, mean)
    s2i <- tapply(x, g, var)
    s2in <- 1 / (n - k) * sum(s2i * (ni - 1))

    compare.stats <- function(i,j) {
        dif <- xi[i] - xi[j]
        A <- (s2i[i] / ni[i] + s2i[j] / ni[j])
        tval <- dif / sqrt(A)
        return(tval)
    }

    PSTAT <- pairwise.table(compare.stats,levels(g),
                            p.adjust.method="none" )

    compare.levels <- function(i,j) {
        dif <- xi[i] - xi[j]
        A <- (s2i[i] / ni[i] + s2i[j] / ni[j])
        tval <- dif / sqrt(A)

        if (welch) {
        df <- A^2 / (s2i[i]^2 / (ni[i]^2 * (ni[i] - 1)) +
                    s2i[j]^2 / (ni[j]^2 * (ni[j] - 1)))
        } else {
          ## checks according to Tamhane (1979, p. 474)
          ok1 <- 9/10 <= ni[i]/ni[j] & ni[i]/ni[j] <= 10/9
          ok2 <- 9/10 <= (s2i[i] / ni[i]) / (s2i[j] / ni[j]) &
            (s2i[i] / ni[i]) / (s2i[j] / ni[j]) <= 10/9
          ok3 <- 4/5 <= ni[i]/ni[j] & ni[i]/ni[j] <= 5/4 &
            1/2 <= (s2i[i] / ni[i]) / (s2i[j] / ni[j]) &
            (s2i[i] / ni[i]) / (s2i[j] / ni[j]) <= 2
          ok4 <- 2/3 <= ni[i]/ni[j] & ni[i]/ni[j] <= 3/2 &
            3/4 <= (s2i[i] / ni[i]) / (s2i[j] / ni[j]) &
            (s2i[i] / ni[i]) / (s2i[j] / ni[j]) <= 4/3
          OK <- any(ok1, ok2, ok3, ok4)
          if (!OK) {
            warning("Sample sizes or standard errors are not balanced. T2 test is recommended.")
          }
          df = ni[i] + ni[j] - 2
        }

        ## for two-sided test, it should be 2 times pt
        pval <- 2 * pt(abs(tval), df = df, lower.tail = FALSE)
        return(pval)
    }

    METHOD <- ifelse(welch,
                     "Tamhane's T2-test for unequal variances",
                     "Tamhane's T2'-test for unequal variances")
    PVAL <-  pairwise.table(compare.levels,
                            level.names = levels(g),
                            p.adjust.method = "none")
    pval <- as.vector(PVAL)

    m <- k * (k - 1) / 2
    padj <- sapply(pval, function(p) {
       min(1, (1 - (1 - p)^m))
    })
    PVAL <- matrix(padj, ncol=(k-1), nrow=(k-1))
    colnames(PVAL) <- colnames(PSTAT)
    rownames(PVAL) <- rownames(PSTAT)
    MODEL <- data.frame(x, g)
    DIST <- "t"
    ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL,
                statistic = PSTAT, p.adjust.method = ifelse(welch, "T2 (Sidak)", "T2' (Sidak)"),
                model = MODEL, dist = DIST, alternative = "two.sided")
    class(ans) <- "PMCMR"
    ans
}

#' @rdname tamhaneT2Test
#' @method tamhaneT2Test formula
#' @aliases tamhaneT2Test.formula
#' @template one-way-formula
#' @export
tamhaneT2Test.formula <-
function(formula, data, subset, na.action, welch = TRUE, ...)
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
    y <- do.call("tamhaneT2Test", c(as.list(mf), welch = welch))
    y$data.name <- DNAME
    y
}

#' @rdname tamhaneT2Test
#' @aliases tamhaneT2Test.aov
#' @method tamhaneT2Test aov
# @param obj A fitted model object, usually an \link[stats]{aov} fit.
#' @export
tamhaneT2Test.aov <- function(x, welch = TRUE, ...) {
    model <- x$model
    DNAME <- paste(names(model), collapse = " by ")
    names(model) <- c("x", "g")
    parms <- c(as.list(model), list(welch = welch))
    y <- do.call("tamhaneT2Test", parms)
    y$data.name <- DNAME
    y
}
