## tamhaneT2Test.R
## Part of the R package: PMCMRplus
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
#' the T2 test (or T2' test) of Tamhane can be performed. A total of \eqn{m = k(k-1)/2}
#' hypotheses can be tested. The null hypothesis
#' H\eqn{_{ij}: \mu_i(x) = \mu_j(x)} is tested in the two-tailed test
#' against the alternative
#' A\eqn{_{ij}: \mu_i(x) \ne \mu_j(x), ~~ i \ne j}.
#'
#' T2 test uses Welch's approximate solution for
#' calculating the degree of freedom. T2' test uses the usual
#' \eqn{df = N - 2} approximation. A warning message appears
#' in the modified T2' test, if none of in Tamhane (1979) given conditions for nearly balanced
#' sample sizes and nearly balanced standard errors is true.
#'
#' The p-values are computed from the t-distribution and adjusted
#' according to Dunn-Sidak.
#'
#' @note
#' T2 test is basically an all-pairs pairwise-t-test. Similar results
#' can be obtained with \code{pairwise.t.test(..., var.equal=FALSE, p.adjust.mehod = FALSE)}.
#'
#' Thanks to Sirio Bolaños for his kind suggestion for adding T2' test
#' into this function.
#' @references
#'  Tamhane, A. C. (1979) A Comparison of Procedures for Multiple Comparisons
#'  of Means with Unequal Variances, \emph{Journal of the American
#'  Statistical Association} \bold{74}, 471--480.
#'
#' @keywords htest
#' @concept allPairsComparisons
#'
#' @examples
#' set.seed(245)
#' mn <- rep(c(1, 2^(1:4)), each=5)
#' sd <- rep(1:5, each=5)
#' x <- mn + rnorm(25, sd = sd)
#' g <- factor(rep(1:5, each=5))
#'
#' fit <- aov(x ~ g)
#' shapiro.test(residuals(fit))
#' bartlett.test(x ~ g) # var1 != varN
#' anova(fit)
#' summary(T2 <- tamhaneT2Test(x, g))
#' T2
#' ## compare with pairwise.t.test
#' WT <- pairwise.t.test(x, g, pool.sd = FALSE, p.adjust.method = "none")
#' p.adj.sidak <- function(p, m) sapply(p, function(p) min(1, 1 - (1 - p)^m))
#' p.raw <- as.vector(WT$p.value)
#' m <- length(p.raw[!is.na(p.raw)])
#' PADJ <- matrix(ans <- p.adj.sidak(p.raw, m),
#'                nrow = 4, ncol = 4)
#' colnames(PADJ) <- colnames(WT$p.value)
#' rownames(PADJ) <- rownames(WT$p.value)
#' PADJ
#'
#' ## same without Welch's approximate solution
#' summary(T2b <- tamhaneT2Test(x, g, welch = FALSE))
#'
#'
#' @seealso
#' \code{\link{dunnettT3Test}}
#' @importFrom stats pt
#' @importFrom stats complete.cases
#' @importFrom stats var
#' @importFrom stats pairwise.table
#' @export
tamhaneT2Test <- function(x, ...) UseMethod("tamhaneT2Test")

#' @rdname tamhaneT2Test
#' @method tamhaneT2Test default
#' @aliases tamhaneT2Test.default
#' @template one-way-parms
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
