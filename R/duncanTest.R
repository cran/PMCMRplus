## duncanTest.R
## Part of the R package: PMCMRplus
##
## Copyright (C) 2018-2020 Thorsten Pohlert
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

## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## ToDo Check p.values as they differ from duncan.test
## in pkg agricolae. Done 2018-07-04

#' @name duncanTest
#' @title Duncan's Multiple Range Test
#' @description
#' Performs Duncan's all-pairs comparisons test for normally distributed
#' data with equal group variances.
#' @details
#' For all-pairs comparisons in an one-factorial layout
#' with normally distributed residuals and equal variances
#' Duncan's multiple range test can be performed. A total of \eqn{m = k(k-1)/2}
#' hypotheses can be tested. The null hypothesis
#' H\eqn{_{ij}: \mu_i(x) = \mu_j(x)} is tested in the two-tailed test
#' against the alternative
#' A\eqn{_{ij}: \mu_i(x) \ne \mu_j(x), ~~ i \ne j}.
#'
#' The p-values are computed from the Tukey-distribution.
#'
#' @template class-PMCMR
#'
#' @references
#' Duncan, D. B. (1955) Multiple range and multiple F tests,
#' \emph{Biometrics} \bold{11}, 1--42.
#'
#' @keywords htest
#' @concept parametric
#' @seealso
#' \code{\link[stats]{Tukey}}, \code{\link[stats]{TukeyHSD}} \code{\link{tukeyTest}}
#' @examples
#' fit <- aov(weight ~ feed, chickwts)
#' shapiro.test(residuals(fit))
#' bartlett.test(weight ~ feed, chickwts)
#' anova(fit)
#'
#' ## also works with fitted objects of class aov
#' res <- duncanTest(fit)
#' summary(res)
#' summaryGroup(res)
#' @export
duncanTest <- function(x, ...) UseMethod("duncanTest")

#' @rdname duncanTest
#' @aliases duncanTest.default
#' @method duncanTest default
#' @template one-way-parms-aov
#' @importFrom stats complete.cases
#' @importFrom stats var
#' @importFrom stats ptukey
#' @export
duncanTest.default <- function(x, g, ...){
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

    ## prepare snk-test
    ni <- tapply(x, g, length)
    n <- sum(ni)
    xi <- tapply(x, g, mean)
    s2i <- tapply(x, g, var)
    s2in <- 1 / (n - k) * sum(s2i * (ni - 1))
    df <- n - k

    ## harmonic mean
    r = k / sum(1/ni)

    ## order means
    o <- order(xi, decreasing = TRUE)
    Xo <- xi[o]
    oName <- names(Xo)

    levNames <- levels(g)
    qval <- matrix(NA, ncol = k, nrow = k)
    colnames(qval) <- levNames
    rownames(qval) <- levNames
    pval <- qval

    ## this is sorted
    for (j in 1:(k-1)) {
      for (i in (j+1):k){

        ## Statistic
        T <- (Xo[j] - Xo[i]) / sqrt(s2in / r)

        ## range
        p <- 1 + abs(j - i)

        ## p-Value
        pp <- ptukey(q = abs(T),
                     nmeans = p,
                     df = df,
                     lower.tail = FALSE)

        ## bonferroni adjustment
        pp <- min(1, 1 - (1 - pp)^(1 / (p - 1)))

        ## assign
        ii <- oName[i]
        jj <- oName[j]
        qval[ii, jj] <- T
        qval[jj, ii] <- T
        pval[ii, jj] <- pp
        pval[jj, ii] <- pp
      }
    }

    pval[upper.tri(pval)] <- NA
    qval[upper.tri(pval)] <- NA

    MODEL <- data.frame(x, g)
    DIST <- "q"
    METHOD <- "Duncan's multiple range test"
    ans <- list(method = METHOD, data.name = DNAME,
                p.value = pval[2:k, 1:(k-1)],
                statistic = qval[2:k, 1:(k-1)],
                p.adjust.method = "duncan",
                model = MODEL, dist = DIST, alternative = "two.sided")
    class(ans) <- "PMCMR"
    ans
}

#' @rdname duncanTest
#' @aliases duncanTest.formula
#' @method duncanTest formula
#' @template one-way-formula
#' @export
duncanTest.formula <-
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
    y <- do.call("duncanTest", c(as.list(mf)))
    y$data.name <- DNAME
    y
}

#' @rdname duncanTest
#' @aliases duncanTest.aov
#' @method duncanTest aov
## @param x A fitted model object, usually an \link[stats]{aov} fit.
#' @export
duncanTest.aov <- function(x, ...) {
  model <- x$model
  DNAME <- paste(names(model), collapse = " by ")
  names(model) <- c("x", "g")
  y <- do.call("duncanTest", as.list(model))
  y$data.name <- DNAME
  y
}
