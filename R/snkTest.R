## snkTest.R
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

#' @name snkTest
#' @title Student-Newman-Keuls Test
#' @description
#' Performs Student-Newman-Keuls all-pairs comparisons test for normally distributed
#' data with equal group variances.
#' @details
#' For all-pairs comparisons in an one-factorial layout
#' with normally distributed residuals and equal variances
#' Student-Newman-Keuls test can be performed. A total of \eqn{m = k(k-1)/2}
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
#' Keuls, M. (1952) The use of the "studentized range"
#' in connection with an analysis of variance,
#' \emph{Euphytica} \bold{1}, 112--122.
#'
#' Newman, D. (1939) The distribution of range in
#' samples from a normal population, expressed in
#' terms of an independent estimate of standard
#' deviation, \emph{Biometrika} \bold{31}, 20--30.
#'
#' Student (1927) Errors of routine analysis,
#' \emph{Biometrika} \bold{19}, 151--164.
#'
#' @concept AllPairsComparison
#' @keywords htest
#' @seealso
#' \code{\link[stats]{Tukey}}, \code{\link[stats]{TukeyHSD}} \code{\link{tukeyTest}}
#' @examples
#' set.seed(245)
#' mn <- rep(c(1, 2^(1:4)), each=5)
#' sd <- rep(1, 25)
#' x <- mn + rnorm(25, sd = sd)
#' g <- factor(rep(1:5, each=5))
#'
#' fit <- aov(x ~ g)
#' shapiro.test(residuals(fit))
#' bartlett.test(x ~ g) # var1 = varN
#' anova(fit)
#' summary(snkTest(x, g))
#' @export
snkTest <- function(x, ...) UseMethod("snkTest")

#' @rdname snkTest
#' @aliases snkTest.default
#' @method snkTest default
#' @template one-way-parms
#' @importFrom stats complete.cases
#' @importFrom stats var
#' @importFrom stats ptukey
#' @export
snkTest.default <- function(x, g, ...){
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
    df <- n -k

    ## order means
    o <- order(xi, decreasing = TRUE)
    levNames <- levels(g)
    qval <- matrix(NA, ncol = k-1, nrow = k-1)
    colnames(qval) <- levNames[1:(k - 1)]
    rownames(qval) <- levNames[2:k]
    pval <- qval

    for (j in 1:(k-1)) {
      for (i in (j+1):k){
        qval[i-1, j] <- (xi[i] - xi[j]) /
          sqrt((s2in / 2) * (1 / ni[i] + 1 / ni[j]))

        p <- 1 + abs(which(o == i, arr.ind = TRUE) -
                       which(o == j, arr.ind = TRUE))

        pval[i-1, j] <- ptukey(q = abs(qval[i-1, j]),
                              nmeans = p,
                              df = df,
                              lower.tail = FALSE)
      }
    }

    MODEL <- data.frame(x, g)
    DIST <- "q"
    METHOD <- "SNK test"
    ans <- list(method = METHOD, data.name = DNAME, p.value = pval,
                statistic = qval, p.adjust.method = "step down",
                model = MODEL, dist = DIST, alternative = "two.sided")
    class(ans) <- "PMCMR"
    ans
}

#' @rdname snkTest
#' @aliases snkTest.formula
#' @method snkTest formula
#' @template one-way-formula
#' @export
snkTest.formula <-
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
    y <- do.call("snkTest", c(as.list(mf)))
    y$data.name <- DNAME
    y
}
