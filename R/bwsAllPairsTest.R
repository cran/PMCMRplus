##  bwsAllPairsTest.R
##
##  Copyright (C) 2017-2020 Thorsten Pohlert
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


#' @title BWS All-Pairs Comparison Test
#' @description Performs Baumgartner-Weiß-Schindler all-pairs comparison test.
#' @details
#' For all-pairs comparisons in an one-factorial layout
#' with non-normally distributed residuals Baumgartner-Weiß-Schindler
#' all-pairs comparison test can be used. A total of \eqn{m = k(k-1)/2}
#' hypotheses can be tested. The null hypothesis
#' H\eqn{_{ij}: F_i(x) = F_j(x)} is tested in the two-tailed test
#' against the alternative
#' A\eqn{_{ij}: F_i(x) \ne F_j(x), ~~ i \ne j}.
#'
#' This function is a wrapper function that sequentially
#' calls \code{\link[BWStest]{bws_test}} for each pair.
#' The default test method (\code{"BWS"}) is the original
#' Baumgartner-Weiß-Schindler test statistic B. For
#' \code{method == "Murakami"} it is the modified BWS statistic
#' denoted B*. The calculated p-values for \code{Pr(>|B|)}
#' or \code{Pr(>|B*|)} can be adjusted to account for Type I error
#' inflation using any method as implemented in \code{\link{p.adjust}}.
#'
#' @name bwsAllPairsTest
#' @template class-PMCMR
#' @keywords htest nonparametric
#' @concept AllPairsComparison
#' @references
#' Baumgartner, W., Weiss, P., Schindler, H. (1998) A nonparametric test for the
#' general two-sample problem, \emph{Biometrics} \bold{54}, 1129--1135.
#'
#' Murakami, H. (2006) K-sample rank test based on modified Baumgartner statistic and its power
#' comparison, \emph{J. Jpn. Comp. Statist.} \bold{19}, 1--13.
#' @examples
#'
#' out <- bwsAllPairsTest(count ~ spray, InsectSprays, p.adjust="holm")
#' summary(out)
#' summaryGroup(out)
#'
#' @seealso
#' \code{\link[BWStest]{bws_test}}.
#' @export
bwsAllPairsTest <- function(x, ...) UseMethod("bwsAllPairsTest")

#' @rdname bwsAllPairsTest
#' @method bwsAllPairsTest default
#' @aliases bwsAllPairsTest.default
#' @template one-way-parms
#' @param p.adjust.method method for adjusting p values (see \code{\link{p.adjust}}).
#' @param method a character string specifying the test statistic to use. Defaults to \code{BWS}.
#' @importFrom stats p.adjust
#' @importFrom stats p.adjust.methods
#' @importFrom stats pairwise.table
#' @importFrom stats complete.cases
#' @importFrom BWStest bws_stat
#' @importFrom BWStest bws_cdf
#' @importFrom BWStest murakami_stat
#' @importFrom BWStest murakami_cdf
#' @export
bwsAllPairsTest.default <-
    function(x, g, method = c("BWS", "Murakami"), p.adjust.method = p.adjust.methods, ...)
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
                                        #
        if (is.null(x$p.adjust.method)){
            p.adjust.method <- p.adjust.methods[1]
        } else {
            p.adjust.method <- x$p.adjust.method
        }
        method <- x$method
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

    N <- length(x)
    if (N < 2)
        stop("not enough observations")

    p.adjust.method <- match.arg(p.adjust.method)
    method <- match.arg(method)
    method <- ifelse(method == "Murakami", "B1", method)

    METHOD <- switch(method,
                     "BWS" = "BWS All-Pairs Test",
                     "B1" = "Murakami's modified All-Pairs BWS Test")

##    compare.levels <- function(i, j) {
##        xi <- x[as.integer(g) == i]
##        xj <- x[as.integer(g) == j]
##        do.call("bws_test", list(x=xi,
##                                 y=xj,
##                                 method=method,
##                                 alternative="two.sided"))$p.value
##    }
##    bws.stats <- function(i, j) {
##        xi <- x[as.integer(g) == i]
##        xj <- x[as.integer(g) == j]
##        bws_stat(x=xi, y=xj)
##
##       # do.call("bws_test",
##       #         list(x=xi,
##       #              y=xj,
##       #              method=method,
##       #              alternative="two.sided"))$statistic
##    }
##    PVAL <- pairwise.table(compare.levels, levels(g), p.adjust.method)
##    STAT <- pairwise.table(compare.stats, levels(g), "none")

    n <- tapply(x, g, length)
    pval <- rep(NA, k * (k-1) / 2)
    stat <- pval
    ij <- 0
    if (method == "BWS") {
        for (i in 2:k)
        {
            for (j in (1:(i-1)))
            {
                ij <- ij + 1
                stat[ij] <- bws_stat(x = x[as.integer(g) == i], y= x[as.integer(g) == j])
            }
        }
        pval <- bws_cdf(b = stat, maxj = 3, lower_tail = FALSE)
        pval <- p.adjust(pval, method = p.adjust.method)

    } else {
        for (i in 2:k)
        {
            for (j in (1:(i-1)))
            {
                ij <- ij + 1
                stat[ij] <- murakami_stat(x = x[as.integer(g) == i],
                                          y = x[as.integer(g) == j],
                                          flavor = 1)
                pval[ij] <- murakami_cdf(B = stat[ij],
                                         n1 = n[i],
                                         n2 = n[j],
                                         flavor = 1,
                                         lower_tail = FALSE)
            }
        }
        pval <- p.adjust(pval, method = p.adjust.method)
    }

    STAT <- matrix(NA, ncol=(k-1), nrow=(k-1))
    rownames(STAT) <- levels(g)[2:k]
    colnames(STAT) <- levels(g)[1:(k-1)]
    PVAL <- STAT

    ij <- 0
    for(i in 2:k){
        for (j in (1:(i-1))){
            ij <- ij + 1
            STAT[(i-1),j] <- stat[ij]
            PVAL[(i-1),j] <- pval[ij]
            }
    }

    ans <- list(method = METHOD,
                data.name = DNAME,
                p.value = PVAL,
                p.adjust.method=p.adjust.method,
                dist = ifelse(method == "BWS", "B", "B*"),
                statistic = STAT,
                alternative = "two.sided",
                model = data.frame(x=x, g=g))
    class(ans) <- "PMCMR"
    ans
}

#' @rdname bwsAllPairsTest
#' @method bwsAllPairsTest formula
#' @aliases bwsAllPairsTest
#' @template one-way-formula
#' @export
bwsAllPairsTest.formula <-
    function(formula, data, subset, na.action,
             method = c("BWS", "Murakami"),
             p.adjust.method = p.adjust.methods, ...)
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
    p.adjust.method <- match.arg(p.adjust.method)
    method = match.arg(method)
    names(mf) <- NULL
    y <- do.call("bwsAllPairsTest",
                 c(as.list(mf), p.adjust.method = p.adjust.method, method=method))
    y$data.name <- DNAME
    y
}
