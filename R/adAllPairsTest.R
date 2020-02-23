##  adAllPairsTest.R
##
##  Copyright (C) 2017, 2018 Thorsten Pohlert
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


#' @title Anderson-Darling All-Pairs Comparison Test
#' @description Performs Anderson-Darling all-pairs comparison test.
#' @details
#' For all-pairs comparisons in an one-factorial layout
#' with non-normally distributed residuals Anderson-Darling's
#' all-pairs comparison test can be used. A total of \eqn{m = k(k-1)/2}
#' hypotheses can be tested. The null hypothesis
#' H\eqn{_{ij}: F_i(x) = F_j(x)} is tested in the two-tailed test
#' against the alternative
#' A\eqn{_{ij}: F_i(x) \ne F_j(x), ~~ i \ne j}.
#'
#' This function is a wrapper function that sequentially
#' calls \code{adKSampleTest} for each pair.
#' The calculated p-values for \code{Pr(>|T2N|)}
#' can be adjusted to account for Type I error multiplicity
#' using any method as implemented in \code{\link{p.adjust}}.
#'
#' @name adAllPairsTest
#' @template class-PMCMR
#' @keywords htest nonparametric
#' @concept AllPairsComparison
#' @references
#' Scholz, F.W., Stephens, M.A. (1987) K-Sample Anderson-Darling Tests.
#' \emph{Journal of the American Statistical Association} \bold{82}, 918--924.
#'
#' @examples
#' adKSampleTest(count ~ spray, InsectSprays)
#'
#' out <- adAllPairsTest(count ~ spray, InsectSprays, p.adjust="holm")
#' summary(out)
#' summaryGroup(out)
#'
#' @seealso
#' \code{\link{adKSampleTest}}, \code{\link{adManyOneTest}},
#' \code{\link[kSamples]{ad.pval}}.
#' @export
adAllPairsTest <- function(x, ...) UseMethod("adAllPairsTest")

#' @rdname adAllPairsTest
#' @method adAllPairsTest default
#' @aliases adAllPairsTest.default
#' @template one-way-parms
#' @param p.adjust.method method for adjusting p values (see \code{\link{p.adjust}}).
#' @importFrom stats p.adjust
#' @importFrom stats p.adjust.methods
#' @importFrom stats pairwise.table
#' @importFrom stats complete.cases
#' @export
adAllPairsTest.default <-
    function(x, g, p.adjust.method = p.adjust.methods, ...)
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

    n <- length(x)
    if (n < 2)
        stop("not enough observations")

    p.adjust.method <- match.arg(p.adjust.method)

    ##   DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(g)))
    ##   g <- factor(g)
    METHOD <- "Anderson-Darling All-Pairs Test"
    compare.levels <- function(i, j) {
        xi <- x[as.integer(g) == i]
        xj <- x[as.integer(g) == j]
        adKSampleTest(list(xi, xj), dist="estimated", ...)$p.value
    }
    compare.stats <- function(i, j) {
        xi <- x[as.integer(g) == i]
        xj <- x[as.integer(g) == j]
        adKSampleTest(list(xi, xj), dist="estimated", ...)$statistic
    }
    PVAL <- pairwise.table(compare.levels, levels(g), p.adjust.method)
    STAT <- pairwise.table(compare.stats, levels(g), "none")
    ans <- list(method = METHOD,
                data.name = DNAME,
                p.value = PVAL,
                p.adjust.method=p.adjust.method,
                dist = "T2N",
                statistic = STAT,
                alternative = "two.sided",
                model = data.frame(x=x, g=g))
    class(ans) <- "PMCMR"
    ans
}


#' @rdname adAllPairsTest
#' @method adAllPairsTest formula
#' @aliases adAllPairsTest
#' @template one-way-formula
#' @export
adAllPairsTest.formula <-
    function(formula, data, subset, na.action,
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
    names(mf) <- NULL
    y <- do.call("adAllPairsTest",
                 c(as.list(mf), p.adjust.method = p.adjust.method))
    y$data.name <- DNAME
    y
}
