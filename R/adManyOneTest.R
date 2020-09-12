#  adManyOneTest.R
#
#  Copyright (C) 2017, 2018 Thorsten Pohlert
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

#' @title Anderson-Darling Many-To-One Comparison Test
#' @description Performs Anderson-Darling many-to-one comparison test.
#' @details
#' For many-to-one comparisons (pairwise comparisons with one control)
#' in an one-factorial layout with non-normally distributed
#' residuals Anderson-Darling's non-parametric test can be performed.
#' Let there be \eqn{k} groups including the control,
#' then the number of treatment levels is \eqn{m = k - 1}.
#' Then \eqn{m} pairwise comparisons can be performed between
#' the \eqn{i}-th treatment level and the control.
#' H\eqn{_i: F_0 = F_i} is tested in the two-tailed case against
#' A\eqn{_i: F_0 \ne F_i, ~~ (1 \le i \le m)}.
#'
#' This function is a wrapper function that sequentially
#' calls \code{adKSampleTest} for each pair.
#' The calculated p-values for \code{Pr(>|T2N|)}
#' can be adjusted to account for Type I error inflation
#' using any method as implemented in \code{\link{p.adjust}}.
#' @name adManyOneTest
#' @template class-PMCMR
#' @keywords htest nonparametric
#'
#' @inherit adAllPairsTest references
#'
#' @seealso
#' \code{\link{adKSampleTest}}, \code{\link{adAllPairsTest}},
#' \code{\link[kSamples]{ad.pval}}.
#' @examples
#' ## Data set PlantGrowth
#' ## Global test
#' adKSampleTest(weight ~ group, data = PlantGrowth)
#'
#' ##
#' ans <- adManyOneTest(weight ~ group,
#'                              data = PlantGrowth,
#'                              p.adjust.method = "holm")
#' summary(ans)
#' @export
adManyOneTest <- function(x, ...) UseMethod("adManyOneTest")

#' @rdname adManyOneTest
#' @method adManyOneTest default
#' @aliases adManyOneTest.default
#' @template one-way-parms
#' @param p.adjust.method method for adjusting
#' p values (see \code{\link{p.adjust}}).
#' @export
adManyOneTest.default <-
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

    p.adjust.method <- match.arg(p.adjust.method)

    n <- length(x)
    if (n < 2)
        stop("not enough observations")
    r <- rank(x)


    METHOD <- "Anderson-Darling Two-Sample Test"
    compare.levels <- function(j) {
        xi <- x[as.integer(g) == 1]
        xj <- x[as.integer(g) == j]
        adKSampleTest(list(xi, xj), dist="estimated", ...)$p.value
    }
    compare.stats <- function(j) {
        xi <- x[as.integer(g) == 1]
        xj <- x[as.integer(g) == j]
        adKSampleTest(list(xi, xj), dist="estimated", ...)$statistic
    }
    pval <- sapply(2:k, function(j) compare.levels(j))
    pval <- p.adjust(pval, method = p.adjust.method)
    stat <- sapply(2:k, function(j) compare.stats(j))

    ## Prepare output
    PVAL <- cbind(pval)
    colnames(PVAL) <- levels(g)[1]
    rownames(PVAL) <- levels(g)[2:k]
    STAT <- cbind(stat)
    colnames(STAT) <- colnames(PVAL)
    rownames(STAT) <- rownames(PVAL)

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

#' @rdname adManyOneTest
#' @method adManyOneTest formula
#' @aliases adManyOneTest.formula
#' @template one-way-formula
#' @export
adManyOneTest.formula <-
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
    y <- do.call("adManyOneTest",
                 c(as.list(mf), p.adjust.method = p.adjust.method))
    y$data.name <- DNAME
    y
}
