##  bwsManyOneTest.R
##
##  Copyright (C) 2017 Thorsten Pohlert
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


#' @title BWS Many-To-One Comparison Test
#' @description Performs Baumgartner-Weiß-Schindler many-to-one comparison test.
#' @details
#' For many-to-one comparisons (pairwise comparisons with one control)
#' in an one-factorial layout with non-normally distributed
#' residuals Baumgartner-Weiß-Schindler's non-parametric test can be performed.
#' Let there be \eqn{k} groups including the control,
#' then the number of treatment levels is \eqn{m = k - 1}.
#' Then \eqn{m} pairwise comparisons can be performed between
#' the \eqn{i}-th treatment level and the control.
#' H\eqn{_i: F_0 = F_i} is tested in the two-tailed case against
#' A\eqn{_i: F_0 \ne F_i, ~~ (1 \le i \le m)}.
#'
#' This function is a wrapper function that sequentially
#' calls \code{\link[BWStest]{bws_stat}} and \code{\link[BWStest]{bws_cdf}}
#' for each pair. For the default test method (\code{"BWS"}) the original
#' Baumgartner-Weiß-Schindler test statistic B and its corresponding Pr(>|B|)
#' is calculated. For \code{method == "BWS"} only a two-sided test is possible.
#'
#' For \code{method == "Murakami"} the modified BWS statistic
#' denoted B* and its corresponding Pr(>|B*|) is computed by sequentially calling
#' \code{\link[BWStest]{murakami_stat}} and \code{\link[BWStest]{murakami_cdf}}.
#' For \code{method == "Murakami"} only a two-sided test is possible.
#'
#' If \code{alternative == "greater"} then the alternative, if one
#' population is stochastically larger than the other is tested:
#' H\eqn{_i: F_0 = F_i} against A\eqn{_i: F_0 \ge F_i, ~~ (1 \le i \le m)}.
#' The modified test-statistic B* according to Neuhäuser (2001) and its
#' corresponding Pr(>B*) or Pr(<B*) is computed by sequentally calling
#' \code{\link[BWStest]{murakami_stat}} and \code{\link[BWStest]{murakami_cdf}}
#' with \code{flavor = 2}.
#'
#' The p-values can be adjusted to account for Type I error
#' inflation using any method as implemented in \code{\link{p.adjust}}.
#'
#' @name bwsManyOneTest
#' @template class-PMCMR
#' @keywords htest nonparametric
#'
#' @references
#' Baumgartner, W., Weiss, P., Schindler, H. (1998) A nonparametric test for the
#' general two-sample problem, \emph{Biometrics} \bold{54}, 1129--1135.
#'
#' Murakami, H. (2006) K-sample rank test based on modified Baumgartner statistic and its power
#' comparison, \emph{J. Jpn. Comp. Statist.} \bold{19}, 1--13.
#'
#' Neuhäuser, M. (2001) One-side two-sample and trend tests based on a modified
#' Baumgartner-Weiss-Schindler statistic.
#' \emph{Journal of Nonparametric Statistics} \bold{13}, 729--739.
#'
#' @examples
#' out <- bwsManyOneTest(weight ~ group, PlantGrowth, p.adjust="holm")
#' summary(out)
#'
#' ## A two-sample test
#' set.seed(1245)
#' x <- c(rnorm(20), rnorm(20,0.3))
#' g <- gl(2, 20)
#' summary(bwsManyOneTest(x ~ g, alternative = "less", p.adjust="none"))
#' summary(bwsManyOneTest(x ~ g, alternative = "greater", p.adjust="none"))
#'
#' \dontrun{
#' ## Check with the implementation in package BWStest
#' BWStest::bws_test(x=x[g==1], y=x[g==2], alternative = "less")
#' BWStest::bws_test(x=x[g==1], y=x[g==2], alternative = "greater")
#' }
#' @seealso
#' \code{\link[BWStest]{murakami_stat}}, \code{\link[BWStest]{murakami_cdf}},
#'  \code{\link[BWStest]{bws_stat}}, \code{\link[BWStest]{bws_cdf}}.
#' @export
bwsManyOneTest <- function(x, ...) UseMethod("bwsManyOneTest")

#' @rdname bwsManyOneTest
#' @method bwsManyOneTest default
#' @aliases bwsManyOneTest.default
#' @template one-way-parms
#' @param alternative the alternative hypothesis. Defaults to \code{two.sided}.
#' @param p.adjust.method method for adjusting p values (see \code{\link{p.adjust}}).
#' @param method a character string specifying the test statistic to use. Defaults to \code{BWS}.
#' @importFrom stats p.adjust
#' @importFrom stats p.adjust.methods
#' @importFrom stats pairwise.table
#' @importFrom stats complete.cases
#' @importFrom BWStest murakami_stat
#' @importFrom BWStest murakami_cdf
#' @importFrom BWStest bws_stat
#' @importFrom BWStest bws_cdf
#' @export
bwsManyOneTest.default <-
    function(x, g, alternative = c("two.sided", "greater", "less"),
             method = c("BWS", "Murakami", "Neuhauser"),
             p.adjust.method = p.adjust.methods, ...)
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
        alternative <- x$alternative
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
    alternative <- match.arg(alternative)
    method <- match.arg(method)

    n <- tapply(x, g, length)

    ## Change method string
    if (alternative != "two.sided"){
        method <- "B2"
    } else if (method == "Murakami"){
        method <- "B1"
    } else if (method == "Neuhauser"){
        method <- "B1"
    }

    METHOD <- switch(method,
                     "BWS" = c("BWS Two-Sample Test",
                               " for multiple comparisons with one control"),
                     "B1" = c("Murakami's modified Two-Sample BWS Test",
                              " for multiple comparisons with one control"),
                     "B2" = c("Neuhauser's modified Two-Sample BWS Test",
                                     " for multiple comparisons with one control"))

    stat <- switch(method,
                   "BWS" = sapply(2:k, function(j)
                       do.call("bws_stat",
                               list(x = x[as.integer(g) == 1],
                                    y = x[as.integer(g) == j])
                               )
                       ),

                   "B1" = sapply(2:k, function(j)
                       do.call("murakami_stat",
                               list(x = x[as.integer(g) == 1],
                                    y = x[as.integer(g) == j],
                                    flavor = 1)
                               )
                       ),

                   "B2" = sapply(2:k, function(j)
                                            do.call("murakami_stat",
                                                    list(x = x[as.integer(g) == 1],
                                                         y = x[as.integer(g) == j],
                                                         flavor = 2)
                                                    )
                                 )
                   )

   ## pval <- switch(method,
    if (method == "BWS") {
        pval <- do.call("bws_cdf",
                        list(b = stat, maxj = 3, lower=FALSE)
                        )

    } else if (method == "B1"){
        pval <- sapply(2:k, function(j)
            do.call("murakami_cdf",
                    list(B = stat[j-1],
                         n1 = n[1],
                         n2 = n[j],
                         flavor = 1,
                         lower = FALSE)
                    )
            )

    } else if (alternative == "greater"){
        pval <- sapply(2:k, function(j)
            do.call("murakami_cdf",
                    list(B = stat[j-1],
                         n1 = n[j],
                         n2 = n[1],
                         flavor = 2,
                         lower=FALSE)
                    )
            )
    } else {
        pval <- sapply(2:k, function(j)
            do.call("murakami_cdf",
                    list(B = stat[j-1],
                         n1 = n[1],
                         n2 = n[j],
                         flavor = 2,
                         lower=TRUE)
                    )
            )
    }

    ## adjust p-values
    pval <- p.adjust(pval, method = p.adjust.method)

    ## Create matrices
    STAT <- cbind(stat)
    colnames(STAT) <- levels(g)[1]
    rownames(STAT) <- levels(g)[-1]
    PVAL <- cbind(pval)
    colnames(PVAL) <- colnames(STAT)
    rownames(PVAL) <- rownames(STAT)

    ans <- list(method = METHOD,
                data.name = DNAME,
                p.value = PVAL,
                p.adjust.method=p.adjust.method,
                dist = ifelse(method == "BWS", "B", "B*"),
                statistic = STAT,
                alternative = alternative,
                model = data.frame(x=x, g=g))
    class(ans) <- "PMCMR"
    ans
}

#' @rdname bwsManyOneTest
#' @method bwsManyOneTest formula
#' @aliases bwsManyOneTest
#' @template one-way-formula
#' @export
bwsManyOneTest.formula <-
    function(formula, data, subset, na.action,
             alternative = c("two.sided", "greater", "less"),
             method = c("BWS", "Murakami", "Neuhauser"),
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
    alternative <- match.arg(alternative)
    p.adjust.method <- match.arg(p.adjust.method)
    method = match.arg(method)
    names(mf) <- NULL
    y <- do.call("bwsManyOneTest",
                 c(as.list(mf), alternative = alternative,
                   p.adjust.method = p.adjust.method, method=method))
    y$data.name <- DNAME
    y
}
