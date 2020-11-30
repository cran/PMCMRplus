## uryWigginsHochbergTest.R
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

#' @name uryWigginsHochbergTest
#' @title Ury, Wiggins, Hochberg Test
#' @description
#' Performs Ury-Wiggins and Hochberg's all-pairs comparison test
#' for normally distributed data with unequal variances.
#'
#' @template class-PMCMR
#'
#' @details
#' For all-pairs comparisons in an one-factorial layout
#' with normally distributed residuals but unequal groups variances
#' the tests of Ury-Wiggins and Hochberg can be performed.
#' Let \eqn{X_{ij}} denote a continuous random variable
#' with the \eqn{j}-the realization (\eqn{1 \le j \le n_i})
#' in the \eqn{i}-th group (\eqn{1 \le i \le k}). Furthermore, the total
#' sample size is \eqn{N = \sum_{i=1}^k n_i}. A total of \eqn{m = k(k-1)/2}
#' hypotheses can be tested: The null hypothesis is
#' H\eqn{_{ij}: \mu_i = \mu_j ~~ (i \ne j)} is tested against the alternative
#' A\eqn{_{ij}: \mu_i \ne \mu_j} (two-tailed). Ury-Wiggins and Hochberg
#' all-pairs test statistics are given by
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
#'  \alpha,
#' }{%
#'  SEE PDF
#' }
#'
#' with Welch's approximate equation for degree of freedom as
#'
#' \deqn{
#'  v_{ij} = \frac{\left( s^2_i / n_i + s^2_j / n_j \right)^2}
#'  {s^4_i / n^2_i \left(n_i - 1\right) + s^4_j / n^2_j \left(n_j - 1\right)}.
#' }{%
#'  SEE PDF
#' }
#'
#' The p-values are computed from the \code{\link[stats]{TDist}}-distribution.
#' The type of test depends
#' on the selected p-value adjustment method (see also \code{\link{p.adjust}}):
#'
#' \describe{
#' \item{bonferroni}{the Ury-Wiggins test is performed with Bonferroni adjusted
#' p-values.}
#' \item{hochberg}{the Hochberg test is performed with Hochberg's adjusted
#' p-values}.
#' }
#'
#' @references
#' Hochberg, Y. (1976) A Modification of the T-Method of Multiple
#' Comparisons for a One-Way Layout With Unequal Variances,
#' \emph{Journal of the American Statistical Association} \bold{71}, 200--203.
#'
#' Ury, H. and Wiggins, A. D. (1971) Large Sample and Other
#' Multiple Comparisons Among Means, \emph{British Journal of
#' Mathematical and Statistical Psychology} \bold{24}, 174--194.
#'
#'
#' @examples
#' fit <- aov(weight ~ feed, chickwts)
#' shapiro.test(residuals(fit))
#' bartlett.test(weight ~ feed, chickwts) # var1 = varN
#' anova(fit)
#'
#' ## also works with fitted objects of class aov
#' res <- uryWigginsHochbergTest(fit)
#' summary(res)
#' summaryGroup(res)
#'
#' @seealso
#' \code{\link{dunnettT3Test}} \code{\link{tamhaneT2Test}} \code{\link[stats]{TDist}}
#' @importFrom stats pt
#' @importFrom stats complete.cases
#' @importFrom stats var
#' @importFrom stats pairwise.table
#' @export
uryWigginsHochbergTest <- function(x, ...) UseMethod("uryWigginsHochbergTest")

#' @rdname uryWigginsHochbergTest
#' @method uryWigginsHochbergTest default
#' @aliases uryWigginsHochbergTest.default
#' @template one-way-parms-aov
#' @param p.adjust.method method for adjusting p values
#'    (see \code{\link{p.adjust}}).
#' @export
uryWigginsHochbergTest.default <-
function(x, g, p.adjust.method = p.adjust.methods, ...){
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
        p.adjust.method <- x$p.adjust.method
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

    p.adjust.method = match.arg(p.adjust.method)

    ## prepare uryWigginsHochberg test
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

    PSTAT <- pairwise.table(compare.stats,levels(g), p.adjust.method="none" )

    compare.levels <- function(i,j) {
        dif <- xi[i] - xi[j]
        A <- (s2i[i] / ni[i] + s2i[j] / ni[j])
        tval <- dif / sqrt(A)
        df <- A^2 / (s2i[i]^2 / (ni[i]^2 * (ni[i] - 1)) +
                    s2i[j]^2 / (ni[j]^2 * (ni[j] - 1)))
        ## two sided
        pval <- 2 * pt(abs(tval), df = df, lower.tail = FALSE)
        pval <- min(1, pval)
        return(pval)
    }

    PVAL <-  pairwise.table(compare.levels, levels(g),
                            p.adjust.method= p.adjust.method)
    if(p.adjust.method == "bonferroni"){
        METHOD <- "Ury and Wiggins test for unequal variances"
    } else if (p.adjust.method == "hochberg"){
        METHOD <- "Hochberg's H test for unequal variances"
    } else {
        METHOD <- "test for unequal variances"
    }
    MODEL <- data.frame(x, g)
    DIST <- "t"
    ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL,
                statistic = PSTAT, p.adjust.method = p.adjust.method,
                model = MODEL, dist = DIST, alternative = "two.sided")
    class(ans) <- "PMCMR"
    ans
}

#' @rdname uryWigginsHochbergTest
#' @method uryWigginsHochbergTest formula
#' @aliases uryWigginsHochbergTest.formula
#' @template one-way-formula
#' @export
uryWigginsHochbergTest.formula <-
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
    names(mf) <- NULL

    p.adjust.method = match.arg(p.adjust.method)

    y <- do.call("uryWigginsHochbergTest", c(as.list(mf),
                                   p.adjust.method = p.adjust.method))
    y$data.name <- DNAME
    y
}

#' @rdname uryWigginsHochbergTest
#' @aliases uryWigginsHochbergTest.aov
#' @method uryWigginsHochbergTest aov
# @param obj A fitted model object, usually an \link[stats]{aov} fit.
#' @export
uryWigginsHochbergTest.aov <- function(x, p.adjust.method = p.adjust.methods, ...) {
    model <- x$model
    DNAME <- paste(names(model), collapse = " by ")
    names(model) <- c("x", "g")
    p.adjust.method <- match.arg(p.adjust.method)
    parms <- c(as.list(model), list(p.adjust.method = p.adjust.method))
    y <- do.call("uryWigginsHochbergTest", parms)
    y$data.name <- DNAME
    y
}

