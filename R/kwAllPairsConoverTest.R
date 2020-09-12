# kwAllPairsConover.test.R
# Part of the R package: PMCMR
#
# Copyright (C) 2015-2018 Thorsten Pohlert
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

#' @name kwAllPairsConoverTest
#' @title Conover's All-Pairs Rank Comparison Test
#'
#' @description
#' Performs Conover's non-parametric all-pairs comparison test
#' for Kruskal-type ranked data.
#'
#' @details
#' For all-pairs comparisons in an one-factorial layout
#' with non-normally distributed residuals Conover's non-parametric test
#' can be performed. A total of \eqn{m = k(k-1)/2}
#' hypotheses can be tested. The null hypothesis
#' H\eqn{_{ij}: \mu_i(x) = \mu_j(x)} is tested in the two-tailed test
#' against the alternative
#' A\eqn{_{ij}: \mu_i(x) \ne \mu_j(x), ~~ i \ne j}.
#'
#' If \code{p.adjust.method == "single-step"} the p-values are computed
#' from the studentized range distribution. Otherwise,
#' the p-values are computed from the t-distribution using
#' any of the p-adjustment methods as included in \code{\link{p.adjust}}.
#'
#' @references
#' Conover, W. J, Iman,  R. L. (1979) \emph{On multiple-comparisons
#'  procedures}, Tech. Rep. LA-7677-MS, Los Alamos Scientific Laboratory.
#'
#' @template class-PMCMR
#' @keywords nonparametric
#'
#' @concept kruskalranks
#'
#' @seealso
#' \code{\link[stats]{Tukey}}, \code{\link[stats]{TDist}},
#' \code{\link[stats]{p.adjust}}, \code{\link{kruskalTest}},
#' \code{\link{kwAllPairsDunnTest}}, \code{\link{kwAllPairsNemenyiTest}}
#'
#' @example examples/kwAllPairsMC.R
#' @export
kwAllPairsConoverTest <- function(x, ...) UseMethod("kwAllPairsConoverTest")

#' @rdname kwAllPairsConoverTest
#' @method kwAllPairsConoverTest default
#' @aliases kwAllPairsConoverTest.default
#' @template one-way-parms
#' @param p.adjust.method method for adjusting p values
#'    (see \code{\link{p.adjust}}).
#' @importFrom stats pt
#' @importFrom stats ptukey
#' @importFrom stats p.adjust.methods
#' @importFrom stats p.adjust
#' @importFrom stats complete.cases
#' @export
kwAllPairsConoverTest.default <-
function(x, g, p.adjust.method = c("single-step", p.adjust.methods), ...){
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
    p.adjust.method <- match.arg(p.adjust.method)
    x.rank <- rank(x)
    R.bar <- tapply(x.rank, g, mean,na.rm=T)
    R.n <- tapply(!is.na(x), g, length)
    g.unique <- unique(g)
    k <- length(g.unique)
    n <- sum(R.n)
    getties <- function(x, n) {
        x.sorted <- sort(x)
        pos <- 1
        tiesum <- 0
        while (pos <= n) {
            val <- x.sorted[pos]
            nt <- length(!is.na(x.sorted[x.sorted==val]))
            pos <- pos + nt
            if (nt > 1){
                tiesum <- tiesum + nt^3  - nt
            }
        }
        C <- 1 - tiesum / (n^3 - n)
        C <- min(c(1,C))
        return(C)
    }
    METHOD <- "Conover's all-pairs test"
    C <- getties(x.rank, n)
    if (C != 1) warning("Ties are present. Quantiles were corrected for ties.")
    ## Kruskal-Wallis statistic
    H <- (12 / (n * (n + 1))) * sum(tapply(x.rank, g, "sum")^2 / R.n) -
        3 * (n + 1)
    H.cor <- H / C

    if (C == 1) {
        S2 <- n * (n + 1) / 12
    } else {
        S2 <-   ( 1 / (n - 1)) * (sum(x.rank^2) - (n * (((n + 1)^2) / 4)))
    }
    compare.stats <- function(i,j) {
        dif <- R.bar[i] - R.bar[j]
        B <- (1 / R.n[i] + 1 / R.n[j])
        D <- (n - 1 - H.cor) / (n - k)
        tval <- dif / sqrt(S2 * B * D)
        return(tval)
    }
    PSTAT <- pairwise.table(compare.stats,levels(g),
                            p.adjust.method="none" )

    compare.levels <- function(i,j) {
        dif <- abs(R.bar[i] - R.bar[j])
        B <- (1 / R.n[i] + 1 / R.n[j])
        D <- (n - 1 - H.cor) / (n - k)
        tval <- dif / sqrt(S2 * B * D)
        pval <- 2 * pt(abs(tval), df=n - k, lower.tail=FALSE)
        return(pval)
    }

    if (p.adjust.method == "single-step"){

        compare.tukey <- function(i,j) {
            dif <- abs(R.bar[i] - R.bar[j])
            B <- (1 / R.n[i] + 1 / R.n[j])
            D <- (n - 1 - H.cor) / (n - k)
            qval <- sqrt(2) * dif / sqrt(S2 * B * D)
            pval <- ptukey(abs(qval),
                           nmeans = k, df=(n - k),
                           lower.tail=FALSE)
            return(pval)
        }
        PVAL <- pairwise.table(compare.tukey,levels(g),
                               p.adjust.method="none" )
        PSTAT <- PSTAT * sqrt(2)
        PARMS <- c(k, n - k)
        names(PARMS) <- c("nmeans", "df")
        DIST <- "q"

    } else {
        PVAL <- pairwise.table(compare.levels,levels(g),
                               p.adjust.method=p.adjust.method )
        PARMS <- n - k
        names(PARMS) <- "df"
        DIST <- "t"
    }
    MODEL <- data.frame(x, g)
    ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL,
                statistic = PSTAT, p.adjust.method = p.adjust.method,
                dist = DIST, model = MODEL,
                parameter = PARMS)
    class(ans) <- "PMCMR"
    ans
}

#' @rdname kwAllPairsConoverTest
#' @method kwAllPairsConoverTest formula
#' @aliases kwAllPairsConoverTest.formula
#' @template one-way-formula
#' @export
kwAllPairsConoverTest.formula <-
function(formula, data, subset, na.action,
         p.adjust.method = c("single-step", p.adjust.methods), ...)
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
    y <- do.call("kwAllPairsConoverTest", c(as.list(mf),
                                                   p.adjust.method))
    y$data.name <- DNAME
    y
}
