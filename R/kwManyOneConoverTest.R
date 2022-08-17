## kwManyOneConoverTest.R
## Part of the R package: PMCMR
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

####
#### Requires package mvtnorm

#' @name kwManyOneConoverTest
#' @title Conover's Many-to-One Rank Comparison Test
#' @description
#' Performs Conover's non-parametric many-to-one comparison
#' test for Kruskal-type ranked data.
#'
#' @details
#' For many-to-one comparisons (pairwise comparisons with one control)
#' in an one-factorial layout with non-normally distributed
#' residuals Conover's non-parametric test can be performed.
#' Let there be \eqn{k} groups including the control,
#' then the number of treatment levels is \eqn{m = k - 1}.
#' Then \eqn{m} pairwise comparisons can be performed between
#' the \eqn{i}-th treatment level and the control.
#' H\eqn{_i: \theta_0 = \theta_i} is tested in the two-tailed case against
#' A\eqn{_i: \theta_0 \ne \theta_i, ~~ (1 \le i \le m)}.
#'
#' If \code{p.adjust.method == "single-step"} is selected,
#' the \eqn{p}-values will be computed
#' from the multivariate \eqn{t} distribution. Otherwise,
#' the \eqn{p}-values are computed from the \eqn{t}-distribution using
#' any of the \eqn{p}-adjustment methods as included in \code{\link{p.adjust}}.
#'
#' @inherit cuzickTest note
#'
#' @inherit kwAllPairsConoverTest references
#'
#' @template class-PMCMR
#'
#' @concept kruskalranks
#' @keywords nonparametric
#'
#' @seealso
#' \code{\link{pmvt}}, \code{\link{TDist}}, \code{\link{kruskalTest}},
#' \code{\link{kwManyOneDunnTest}}, \code{\link{kwManyOneNdwTest}}
#' @example examples/kwManyOneMC.R
#' @export
kwManyOneConoverTest <- function(x, ...) UseMethod("kwManyOneConoverTest")

#' @rdname kwManyOneConoverTest
#' @method kwManyOneConoverTest default
#' @aliases kwManyOneConoverTest.default
#' @template one-way-parms
#' @param alternative the alternative hypothesis. Defaults to \code{two.sided}.
#' @param p.adjust.method  method for adjusting p values
#'    (see \code{\link{p.adjust}}).
#' @importFrom stats pt
#' @importFrom stats p.adjust.methods
#' @importFrom stats p.adjust
#' @importFrom stats complete.cases
#' @importFrom mvtnorm pmvt
#' @export
kwManyOneConoverTest.default <-
function(x, g, alternative = c("two.sided", "greater", "less"), p.adjust.method = c("single-step" ,p.adjust.methods),...){
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
        if (is.null(x$alternative)){
            alternative <- "two.sided"
        } else {
            alternative <- x$alternative
        }
        if(is.null(x$p.adjust.method)){
            p.adjust.method <- "single-step"
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
    alternative <- match.arg(alternative)
    p.adjust.method <- match.arg(p.adjust.method)

    ## rank version
    x.rank <- rank(x)
    R.bar <- tapply(x.rank, g, mean, na.rm = TRUE)
    R.n <- tapply(!is.na(x), g, length)
    k <- nlevels(g)
    n <- sum(R.n)

    METHOD <- "Conover's many-to-one test"

    ## Kruskal-Wallis statistic
    H <- HStat(x.rank, g)
    C <- gettiesKruskal(x.rank)
    H.cor <- H / C

    if (C == 1) {
        S2 <- n * (n + 1) / 12
    } else {
        warning("Ties are present. Quantiles were corrected for ties.")
        S2 <-   ( 1 / (n - 1)) * (sum(x.rank^2) - (n * (((n + 1)^2) / 4)))
    }

    compare.stats <- function(i) {
        dif <- R.bar[i] - R.bar[1]
        B <- (1 / R.n[i] + 1 / R.n[1])
        D <- (n - 1 - H.cor) / (n - k)
        tval <- dif / sqrt(S2 * B * D)
        return(tval)
    }

    if (p.adjust.method != "single-step"){
        df <- n - k
        STATISTIC <- rep(NA, k - 1)
        for (j in 2:k) {
            STATISTIC[j-1] <- compare.stats(j)
        }

        PVAL <- switch(alternative,
                       "two.sided" =
                           2 * pt(abs(STATISTIC), df=df, lower.tail=FALSE),
                       "greater" =
                           pt(STATISTIC, df=df, lower.tail=FALSE),
                       "less" =
                           pt(STATISTIC, df=df)
                       )

        PVAL <- p.adjust(PVAL, method = p.adjust.method)
        PARMS <- c(df = df)
        DIST <- "t"

    } else {
        ## correlation matrix

        n0 <- R.n[1]
        nn <- R.n[2:k]
        kk <- k - 1
        corr <- matrix(0, nrow = kk, ncol = kk)
        corr <- diag(kk)
        for ( i in 1:(kk-1)){
            for (j in (i+1):kk){
                corr[i,j] <- ((nn[i] * nn[j]) /
                              ((nn[i] + n0) * (nn[j] + n0)))^(1/2)
                corr[j,i] <- corr[i, j]
            }
        }

        df <- length(x) - k
        STATISTIC <- rep(NA, k - 1)

        ## Get statistic values
        for (j in 2:k) {
            STATISTIC[j-1] <- compare.stats(j)
        }

        ## Get p-values from multivariate t-distribution
        if (alternative == "two.sided") {
            PVAL <- sapply(STATISTIC,
                           function(x)  1 - pmvt(lower= -rep(abs(x),kk),
                                                 upper=rep(abs(x), kk), df=df,
                                                 corr=corr))
        } else if (alternative == "greater"){
            PVAL <- sapply(STATISTIC,
                           function(x)  1 - pmvt(lower= -Inf,
                                                 upper=rep(x, kk), df=df,
                                                 corr=corr))
        } else {
            PVAL <- sapply(STATISTIC,
                           function(x)  1 - pmvt(lower= rep(x, kk),
                                                 upper=Inf, df=df, corr=corr))
        }

        ## Names

        PARMS <- c(k = kk, df = df)
        DIST <- "t"
    }
    LNAME <- levels(g)[2:k]
    PSTAT <- matrix(data=STATISTIC, nrow = (k-1), ncol = 1,
                    dimnames = list(LNAME, levels(g)[1]))
    PVAL <- matrix(data=PVAL, nrow = (k-1), ncol = 1,
                   dimnames = list(LNAME, levels(g)[1]))
    MODEL <- data.frame(x, g)
    ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL,
                statistic = PSTAT, p.adjust.method = p.adjust.method,
                dist = DIST, model = MODEL, alternative = alternative,
                parameter = PARMS)
    class(ans) <- "PMCMR"
    ans
}

#' @rdname kwManyOneConoverTest
#' @method kwManyOneConoverTest formula
#' @aliases kwManyOneConoverTest.formula
#' @template one-way-formula
#' @export
kwManyOneConoverTest.formula <-
function(formula, data, subset, na.action, alternative = c("two.sided", "greater", "less"), p.adjust.method = c("single-step" ,p.adjust.methods),...)
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
    names(mf) <- NULL
    y <- do.call("kwManyOneConoverTest", c(as.list(mf),
                                        alternative = alternative,
                                        p.adjust.method=p.adjust.method))
    y$data.name <- DNAME
    y
}
