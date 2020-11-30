## kwManyOneDunnTest.R
## Part of the R package: PMCMR
##
## Copyright (C) 2015-2018 Thorsten Pohlert
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

#' @rdname kwManyOneDunnTest
#' @title Dunn's Many-to-One Rank Comparison Test
#' @description
#' Performs Dunn's non-parametric many-to-one comparison
#' test for Kruskal-type ranked data.
#'
#' @details
#' For many-to-one comparisons (pairwise comparisons with one control)
#' in an one-factorial layout with non-normally distributed
#' residuals Dunn's non-parametric test can be performed.
#' Let there be \eqn{k} groups including the control,
#' then the number of treatment levels is \eqn{m = k - 1}.
#' Then \eqn{m} pairwise comparisons can be performed between
#' the \eqn{i}-th treatment level and the control.
#' H\eqn{_i: \theta_0 = \theta_i} is tested in the two-tailed case against
#' A\eqn{_i: \theta_0 \ne \theta_i, ~~ (1 \le i \le m)}.
#'
#' If \code{p.adjust.method == "single-step"} is selected,
#' the \eqn{p}-values will be computed
#' from the multivariate normal distribution. Otherwise,
#' the \eqn{p}-values are computed from the standard normal distribution using
#' any of the \eqn{p}-adjustment methods as included in \code{\link{p.adjust}}.
#'
#' @inherit kwAllPairsDunnTest references
#'
#' @template class-PMCMR
#' @concept kruskalranks
#' @keywords nonparametric
#' @seealso
#' \code{\link{pmvnorm}}, \code{\link{TDist}}, \code{\link{kruskalTest}},
#' \code{\link{kwManyOneConoverTest}}, \code{\link{kwManyOneNdwTest}}
#' @example examples/kwManyOneMC.R
#' @export
kwManyOneDunnTest <- function(x, ...) UseMethod("kwManyOneDunnTest")

#' @rdname kwManyOneDunnTest
#' @aliases kwManyOneDunnTest.default
#' @method kwManyOneDunnTest default
#' @template one-way-parms
#' @param alternative the alternative hypothesis. Defaults to \code{two.sided}.
#' @param p.adjust.method  method for adjusting p values
#'    (see \code{\link{p.adjust}}).
#' @importFrom mvtnorm pmvnorm
#' @importFrom stats pnorm
#' @importFrom stats p.adjust.methods
#' @importFrom stats p.adjust
#' @importFrom stats complete.cases
#' @export
kwManyOneDunnTest.default <-
    function(x, g, alternative = c("two.sided", "greater", "less"),
             p.adjust.method = c("single-step", p.adjust.methods), ...)
{

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
        ##
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

    x.rank <- rank(x)
    R.bar <- tapply(x.rank, g, mean,na.rm=T)
    R.n <- tapply(!is.na(x), g, length)
    g.unique <- unique(g)
    k <- length(g.unique)
    n <- sum(R.n)

    METHOD <- "Dunn's many-to-one test"

    ## get the ties
    C <- gettiesDunn(x.rank)

    if (C != 0) warning("Ties are present. z-quantiles were corrected for ties.")
    ## mean Rsum of controll is in R.bar[1]
    compare.stats <- function(j) {
            ##dif <- abs(R.bar[1] - R.bar[j])
        dif <- R.bar[j] - R.bar[1]
        A <- n * (n+1) / 12
        B <- (1 / R.n[1] + 1 / R.n[j])
        zval <- dif / sqrt((A - C) * B)
        return(zval)
    }

    PSTATv <- rep(NA, k-1)
    for (j in 2:k) {PSTATv[j-1] <- compare.stats(j)}

    if (p.adjust.method != "single-step"){
        ## unadjusted p-values
        if (alternative == "two.sided"){
            PVALv <- 2 * pnorm(abs(PSTATv), lower.tail = FALSE)
        } else if (alternative == "greater"){
            PVALv <- pnorm(PSTATv, lower.tail = FALSE)
        } else {
            PVALv <- pnorm(PSTATv)
        }
        ## adjusted p-values
        PADJv <- p.adjust(PVALv, method = p.adjust.method)
    } else {
        ## use function pmvt of package mvtnorm
        pstat <- as.vector(PSTATv)
        m <- k - 1
        df <- n - k
        ## correlation matrix
        ni <- tapply(x, g, length)
        n0 <- ni[1]
        nn <- ni[2:k]
        cr <- diag(m)
        for (i in 1:(m-1)){
            for (j in (i+1):m){
				 #cr[i,j] <- 0.5
				 #cr[j,i] <- 0.5
                cr[i,j] <- ((nn[i] * nn[j]) /
                            ((nn[i] + n0) * (nn[j] + n0)))^(1/2)
                cr[j,i] <- cr[i, j]
                cr[j,i] <- cr[i, j]

            }
        }
        if (alternative == "two.sided"){
            pvalv <- sapply(pstat, function(x)
                1 - pmvnorm(lower = -rep(abs(x), m),
                            upper = rep(abs(x), m),
                            corr = cr))
        } else if (alternative == "greater"){
            pvalv <- sapply(pstat, function(x)
                1 - pmvnorm(lower = -Inf,
                            upper = rep(x, m),
                            corr = cr))
        } else {
            pvalv <- sapply(pstat, function(x)
                1 - pmvnorm(lower = rep(x, m),
                            upper = Inf,
                            corr = cr))
        }
        PADJv <- pvalv
    }
    LNAME <- levels(g)[2:k]
    MODEL <- data.frame(x = x, g = g)
    ## build matrix
    PSTAT <- matrix(data=PSTATv, nrow = (k-1), ncol = 1,
                    dimnames = list(LNAME, levels(g)[1]))
    PVAL <- matrix(data=PADJv, nrow = (k-1), ncol = 1,
                   dimnames = list(LNAME, levels(g)[1]))
    ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL,
                statistic = PSTAT, p.adjust.method = p.adjust.method,
                alternative = alternative, dist = "z", model = MODEL)
    class(ans) <- "PMCMR"
    ans
}

#' @rdname kwManyOneDunnTest
#' @aliases kwManyOneDunnTest.formula
#' @method kwManyOneDunnTest formula
#' @template one-way-formula
#' @export
kwManyOneDunnTest.formula <-
    function(formula, data, subset, na.action, alternative = c("two.sided", "greater", "less"),
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
    names(mf) <- NULL
    alternative <- match.arg(alternative)
    p.adjust.method <- match.arg(p.adjust.method)
    y <- do.call("kwManyOneDunnTest", c(as.list(mf), alternative = alternative, p.adjust.method = p.adjust.method))
    y$data.name <- DNAME
    y
}

