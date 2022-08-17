##  normalScoresManyOneTest.R
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

#' @title Lu-Smith Many-One Comparisons Normal Scores Test
#' @description Performs Lu-Smith multiple comparison
#' normal scores test with one control.
#' @details
#' For many-to-one comparisons in an one-factorial layout
#' with non-normally distributed residuals Lu and Smith's
#' normal scores transformation can be used prior to
#' a many-to-one comparison test. A total of \eqn{m = k-1}
#' hypotheses can be tested. The null hypothesis
#' H\eqn{_{i}: F_0(x) = F_i(x)} is tested in the two-tailed test
#' against the alternative
#' A\eqn{_{i}: F_0(x) \ne F_i(x), ~~ 1 \le i \le k-1}.
#' For \code{p.adjust.method = "single-step"} the
#' multivariate t distribution is used to calculate
#' p-values (see \code{\link[mvtnorm]{pmvt}}). Otherwise, the
#' t-distribution is used for the calculation of p-values
#' with a latter p-value adjustment as
#' performed by \code{\link{p.adjust}}.
#'
#' @name normalScoresManyOneTest
#' @template class-PMCMR
#' @keywords htest nonparametric
#' @concept normalscores
#'
#' @inherit cuzickTest note
#' @inherit normalScoresAllPairsTest references
#'
#' @seealso
#' \code{\link{normalScoresTest}}, \code{\link{normalScoresAllPairsTest}}, \code{\link[SuppDists]{normOrder}}, \code{\link[mvtnorm]{pmvt}}.
#'
#' @examples
#' ## Data set PlantGrowth
#' ## Global test
#' normalScoresTest(weight ~ group, data = PlantGrowth)
#'
#' ## Lu-Smith's many-one comparison test
#' ans <- normalScoresManyOneTest(weight ~ group, data = PlantGrowth, p.adjust.method = "holm")
#' summary(ans)
#' @export
normalScoresManyOneTest <- function(x, ...) UseMethod("normalScoresManyOneTest")

#' @rdname normalScoresManyOneTest
#' @method normalScoresManyOneTest default
#' @aliases normalScoresManyOneTest.default
#' @template one-way-parms
#' @param alternative the alternative hypothesis. Defaults to \code{two.sided}.
#' @param p.adjust.method method for adjusting p values (see \code{\link{p.adjust}}).
#' @importFrom stats pt
#' @importFrom mvtnorm pmvt
#' @export
normalScoresManyOneTest.default <-
    function(x, g, alternative = c("two.sided", "greater", "less"),
             p.adjust.method = c("single-step", p.adjust.methods), ...)
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
        if (is.null(x$alternative)){
            alternative <- "two.sided"
        } else {
            alternative <- x$alternative
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
    alternative <- match.arg(alternative)

    n <- length(x)
    if (n < 2)
        stop("not enough observations")
    r <- rank(x, ties.method = "random")

    # transform to z-scores
    eij <- normOrder(n)
    zscores <- eij[r]
    ei <- tapply(zscores, g, sum)
    ni <- tapply(zscores, g, length)
    s2 <- (1 / (n - 1)) * sum(zscores^2)
    STATISTIC <- (1 / s2) * sum(ei^2 / ni)
    PARAMETER <- k - 1
    A.mn <- ei / ni

    compare.stats <- function(i) {
        dif <- A.mn[i] - A.mn[1]
        B <- (1 / ni[i] + 1 / ni[1])
        tval <- dif / sqrt(s2 * (n-1-STATISTIC)/(n-k) * B)
        return(tval)
    }


    pstat <- sapply(2:k, function(i) compare.stats(i))

    if (p.adjust.method != "single-step"){
        if (alternative == "two.sided"){
            pval <- 2 * pt(abs(pstat), df=n-k, lower.tail = FALSE)
        } else if (alternative == "greater"){
            pval <- pt(pstat, df=n-k, lower.tail = FALSE)
        } else {
            pval <- pt(pstat, df=n-k)
        }
        padj <- p.adjust(pval, method = p.adjust.method)
    } else {
        ## use function pmvt of package mvtnorm
        m <- k - 1
        df <- n - k
        ## correlation matrix
        #ni <- tapply(x, g, length)
        n0 <- ni[1]
        nn <- ni[2:k]
        cr <- diag(m)
        ## cr also for unequal sample sizes
        for (i in 1:(m-1)){
            for (j in (i+1):m){
                cr[i,j] <- ((nn[i] * nn[j]) /
                            ((nn[i] + n0) * (nn[j] + n0)))^(1/2)
                cr[j,i] <- cr[i, j]
                cr[j,i] <- cr[i, j]
            }
        }
        if (alternative == "two.sided"){
            padj <- sapply(pstat, function(x)
                1 - pmvt(lower = -rep(abs(x), m),
                         upper = rep(abs(x), m),
                         df = df,
                         corr = cr))
        } else if (alternative == "greater"){
            padj <- sapply(pstat, function(x)
                1 - pmvt(lower = -Inf,
                         upper = rep(x, m),
                         df = df,
                         corr = cr))
        } else {
            padj <- sapply(pstat, function(x)
                1 - pmvt(lower = rep(x, m),
                         upper = Inf,
                         df = df,
                         corr = cr))
        }
    }

    GRPNAME <- levels(g)
    PVAL <- cbind(padj)
    PSTAT <- cbind(pstat)
    colnames(PVAL) <- GRPNAME[1]
    colnames(PSTAT) <- GRPNAME[1]
    rownames(PVAL) <- GRPNAME[2:k]
    rownames(PVAL) <- GRPNAME[2:k]
    METHOD <- paste0("Lu's normal scores test\n",
                     "\t\tfor multiple comparisons with one control")
    MOD <- data.frame(x = x, g= g)
    ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL,
                statistic = PSTAT, p.adjust.method = p.adjust.method,
                dist = "t", model = MOD, alternative = alternative)
    class(ans) <- "PMCMR"
    return(ans)
}

#' @rdname normalScoresManyOneTest
#' @method normalScoresManyOneTest formula
#' @aliases normalScoresManyOneTest.formula
#' @template one-way-formula
#' @export
normalScoresManyOneTest.formula <-
    function(formula, data, subset, na.action,
             alternative = c("two.sided", "greater", "less"),
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
    alternative <- match.arg(alternative)
    names(mf) <- NULL
    y <- do.call("normalScoresManyOneTest",
                 c(as.list(mf), alternative = alternative,
                   p.adjust.method = p.adjust.method))
    y$data.name <- DNAME
    y
}
