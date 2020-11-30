## chaAllPairsNashimotoTest.R
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

#' @name chaAllPairsNashimotoTest
#' @title All-Pairs Comparisons for Simply Ordered Mean Ranksums
#'
#' @description
#' Performs Nashimoto and Wright's all-pairs comparison procedure
#' for simply ordered mean ranksums (NPT'-test and NPY'-test).
#'
#' According to the authors, the procedure shall only be
#' applied after Chacko's test (see \code{\link{chackoTest}}) indicates
#' global significance.
#'
#' @details
#' The modified procedure uses the property of a simple order,
#' \eqn{\theta_m' - \theta_m \le \theta_j - \theta_i \le \theta_l' - \theta_l
#' \qquad (l \le i \le m~\mathrm{and}~ m' \le j \le l')}.
#' The null hypothesis H\eqn{_{ij}: \theta_i = \theta_j} is tested against
#' the alternative A\eqn{_{ij}: \theta_i < \theta_j} for any
#' \eqn{1 \le i < j \le k}.
#'
#' Let \eqn{R_{ij}} be the rank of \eqn{X_{ij}},
#' where \eqn{X_{ij}} is jointly ranked
#' from \eqn{\left\{1, 2, \ldots, N \right\}, ~~ N = \sum_{i=1}^k n_i},
#' then the test statistics for all-pairs comparisons
#' and a balanced design is calculated as
#' \deqn{
#'  \hat{T}_{ij} = \max_{i \le m < m' \le j}
#'  \frac{\left(\bar{R}_{m'} - \bar{R}_m \right)}
#'  {\sigma_a / \sqrt{n}},
#' }{%
#'  SEE PDF
#' }
#'
#' with \eqn{n = n_i; ~ N = \sum_i^k n_i ~~ (1 \le i \le k)}, \eqn{\bar{R}_i}
#' the mean rank for the \eqn{i}th group,
#' and the expected variance (without ties) \eqn{\sigma_a^2 = N \left(N + 1 \right) / 12}.
#'
#' For the NPY'-test (\code{dist = "h"}), if \eqn{T_{ij} > h_{k-1,\alpha,\infty}}.
#'
#' For the unbalanced case with moderate imbalance the test statistic is
#' \deqn{
#'  \hat{T}_{ij} = \max_{i \le m < m' \le j} \frac{\left(\bar{R}_{m'} - \bar{R}_m \right)}
#'  {\sigma_a \left(1/n_m + 1/n_{m'}\right)^{1/2}},
#' }{%
#'  SEE PDF
#' }
#'
#' For the NPY'-test (\code{dist="h"}) the null hypothesis is rejected in an unbalanced design,
#' if \eqn{\hat{T}_{ij} > h_{k,\alpha,\infty} / \sqrt{2}}.
#' In case of a NPY'-test, the function does not return p-values. Instead the critical h-values
#' as given in the tables of Hayter (1990) for \eqn{\alpha = 0.05} (one-sided)
#' are looked up according to the number of groups (\eqn{k-1}) and
#' the degree of freedoms (\eqn{v = \infty}).
#'
#' For the NPT'-test (\code{dist = "Normal"}), the null hypothesis is rejected, if
#' \eqn{T_{ij} > \sqrt{2} t_{\alpha,\infty} = \sqrt{2} z_\alpha}. Although Nashimoto and Wright (2005) originally did not use any p-adjustment,
#' any method as available by \code{\link{p.adjust.methods}} can
#' be selected for the adjustment of p-values estimated from
#' the standard normal distribution.
#'
#' @return
#' Either a list of class \code{"osrt"} if \code{dist = "h"}  or a list
#' of class \code{"PMCMR"} if \code{dist = "Normal"}.
#' @template returnOsrt
#' @template class-PMCMR
#'
#' @note
#' The function will give a warning for the unbalanced case and returns the
#' critical value \eqn{h_{k-1,\alpha,\infty} / \sqrt{2}} if applicable.
#'
#' @references
#' Hayter, A. J.(1990) A One-Sided Studentised Range
#' Test for Testing Against a Simple Ordered Alternative,
#' \emph{J Amer Stat Assoc} \bold{85}, 778--785.
#'
#' Nashimoto, K., Wright, F.T. (2007)
#' Nonparametric Multiple-Comparison Methods for Simply Ordered Medians.
#' \emph{Comput Stat Data Anal} \bold{51}, 5068--5076.
#'
#' @keywords htest nonparametric
#'
#' @seealso
#'  \code{\link{Normal}}, \code{\link{chackoTest}},
#'  \code{\link{NPMTest}}
#'
#' @example examples/shirleyEx.R
# @examples
# ## Example from Shirley (1977)
# ## Reaction times of mice to stimuli to their tails.
# y <- c(2.4, 3, 3, 2.2, 2.2, 2.2, 2.2, 2.8, 2, 3,
# 2.8, 2.2, 3.8, 9.4, 8.4, 3, 3.2, 4.4, 3.2, 7.4, 9.8, 3.2, 5.8,
# 7.8, 2.6, 2.2, 6.2, 9.4, 7.8, 3.4, 7, 9.8, 9.4, 8.8, 8.8, 3.4,
# 9, 8.4, 2.4, 7.8)
# g <- gl(4, 10)
#
# chackoTest(x , g)
#
# ## Default is standard normal distribution (NPT'-test)
# summary(chaAllPairsNashimotoTest(x, g, p.adjust.method = "none"))
#
# ## same but h-distribution (NPY'-test)
# chaAllPairsNashimotoTest(x, g, dist = "h")
#
#' @export
chaAllPairsNashimotoTest <- function(x, ...)
    UseMethod("chaAllPairsNashimotoTest")

#' @rdname chaAllPairsNashimotoTest
#' @aliases chaAllPairsNashimotoTest.default
#' @method chaAllPairsNashimotoTest default
#' @template one-way-parms
#' @param p.adjust.method method for adjusting p values. Ignored if \code{dist = "h"}.
#' @param alternative the alternative hypothesis. Defaults to \code{greater}.
#' @param dist the test distribution. Defaults to \code{Normal}.
#' @importFrom stats pnorm
#' @importFrom stats complete.cases
#' @importFrom stats p.adjust
#' @importFrom stats p.adjust.methods
#' @export
chaAllPairsNashimotoTest.default <-
function(x, g, p.adjust.method = c(p.adjust.methods),
         alternative = c("greater", "less"),
         dist = c("Normal", "h"), ...){
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

    p.adjust.method <- match.arg(p.adjust.method)
    alternative <- match.arg(alternative)
    dist <- match.arg(dist)

    if (alternative == "less") {
        x <- -x
    }

    rij <- rank(x)
    Ri <- tapply(rij, g, mean)
    ni <- tapply(x, g, length)
    k <- nlevels(g)
    N <- length(x)
    df <- Inf

    sigma <- sqrt(N * (N + 1) / 12)

    ## balanced design
    n <- ni[1]
    ## check for all equal
    ok <- sapply(2:k, function(i) ni[i] == n)
    is.balanced <- all(ok)

    STAT <- matrix(NA, ncol=k-1, nrow=k-1)

    if (is.balanced) {
    for (i in 1:(k-1)){
        for(j in (i+1):k){
            u <- j
            m <- i:(u-1)
            tmp <- sapply(m, function(m) {
                    (Ri[u] - Ri[m]) /
                        (sigma / sqrt(n))
            })
            STAT[j-1,i] <- max(tmp)
        }
    }
    } else {
        for (i in 1:(k-1)){
            for(j in (i+1):k){
                u <- j
                m <- i:(u-1)
                tmp <- sapply(m, function(m) {
                    (Ri[u] - Ri[m]) /
                        (sigma * sqrt(1/ni[m] + 1/ni[u]))
                })
                STAT[j-1,i] <- max(tmp)
            }
        }
    }
    colnames(STAT) <- levels(g)[1:(k-1)]
    rownames(STAT) <- levels(g)[2:k]


    if (dist == "Normal") {
        ## modify STAT
        SQRT2 <- sqrt(2)
        STAT <- STAT / SQRT2
        PVAL <- pnorm(STAT, lower.tail = FALSE)
        DIST <- "z"
        METHOD <- "Pairwise comparisons using Nashimoto-Wright's NPT'-Test"
        p <- as.vector(PVAL)
        pad <-
            p.adjust(p, method = p.adjust.method, n = k * (k - 1) / 2)
        PVAL <- matrix(pad, ncol = (k - 1), nrow = (k - 1))

        colnames(PVAL) <- colnames(STAT)
        rownames(PVAL) <- rownames(STAT)
        MODEL <- data.frame(x, g)
        ans <- list(
            method = METHOD,
            data.name = DNAME,
            p.value = PVAL,
            statistic = STAT,
            p.adjust.method = p.adjust.method,
            model = MODEL,
            dist = DIST,
            alternative = "greater"
        )
        class(ans) <- "PMCMR"
        return(ans)

    } else {
        ## h value
        ## get critical h-value with k = k and v = Inf
        nrows <- nrow(TabCrit$hayter.h005)
        kk <- as.numeric(colnames(TabCrit$hayter.h005))

        ##
        k <- k - 1
        ## check for kk
        if (k > max(kk) | k < min(kk)) stop("No critical values for k = ", k)

        hCrit <- unlist(TabCrit$hayter.h005[nrows, paste0(k)])
        hCrit <- adjust.hCrit(hCrit, is.balanced)

        METHOD <- "Pairwise comparisons using Nashimoto-Wright's NPY'-Test"
        parameter = c(k, Inf)
        names(parameter) <- c("k", "df")

        ans <- list(
            method = METHOD,
            data.name = DNAME,
            crit.value = hCrit,
            statistic = STAT,
            parameter = parameter,
            alternative = alternative,
            dist = "h"
        )
        class(ans) <- "osrt"
        return(ans)
    }
}

#' @rdname chaAllPairsNashimotoTest
#' @aliases chaAllPairsNashimotoTest.formula
#' @method chaAllPairsNashimotoTest formula
#' @template one-way-formula
#' @export
chaAllPairsNashimotoTest.formula <-
function(formula, data, subset, na.action,
         p.adjust.method = c(p.adjust.methods),
         alternative = c("greater", "less"),
         dist = c("Normal", "h"),
         ...)
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
    dist <- match.arg(dist)
    names(mf) <- NULL
    y <- do.call("chaAllPairsNashimotoTest",
                 c(as.list(mf),
                   p.adjust.method = p.adjust.method,
                   alternative = alternative,
                   dist = dist))
    y$data.name <- DNAME
    y
}
