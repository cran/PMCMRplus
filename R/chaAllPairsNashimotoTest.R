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

#' @title All-Pairs Comparisons for Simply Ordered Mean Ranksums
#'
#' @description
#' Performs Nashimoto and Wright's all-pairs comparison procedure
#' for simply ordered mean ranksums (NPY-test).
#' According to the authors, the procedure shall only be
#' applied after Chacko's test (see \code{\link{chackoTest}}) indicates
#' global significance.
#'
#' The modified procedure uses the property of a simple order,
#' \eqn{\theta_m' - \theta_m \le \theta_j - \theta_i \le \theta_l' - \theta_l
#' \qquad (l \le i \le m~\mathrm{and}~ m' \le j \le l')}.
#' The null hypothesis H\eqn{_{ij}: \theta_i = \theta_j} is tested against
#' the alternative A\eqn{_{ij}: \theta_i < \theta_j} for any
#' \eqn{1 \le i < j \le k}.
#'
#' In the NPT test the p-values
#' are estimated from the standard normal distribution.
#'
#' @details
#' Although Nashimoto and Wright (2005) originally did not use any p-adjustment,
#' any method as available by \code{\link{p.adjust.methods}} can
#' be selected for the adjustment of p-values estimated from
#' the standard normal distribution.
#'
#' @name chaAllPairsNashimotoTest
#' @references
#' Nashimoto, K., Wright, F.T. (2007)
#'  Nonparametric Multiple-Comparison Methods for Simply Ordered Medians.
#'  \emph{Comput Stat Data Anal} \bold{51}, 5068–5076.
#'
#' @keywords htest nonparametric
#'
#' @seealso
#'  \code{\link{Normal}}, \code{\link{chackoTest}}
#' @template class-PMCMR
#' @examples
#' ## Example from Sachs (1997, p. 402)
#' x <- c(106, 114, 116, 127, 145,
#'        110, 125, 143, 148, 151,
#'        136, 139, 149, 160, 174)
#' g <- gl(3,5)
#' levels(g) <- c("A", "B", "C")
#' chackoTest(x , g)
#' chaAllPairsNashimotoTest(x, g, p.adjust.method = "none")
#' @export
chaAllPairsNashimotoTest <- function(x, ...)
    UseMethod("chaAllPairsNashimotoTest")

#' @rdname chaAllPairsNashimotoTest
#' @aliases chaAllPairsNashimotoTest.default
#' @method chaAllPairsNashimotoTest default
#' @template one-way-parms
#' @param p.adjust.method method for adjusting p values
#' @importFrom stats pnorm
#' @importFrom stats ptukey
#' @importFrom stats complete.cases
#' @importFrom stats p.adjust
#' @importFrom stats p.adjust.methods
#' @export
chaAllPairsNashimotoTest.default <-
function(x, g, p.adjust.method = c(p.adjust.methods), ...){
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

    rij <- rank(x)
    Ri <- tapply(rij, g, mean)
    ni <- tapply(x, g, length)
    k <- nlevels(g)
    N <- length(x)
    df <- Inf

    sigma <- sqrt(N * (N + 1) / 12)


    n <- ni[1]
    ## check for all equal
    ok <- sapply(2:k, function(i) ni[i] == n)
    if (!all(ok)) {
        warning("NPM-test is for balanced designs only. Using n = Mean(ni).")
        n <- round(mean(ni), 0)
    }


    STAT <- matrix(NA, ncol=k-1, nrow=k-1)
    for (i in 1:(k-1)){
        for(j in (i+1):k){
            u <- j
            m <- i:(u-1)
            tmp <- sapply(m, function(m) {
      #          if (p.adjust.method != "single step"){
      #              (Ri[u] - Ri[m]) / (sqrt(2) * sigma *
      #                                 sqrt(1 / ni[m] + 1 /ni[u]))
      #          } else {
                    (Ri[u] - Ri[m]) /
                        (sigma * sqrt(2) / sqrt(n))
       #         }
            })
            STAT[j-1,i] <- max(tmp)
        }
    }

    colnames(STAT) <- levels(g)[1:(k-1)]
    rownames(STAT) <- levels(g)[2:k]

#    if (p.adjust.method == "single-step"){
#        PVAL <- ptukey(STAT, nmeans = (k-1),
#                       df = df, lower.tail=FALSE)
#        DIST <- "q"
#        METHOD <- "Nashimoto-Wright NPY'-Test for ordered means \n\t\t of non-normal data"
#    } else {
        PVAL <- pnorm(STAT, lower.tail=FALSE)
        DIST <- "z"
        METHOD <- "Nashimoto-Wright NPT-Test for ordered means \n\t\t of non-normal data"
        p <- as.vector(PVAL)
        pad <- p.adjust(p, method = p.adjust.method, n = k * (k - 1) / 2)
        PVAL <- matrix(pad, ncol=(k-1), nrow=(k-1))
#    }
    colnames(PVAL) <- colnames(STAT)
    rownames(PVAL) <- rownames(STAT)
    MODEL <- data.frame(x, g)
    ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL,
                statistic = STAT, p.adjust.method = p.adjust.method,
                model = MODEL, dist = DIST, alternative = "greater")
    class(ans) <- "PMCMR"
    ans
}

#' @rdname chaAllPairsNashimotoTest
#' @aliases chaAllPairsNashimotoTest.formula
#' @method chaAllPairsNashimotoTest formula
#' @template one-way-formula
#' @export
chaAllPairsNashimotoTest.formula <-
function(formula, data, subset, na.action,
         p.adjust.method = c(p.adjust.methods), ...)
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
    y <- do.call("chaAllPairsNashimotoTest", c(as.list(mf),
                                     p.adjust.method))
    y$data.name <- DNAME
    y
}
