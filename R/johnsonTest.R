##  johnsonTest.R
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
##
#' @name johnsonTest
#' @title Testing against Ordered Alternatives (Johnson-Mehrotra Test)
#' 
#' @description
#' Performs the Johnson-Mehrotra test for testing against ordered alternatives
#' in a balanced one-factorial sampling design.
#' 
#' @details
#' The null hypothesis, H\eqn{_0: \theta_1 = \theta_2 = \ldots = \theta_k}
#' is tested against a simple order hypothesis,
#' H\eqn{_\mathrm{A}: \theta_1 \le \theta_2 \le \ldots \le
#' \theta_k,~\theta_1 < \theta_k}.
#'
#' The p-values are estimated from the standard normal distribution.
#' 
#' @template class-htest
#' @template trendTests
#' @references
#' Bortz, J. (1993). \emph{Statistik für Sozialwissenschaftler} (4th ed.).
#' Berlin: Springer.
#' 
#' Johnson, R. A., & Mehrotra, K. G. (1972). Some c-sample
#' nonparametric tests for ordered alternatives.
#' \emph{Journal of the Indian Statistical Association}, 9, 8--23.
#'
#' @export johnsonTest
johnsonTest <- function(x, ...) UseMethod("johnsonTest")

#' @rdname johnsonTest
#' @method johnsonTest default
#' @aliases johnsonTest.default
#' @template one-way-parms
#' @param alternative the alternative hypothesis. Defaults to \code{"two.sided"}.
#' @importFrom stats pnorm complete.cases
#' @importFrom SuppDists normOrder
#' @export
johnsonTest.default <-
    function(x, g, alternative = c("two.sided", "greater", "less"), 
             ...)
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
        if(!is.null(x$alternative)){
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
    alternative <- match.arg(alternative)
    N <- length(x)
    if (N < 2)
        stop("not enough observations")
    n <- tapply(x , g, length)
    if (any (n != N / k))
        stop("Johnson test is only valid for 'k' equal sample sizes")
    n <- unique(n)
    r <- rank(x)
    
    ## transform to expected normal order scores
    sco <- normOrder(N)
    lo <- floor(r)
    up <- ceiling(r)
	
    ## mid-expected normal order scores for mid-ranks
    zscores <- (sco[up] + sco[lo]) / 2

    ## mean scores per group
    etai <- tapply(zscores, g, mean)

    aux <- 1:k
    Mi <- sqrt((aux - 1) * ( 1 - (aux - 1) / k)) -
        sqrt(aux * (1 - aux / k))

    ## T value
    T <- sum(Mi * etai)

    ## test value
    PSTAT <- T / sqrt(sum(zscores^2) / (n * (N - 1)) * sum(Mi^2))
	
    ## get p.value
    if (alternative == "two.sided"){
        PVAL <- 2 * min(0.5, pnorm(abs(PSTAT), lower.tail = FALSE))
    } else if (alternative == "greater") {
        PVAL <- pnorm(PSTAT, lower.tail = FALSE)
    } else {
        PVAL <- pnorm(PSTAT)
    }
	
    names(PSTAT) <- "z"
    names(T) <- "T*"
    METHOD <- paste("Johnson-Mehrotra test")  
    ans <- list(method = METHOD, p.value = PVAL,
                statistic = PSTAT, estimate = T, 
                data.name = DNAME, alternative = alternative)
    class(ans) <- "htest"
    return(ans)
}

#' @rdname johnsonTest
#' @method johnsonTest formula
#' @aliases johnsonTest.formula
#' @template one-way-formula
#' @export
johnsonTest.formula <-
    function(formula, data, subset, na.action,
         alternative = c("two.sided", "greater", "less"), ...)
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
    names(mf) <- NULL
    y <- do.call("johnsonTest", c(as.list(mf), alternative = alternative))
    y$data.name <- DNAME
    y
}
