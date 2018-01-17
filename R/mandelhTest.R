## mandelhTestR
## Part of the R package: PMCMRplus
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
##
#' @name mandelhTest
#' @title Mandel's h test according to E 651 ASTM
#' @description The function calculates the 
#'  consistency statistics h and corresponding
#'  p-values for each group (lab) according to
#' Practice E 691 ASTM.
#'
#' @template class-mandel
#' 
#' @seealso
#' \code{\link{qmandelh}} \code{\link{pmandelh}} 
#' @references
#' 
#' Practice E 691, 2005, \emph{Standard Practice for
#' Conducting an Interlaboratory Study to Determine the
#' Precision of a Test Method}, ASTM International.
#' 
#' @importFrom stats pt qt complete.cases sd
#' @keywords htest
#' @examples
#' data(Pentosan)
#' mandelhTest(value ~ lab, data=Pentosan, subset=(material == "A"))
#' @export 
mandelhTest <- function(x, ...) UseMethod("mandelhTest")

#' @rdname mandelhTest
#' @method mandelhTest default
#' @aliases mandelhTest.default
#' @template one-way-parms
#' @export 
mandelhTest.default <- function(x, g, ...)
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
	
    ## Cell averages
    xbar <- tapply(x, g, mean)

    ## Cell standard deviation
    s <- tapply(x, g, sd)
    ## mean for unbalanced designs
    n <- mean(tapply(x, g, length))

    xgrandmean <- mean(x)
    d <- sapply(xbar, function(i) i - xgrandmean)

    ## standard deviation of cell averages
    sxfn <- function(x){
        p <- length(x)
        grandmean <- mean(x)
        d <- sapply(x, function(i) i - grandmean)
        sx <- sqrt(sum(d^2/(p-1)))
        return(sx)
    }

    sx <- sxfn(xbar)
    ## h-statistic
    h <- d /sx

    ## two.sided upper quantile
    pval <- sapply(h, function(h)
        2 * min(0.5, pmandelh(q=abs(h), k=k, lower.tail=FALSE)))
    
    METHOD <- "Mandel's h-test"

    ans <- list(method = METHOD,
                data.name = DNAME,
	        p.value = pval,
                statistics = h,
                grouplev = levels(g),
                nrofrepl = n,
                alternative = "two.sided")
    class(ans) <- "mandel"
    ans
}

#' @rdname mandelhTest
#' @method mandelhTest formula
#' @aliases mandelhTest.formula
#' @template one-way-formula
#' @export
mandelhTest.formula <-
    function(formula, data, subset, na.action, ...)
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
    y <- do.call("mandelhTest", c(as.list(mf)))
    y$data.name <- DNAME
    y
}
