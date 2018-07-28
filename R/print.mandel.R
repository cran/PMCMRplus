## print.mandel.R
## Part of the R package: PMCMRplus
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
##
#' @title Mandel Printing
#' @description
#' \code{print.mandel} is the \emph{mandel} method of the generic
#' \code{\link{print}} function which prints its argument
#' and returns it \emph{invisibly} (via \code{\link{invisible}(x)}).
#' @param x an object used to select a method.
#' @param \ldots further arguments. Currently ignored.
#' @seealso
#' \code{\link{mandelhTest}}, \code{\link{mandelkTest}}
#' @aliases print print.mandel
#' @method print mandel
#' @keywords print
#' @importFrom stats symnum
#' @export
print.mandel <- function(x, ...)
{
    if(!inherits(x, "mandel")){
        stop("'x' is not an object of class 'mandel'")
    }
    if(grepl("h", x$method)){
        P <- "Pr(>|h|)"

    } else {
        P <- "Pr(>k)"
    }
    message("\n\t", x$method)
    message(x$data.name)
    message("Nr. of groups (labs): ", length(x$statistic))
    message("Nr. of replicates: ", x$n)
    symp <- symnum(x$p.value, corr=FALSE,
                   cutpoints = c(0,  .001,.01,.05, .1, 1),
                   symbols = c("***","**","*","."," "))
    out <- data.frame(format.pval(x$p.value), symp)
    rownames(out) <- x$grouplev
    colnames(out) <- c(P,"")
    message("Probabilities")
    print(out)
    message("---")
    message("Signif. codes: ", attr(symp, 'legend'))
    invisible(x)
}

## summary function
#' @title Object Summary for class \code{"mandel"}
#' @description \code{summary.mandel} is a function
#' used to produce result summaries of the results of
#' the functions \code{\link{mandelhTest}} or \code{\link{mandelkTest}}.
#' @param object an object of class \code{"mandel"} for
#' which a summary is desired.
#' @param \ldots further arguments. Currently ignored.
#' @seealso
#' \code{\link{mandelhTest}}, \code{\link{mandelkTest}}
#' @name summary.mandel
#' @aliases summary.mandel
#' @method summary mandel
#' @importFrom stats symnum
#' @keywords methods
#' @export
summary.mandel <- function(object, ...)
{
    if(!inherits(object, "mandel")){
        stop("'object' is not an object of class 'mandel'")
    }

    symp <- symnum(object$p.value, corr=FALSE,
                   cutpoints = c(0,  .001,.01,.05, .1, 1),
                   symbols = c("***","**","*","."," "))

    if(grepl("h-test", object$method)){
        P <- "Pr(>|h|)"
        V <- "h"
        rwn <- sapply(object$grouplev, function(g)
            paste0(g, " = Grand Mean"))
        message("\n\t", object$method)
        message("\tBetween-laboratory consistency")
    } else {
        P <- "Pr(>k)"
        V <- "k"
        rwn <- sapply(object$grouplev, function(g)
            paste0(g, " = 'Repeatability SD'"))
        message("\n\t", object$method)
        message("\tWithin-laboratory consistency")
    }

    message("\nData:")
    message(object$data.name)
    message("\nNr. of groups (labs): ", length(object$statistics))
    message("Nr. of replicates: ", object$nrofrepl, "\n")
    out <- data.frame(formatC(object$statistics, digits=4, format="f"),
                      format.pval(object$p.value), symp)
    rownames(out) <- rwn
    colnames(out) <- c(V, P,"")
    message("Null hypothesis")
    print(out)
    message("---")
    message("Signif. codes: ", attr(symp, 'legend'),"\n")
    invisible(object)
}
