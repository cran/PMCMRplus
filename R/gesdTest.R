## gesdTest.R
## Part of the R package: PMCMRplus
##
## Copyright (C) 2017 Thorsten Pohlert
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
#' @title Generalized Extreme Studentized Deviate Many-Outlier Test
#' @description
#' Performs Rosner's generalized extreme studentized deviate
#' procedure to detect up-to \code{maxr} outliers in a
#' univariate sample that follows an approximately normal distribution.
#' 
#' @references
#' Rosner, B. (1983), Percentage Points for a Generalized ESD
#' Many-Outlier Procedure, \emph{Technometrics}, 25, 165--172.
#'
#' @examples
#' ## Taken from Rosner (1983):
#' x <- c(-0.25,0.68,0.94,1.15,1.20,1.26,1.26,
#' 1.34,1.38,1.43,1.49,1.49,1.55,1.56,
#' 1.58,1.65,1.69,1.70,1.76,1.77,1.81,
#' 1.91,1.94,1.96,1.99,2.06,2.09,2.10,
#' 2.14,2.15,2.23,2.24,2.26,2.35,2.37,
#' 2.40,2.47,2.54,2.62,2.64,2.90,2.92,
#' 2.92,2.93,3.21,3.26,3.30,3.59,3.68,
#' 4.30,4.64,5.34,5.42,6.01)
#'
#' out <- gesdTest(x, 10)
#'
#' ## print method
#' out
#'
#' ## summary method
#' summary(out)
#' 
#' @importFrom stats qt
#' @importFrom stats sd
#' @keywords htest univariate
#' @concept outliers
#' @param x a numeric vector of data.
#' @param maxr the maximum number of outliers to be tested.
#' @export
gesdTest <- function(x, maxr){

    x <- na.omit(x)
    n <- length(x)
    if (n < 25 & n >=15) {
        warning("Due to sample-size, results are 'reasonable'")
    } else if (n <15) {
        warning("Due to sample-size, results are 'not reasonable'")
    }
        
    if(maxr > n){
        stop("Number of potential outliers > sample-size. Reduce 'maxr'")
    }
    oldx <- x
    ix <- rep(NA, maxr)
    PVAL <- rep(NA, maxr)
    R <- rep(NA, maxr)
    ## repeated single outlier Grubb's test
    for (i in (1:maxr)){
        out <- grubbsTest(x, alternative = "two.sided")
        o <- out$estimate[1]

        ## Danger!! This does not select
        ## the correct number, as x changes
        ix[i] <- which(out$estimate[2] == oldx)
        PVAL[i] <- out$p.value
        R[i] <- out$statistic
        x <- x[-o]
    }

    ans <- list(method = "GESD multiple outlier test",
                statistic = R,
                p.value = PVAL,
                ix = ix,
                alternative = "two.sided")
    class(ans) <- "gesdTest"
    return(ans)
}

#' @title gesdTest Printing
#' @description
#' \code{print.gesdTest} is the \emph{gesdTest} method of the generic
#' \code{\link{print}} function which prints its argument
#' and returns it \emph{invisibly} (via \code{\link{invisible}(x)}).
#' @param x an object used to select a method.
#' @param \ldots further arguments. Currently ignored.
#' @aliases print.gesdTest
#' @method print gesdTest
#' @keywords print
#' @export
print.gesdTest <- function(x, ...)
{
    cat("\n\t", x$method, "\n\n")
    n <- length(x$p.value)
    mat <- matrix(x$p.value, ncol=1, nrow= n,
                  dimnames = list(c(1:n),
                                  c("p-value")))
    cat("Nr. of outliers tested:\n")
    print(mat)
    cat("\nalternative hypothesis: ", x$alternative, "\n\n")
    invisible(x)
}

#' @title Summarize an gesdTest Object
#' @description
#' Summarize an object of class \emph{gesdTest}.
#' @method summary gesdTest
#' @aliases summary.gesdTest
#' @param object an object of class \code{"gesdTest"}.
#' @param \dots further arguments. Currenly ignored.
#' @keywords models
#' @importFrom stats symnum
#' @export
summary.gesdTest <- function(object, ...)
{
    
    cat("\n\t", object$method, "\n\n")
    n <- length(object$p.value)

    symp <- symnum(object$p.value, corr=FALSE,
                   cutpoints = c(0,  .001,.01,.05, .1, 1),
                   symbols = c("***","**","*","."," "))

    outl <- character(n)
    for (i in 1:n){
        tmp <- ""
        for (j in 1:i){
            tmp <- paste0(tmp, " ", object$ix[j])
        }
        outl[i] <- tmp
    }
    
    df <- data.frame(outl, object$statistic,
                     format.pval(object$p.value), symp)
    colnames(df) <- c("i", "R", "Pr(>|R|)", "")
    
    cat("Outliers tested:\n")
    print(df)
    cat("\nalternative hypothesis: ", object$alternative, "\n")
    cat("---")
    cat("\nSignif. codes: ", attr(symp, 'legend'), "\n\n")
    invisible(object)
}
