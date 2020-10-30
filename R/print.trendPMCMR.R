##  print.trendPMCMR.R
##
##  Copyright (C) 2015-2018 Thorsten Pohlert
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
## http://www.r-project.org/Licenses/

#' @title Summarize an trendPMCMR Object
#' @description
#' Summarize an object of class \emph{trendPMCMR}.
#' @method summary trendPMCMR
#' @aliases summary.trendPMCMR
#' @param object an object of class \code{"trendPMCMR"}.
#' @param \dots further arguments. Currenly ignored.
#' @keywords methods
#' @return
#' A detailed output of all pairwise hypotheses,
#' the test statistics, the corresponding p-values and
#' symbols that indicates the level of significance.
#'
#' @seealso
#' \code{\link{print.trendPMCMR}}
#' @importFrom stats symnum
#' @export
summary.trendPMCMR <- function(object, ...)
{
    OK <- inherits(object, c("trendPMCMR"))
    if (!OK)
        stop ("Not an object of class trendPMCMR")
    if (!is.matrix(object$statistic))
        stop ("Matrix object$statistic not found.")
    pval <- as.numeric(object$p.value)
    stat <- as.numeric(object$statistic)
    grp1 <- as.numeric(c(col(object$p.value)))
    cnam <- colnames(object$p.value)
    grp2 <- as.numeric(c(row(object$p.value)))
    rnam <- rownames(object$p.value)
    STAT <- object$dist

    #c <- 0

    #  H0 <- paste(rnam[grp2], "-", cnam[grp1], "==", "0")
    if (!is.null(object$alternative)) {
        PVAL <- switch(
            object$alternative,
            less = paste("Pr(<", STAT, ")", sep = ""),
            greater = paste("Pr(>", STAT, ")", sep = ""),
            two.sided = paste("Pr(>|", STAT, "|)", sep = "")
        )

        if (object$alternative == "less") {
            H0 <- sapply(grp2, function(i) {
                #c <- c + 1
                if (i == 1) {
                    paste0(rnam[i], " >= ", cnam[1])
                } else if (i == 2) {
                    paste0(rnam[i], " >= ", rnam[i-1],
                           " >= ", cnam[1])
                } else {
                    paste0(rnam[i], " >= ", rnam[i-1],
                           " >= ... >= ", cnam[1])
                }
            })
        } else if (object$alternative == "greater") {
            H0 <- sapply(grp2, function(i) {
                #c <- c + 1
                if (i == 1) {
                    paste0(rnam[i], " <= ", cnam[1])
                } else if (i == 2) {
                    paste0(rnam[i], " <= ", rnam[i-1],
                           " <= ", cnam[1])
                } else {
                    paste0(rnam[i], " <= ", rnam[i-1],
                           " <= ... <= ", cnam[1])
                }
            })
        } else {
            H0 <- sapply(grp2, function(i) {
                #c <- c + 1
                if (i == 1) {
                    paste0(rnam[i], " == ", cnam[1])
                } else if (i == 2) {
                    paste0(rnam[i], " == ", rnam[i-1],
                           " == ", cnam[1])
                } else {
                    paste0(rnam[i], " == ", rnam[i-1],
                           " == ... == ", cnam[1])
                }
            })
        }
    } else {
        PVAL <- paste("Pr(>|", STAT, "|)", sep = "")
        H0 <- sapply(grp2, function(i) {
            #c <- c + 1
            if (i == 1) {
                paste0(rnam[i], " == ", cnam[1])
            } else if (i == 2) {
                paste0(rnam[i], " == ", rnam[i-1],
                       " == ", cnam[1])
            } else {
                paste0(rnam[i], " == ", rnam[i-1],
                       " == ... == ", cnam[1])
            }
        })
    }

    STAT2 <- paste0(STAT, " value")
    OK <- !is.na(pval)
    cat("\n\t", object$method, "\n")
    cat("data: ", object$data.name, "\n")
    if (!is.null(object$alternative)) {
        cat("alternative hypothesis: ", object$alternative, "\n")
    }
    cat("P value adjustment method: ", object$p.adjust.method, "\n")
    ## Symbols
    symp <- symnum(
        pval[OK],
        corr = FALSE,
        cutpoints = c(0,  .001, .01, .05, .1, 1),
        symbols = c("***", "**", "*", ".", " ")
    )

    xdf <- data.frame(statistic = round(stat[OK], 3),
                      p.value = format.pval(pval[OK]),
                      symp)
    rownames(xdf) <- H0[OK]
    names(xdf) <- c(STAT2, PVAL, "")
    cat("H0\n")
    ##
    print(xdf)
    cat("---\n")
    cat("Signif. codes: ", attr(symp, 'legend'), "\n")
    invisible(object)
}

#' @title trendPMCMR Printing
#' @description
#' \code{print.trendPMCMR} is the \emph{trendPMCMR} method of the generic
#' \code{\link{print}} function which prints its argument
#' and returns it \emph{invisibly} (via \code{\link{invisible}(x)}).
#' @param x an object used to select a method.
#' @param \ldots further arguments. Currently ignored.
#' @aliases print.trendPMCMR
#' @method print trendPMCMR
#' @keywords print
#' @export
print.trendPMCMR <-
    function(x, ...)
    {
        message("\n\t", x$method, "\n")
        message("data: ", x$data.name, "\n")
        pp <- format.pval(x$p.value, 2, na.form = "-")
        attributes(pp) <- attributes(x$p.value)
        print(pp, quote = FALSE, ...)
        message("\nP value adjustment method: ", x$p.adjust.method)
        if (!is.null(x$alternative)) {
            message("alternative hypothesis: ", x$alternative)
        }
        invisible(x)
    }
