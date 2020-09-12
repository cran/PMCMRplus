##  print.steel.R
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

#' @title Summarize a steel Object
#' @description
#' Summarize an object of class \emph{steel}.
#' @method summary steel
#' @aliases summary.steel
#' @param object an object of class \code{"steel"}.
#' @param \dots further arguments. Currenly ignored.
#' @keywords methods
#' @return
#' A detailed output of all pairwise hypotheses,
#' the test statistics, the corresponding p-values and
#' symbols that indicates the level of significance.
#' @examples
#' ans <- vanWaerdenAllPairsTest(count ~ spray, InsectSprays)
#' summary(ans)
#' @seealso
#' \code{\link{print.steel}}, \code{\link{summaryGroup}}.
#' @importFrom stats symnum
#' @export
summary.steel <- function(object, ...)
{
    OK <- inherits(object, c("steel"))
    if (!OK)
        stop ("Not an object of class steel")
    if (!is.matrix(object$statistic))
        stop ("Matrix object$statistic not found.")
    tval <- as.numeric(object$R.crit)
    stat <- as.numeric(object$statistic)
    grp1 <- as.numeric(c(col(object$R.crit)))
    cnam <- colnames(object$R.crit)
    grp2 <- as.numeric(c(row(object$R.crit)))
    rnam <- rownames(object$R.crit)
    STAT <- object$dist

 #   if (object$alternative == "greater") {
#      dec <- ifelse(stat > tval, "reject", "accept")
#    } else {
      dec <- ifelse(stat <= tval, "reject", "accept")
 #   }

    TVAL <- "R crit"
    H0 <- paste0(rnam[grp2], " = ", cnam[grp1])
    STAT2 <- paste0(STAT)

    k <- nlevels(object$model$g) - 1 # ex control
    n <- nrow(object$model) / (k + 1) # incl control

    #OK <- !is.na(pval)
    message("\n\t", object$method, "\n")
    message("data: ", object$data.name)
    if (!is.null(object$alternative)){
        message("alternative hypothesis: ", object$alternative)
    }
    message("alpha: 0.05")
    message("k: ", k)
    message("n: ", n)
    xdf <- data.frame(statistic = round(stat[OK], 3),
                      Rcrit = round(tval[OK], 3),
                      dec = dec)
     #                 symp)
    rownames(xdf) <- H0[OK]
    names(xdf) <- c(STAT2, TVAL, "decision")
    message("H0")
    ##
    print(xdf)
    message("---")
    #message("Signif. codes: ", attr(symp, 'legend'))
    invisible(object)
}

#' @title steel Printing
#' @description
#' \code{print.steel} is the \emph{steel} method of the generic
#' \code{\link{print}} function which prints its argument
#' and returns it \emph{invisibly} (via \code{\link{invisible}(x)}).
#' @param x an object used to select a method.
#' @param \ldots further arguments. Currently ignored.
#' @aliases print.steel
#' @method print steel
#' @keywords print
#' @export
print.steel <-
function(x, ...)
{
  summary.steel(x)
    invisible(x)
}
