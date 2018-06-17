##  print.williams.R
##
##  Copyright (C) 2015-2017 Thorsten Pohlert
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

#' @title Summarize an williams Object
#' @description
#' Summarize an object of class \emph{williams}.
#' @method summary williams
#' @aliases summary.williams
#' @param object an object of class \code{"williams"}.
#' @param \dots further arguments. Currenly ignored.
#' @keywords models
#' @return
#' A detailed output of all pairwise hypotheses,
#' the test statistics, the corresponding p-values and
#' symbols that indicates the level of significance.
#' @examples
#' ans <- vanWaerdenAllPairsTest(count ~ spray, InsectSprays)
#' summary(ans)
#' @seealso
#' \code{\link{print.williams}}, \code{\link{summaryGroup}}.
#' @importFrom stats symnum
#' @export
summary.williams <- function(object, ...)
{
    OK <- inherits(object, c("williams"))
    if (!OK)
        stop ("Not an object of class williams")
    if (!is.matrix(object$statistic))
        stop ("Matrix object$statistic not found.")
    tval <- as.numeric(object$t.value)
    stat <- as.numeric(object$statistic)
    grp1 <- as.numeric(c(col(object$t.value)))
    cnam <- colnames(object$t.value)
    grp2 <- as.numeric(c(row(object$t.value)))
    rnam <- rownames(object$t.value)
    STAT <- object$dist

    if (object$alternative == "greater") {
      dec <- ifelse(stat > tval, "reject", "accept")
    } else {
      dec <- ifelse(stat < tval, "reject", "accept")
    }

    TVAL <- "t\' crit"
    H0 <- paste(rnam[grp2], "-", cnam[grp1], "==", "0")
    STAT2 <- paste0(STAT, " value")

    #OK <- !is.na(pval)
    message("\n\tMany-to-one comparisons using ", object$method, "\n")
    message("data: ", object$data.name)
    if (!is.null(object$alternative)){
        message("alternative hypothesis: ", object$alternative)
    }
    message("degree of freedom: ", object$df)
    message("alpha: 0.05\n")

    xdf <- data.frame(statistic = round(stat[OK], 3),
                      t.value = round(tval[OK], 3),
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

#' @title williams Printing
#' @description
#' \code{print.williams} is the \emph{williams} method of the generic
#' \code{\link{print}} function which prints its argument
#' and returns it \emph{invisibly} (via \code{\link{invisible}(x)}).
#' @param x an object used to select a method.
#' @param \ldots further arguments. Currently ignored.
#' @aliases print.williams
#' @method print williams
#' @keywords print
#' @export
print.williams <-
function(x, ...)
{
  summary.williams(x)
    invisible(x)
}
