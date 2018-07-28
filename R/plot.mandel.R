## plot.mandel
## Part of the R package: PMCMRplut
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
##
#' @title Plotting mandel Objects
#' @description
#' Plotting method for objects inheriting from class \code{"mandel"}.
#'
#' @method plot mandel
#' @aliases plot.mandel
#' @param x an object with class \code{"mandel"}.
#' @param alpha level of significance. Defaults to \code{0.005}.
#' @param \ldots further arguments, currently ignored.
#'
#'
#' @importFrom graphics abline barplot par text title
#' @seealso
#' \code{demo(Pentosan)}
#' @examples
#' ##
#' \dontrun{
#' data(Pentosan)
#' md <- mandelkTest(value ~ lab, Pentosan, subset = (material == "B"))
#' plot(md)
#' }
#' @importFrom graphics text
#' @importFrom graphics barplot
#' @importFrom graphics abline
#' @export
plot.mandel <- function(x, alpha = 0.005, ...)
{

    xnm <- deparse(substitute(x))
    if(!inherits(x, "mandel")) stop(xnm, "is not an object of type mandel")

    m <- match.call(expand.dots=TRUE)
    m <- as.list(m)

    ## delete not needed stuff from m
    m[[1]] <- NULL
    m$x <- NULL
    m$alpha <- NULL

    ## Check, whether h or k test
    is.hTest <- grepl(pattern="h", x$method)
    statistics <- x$statistics

    ## critical h or k value
    if (is.hTest) {
        k <- length(x$grouplev)
        crit <- qmandelh(alpha / 2, k, lower.tail = FALSE)
        ylab <- "h - statistics"
        stats <- "h"
    } else {
        k <- length(x$grouplev)
        n <- x$nrofrepl
        crit <- qmandelk(alpha, k, n, lower.tail = FALSE)
        ylab <- "k - statistics"
        stats <- "k"
    }

    ylmx <- 1.2 * max(abs(statistics), crit)
    ylmn <-  ifelse(is.hTest, -ylmx, 0)

    names.arg <- x$grouplev
    main <- x$method
    inp <- list(height = statistics,
                ylim = c(ylmn, ylmx),
                ylab = ylab,
                names.arg = names.arg,
                main = main)
    inp <- c(inp, m)

    ## plot
    mp <- do.call("barplot", inp)

    abline(a = crit, b = 0, lty = 3)
    if (is.hTest) abline(a = -crit, b = 0, lty = 3)

    txt <- paste0(stats, "-crit, p = ", format.pval(alpha))
    text(x = mean(mp),
         y = crit,
         labels = txt,
         pos = 3,
         offset = 0.5)

    invisible(x)
}
