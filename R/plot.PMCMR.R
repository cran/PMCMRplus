##  plot.PMCMR.R
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

#' @title Plotting PMCMR Objects
#' @description
#' Plots a bar-plot for objects of class \code{"PMCMR"}.
#' @aliases barPlot
#' @param x an object of class \code{"PMCMR"}.
#' @param alpha the selected alpha-level. Defaults to 0.05.
#' @param \dots further arguments for method \code{barplot}.
#' @return
#' A barplot where the height of the bars corresponds to the arithmetic
#' mean. The extend of the whiskers are \eqn{\pm z_{(1-\alpha/2)}
#' \times s_{\mathrm{E},i}}, where the latter denotes the standard error
#' of the \eqn{i}th group. Symbolic letters are depicted on top of the bars,
#' whereas different letters indicate significant differences between
#' groups for the selected level of alpha.
#' @note
#' The barplot is strictly spoken only valid for normal data, as
#' the depicted significance intervall implies symetry.
#' @importFrom multcompView multcompLetters
#' @importFrom graphics text
#' @importFrom graphics barplot
#' @importFrom graphics segments
#' @importFrom stats sd
#' @importFrom stats qnorm
#' @keywords plot
#' @examples
#' ## data set chickwts
#' ans <- tukeyTest(weight ~ feed, data = chickwts)
#' barPlot(ans)
#' @export
barPlot <-
function(x, alpha = 0.05, ...)
{
		 
    OK <- inherits(x, c("PMCMR"))
    if (!OK)
        stop ("Not an x of class PMCMR or pairwise.htest")
    if (!is.matrix(x$p.value))
        stop ("Matrix x$p.value not found.")
    ppval <- get.pvalues(x)
    out.mcv <- multcompLetters(ppval, threshold = alpha)$Letters
    dat <- x$model
    groups <- dat$g
    levels(groups) <-  c(colnames(x$statistic)[1], rownames(x$statistic))
    resp <- dat$x
    
    m <- match.call(expand.dots=TRUE)
    m <- as.list(m)
	
    ## delete not needed stuff from m
    m[[1]] <- NULL
    m$x <- NULL
    m$type <- NULL
    m$alpha <- NULL
    
    labs <- strsplit(x$data.name, " ")[[1]]
 
    k <- nlevels(groups)
    n <- tapply(resp, groups, length)
    xmean <- tapply(resp, groups, mean)
    xse <- tapply(resp, groups, sd) / sqrt(n)

    qn <- 2 * qnorm(1 - alpha)
    up <- xmean + qn * xse
    lo <- xmean - qn * xse

    ## put everything into inp
    inp <- list(height = xmean)

    if(is.null(m$ylim)){ 
        inp$ylim <- c(#1.3 * abs(min(resp)) * sign(min(resp)),
                                        #min(resp),
            ifelse(min(resp) < 0, min(resp), 0),
            1.3 * abs(max(up)) * sign(max(up)))
    }
    inp <- c(inp, m)
    mp <- do.call("barplot", inp)
	
    text(mp, 1.2 * max(up), out.mcv)

    eps <- 0.02
    segments(mp, up, mp, lo)
    segments(mp-eps, lo, mp+eps, lo)
    segments(mp-eps, up, mp+eps, up)
}	

#' @title Plotting PMCMR Objects
#' @description
#' Plotting method for objects inheriting from class \code{"PMCMR"}.
#' @aliases plot.PMCMR
#' @method plot PMCMR
#' @param x an object of class \code{"PMCMR"}.
#' @param alpha the selected alpha-level. Defaults to 0.05.
#' @param \dots further arguments for method \code{boxplot}.
#' @return
#' A box-whisker plot for each factor level. The range of the whiskers indicate
#' the extremes (\code{boxplot = x, ..., range=0}). Letter symbols are depicted on top of each box.
#' Different letters indicate significant
#' differences between groups on the selected level of alpha.
#' @importFrom graphics boxplot
#' @importFrom graphics text
#' @importFrom graphics par
#' @importFrom multcompView multcompLetters
#' @keywords plot
#' @examples
#' ## data set InsectSprays
#' ans <- kwAllPairsNemenyiTest(count ~ spray, data = InsectSprays)
#' plot(ans)
#' plot(ans, col="red",main="My title", xlab="Spray", "Count")
#' @seealso
#' \code{\link{boxplot}}
#' @export
plot.PMCMR <- function(x, alpha = 0.05, ...)
{
    OK <- inherits(x, c("PMCMR"))
    if (!OK)
        stop ("Not an x of class PMCMR")
    if (!is.matrix(x$p.value))
        stop ("Matrix x$p.value not found.")
    ppval <- get.pvalues(x)
    out.mcv <- multcompLetters(ppval, threshold = alpha)$Letters
    dat <- x$model
    groups <- dat$g
    levels(groups) <-  c(colnames(x$statistic)[1], rownames(x$statistic))
    resp <- dat$x
    m <- match.call(expand.dots=TRUE)
    m <- as.list(m)
	
    ## delete not needed stuff from m
    m[[1]] <- NULL
    m$x <- NULL
    m$type <- NULL
    m$alpha <- NULL
    
    labs <- strsplit(x$data.name, " ")[[1]]
    
##    if (type == "boxplot"){
        ## get the data
    inbxp <- split(resp, groups)
    names(inbxp) <- NULL
	
    ## check for parameters given for boxplot
    if (is.null(m$range)) inbxp$range <- 0
    if (is.null(m$width)) inbxp$width <- NULL
    if (is.null(m$varwidth)) inbxp$varwidth <- FALSE
    if (is.null(m$notch)) inbxp$notch <- FALSE
    if (is.null(m$outline)) inbxp$outline <- TRUE
    if (is.null(m$names)) inbxp$names <- levels(groups)
    if (is.null(m$plot)) inbxp$plot <- TRUE
    if (is.null(m$border)) inbxp$border <- par("fg")
    if (is.null(m$col)) inbxp$col <- NULL
    if (is.null(m$log)) inbxp$log <- ""
    if (is.null(m$horizontal)) inbxp$horizontal <- FALSE
    if (is.null(m$add)) inbxp$add <- NULL
    
    ## check for graphic parameters	   
    pp <- list(NULL)
    if (is.null(m$xlab)){
        pp$xlab <- labs[3]
    } else {
        pp$xlab <- m$xlab
        m$xlab <- NULL
    }
    if (is.null(m$ylab)){
        pp$ylab <- labs[1]
    } else {
        pp$ylab <- m$ylab
            m$ylab <- NULL
    }
    if (is.null(m$main)){
        pp$main <- x$method
    } else {
        pp$main <- m$main
        m$main <- NULL
    }
    if (is.null(m$ylim)){
        pp$ylim <- c(#1.3 * abs(min(resp)) * sign(min(resp)),
            min(resp),
            1.3 * abs(max(resp)) * sign(max(resp)))
    } else {
        pp$ylim <- m$ylim
        m$ylim <- NULL
    }
    if (is.null(m$pars)){
        inbxp <- c(inbxp, m)
        inbxp$pars <- pp
    } else {
        ppp <- as.list(m$pars)
        inbxp <- c(inbxp, m)
        inbxp$pars <- c(pp, ppp)
    }
	    
    (mp <- do.call("boxplot", inbxp))
    k <- nlevels(groups)
    text(rep(1.2 * max(resp), k), out.mcv)    
} 
