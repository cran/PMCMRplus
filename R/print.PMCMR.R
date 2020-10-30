##  print.PMCMR.R
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

#' @title Summarize an PMCMR Object
#' @description
#' Summarize an object of class \emph{PMCMR}.
#' @method summary PMCMR
#' @aliases summary.PMCMR
#' @param object an object of class \code{"PMCMR"}.
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
#' \code{\link{print.PMCMR}}, \code{\link{summaryGroup}}.
#' @importFrom stats symnum
#' @export
summary.PMCMR <- function(object, ...)
{
    OK <- inherits(object, c("PMCMR"))
    if (!OK)
        stop ("Not an object of class PMCMR")
    if (!is.matrix(object$statistic))
        stop ("Matrix object$statistic not found.")
    pval <- as.numeric(object$p.value)
    stat <- as.numeric(object$statistic)
    grp1 <- as.numeric(c(col(object$p.value)))
    cnam <- colnames(object$p.value)
    grp2 <- as.numeric(c(row(object$p.value)))
    rnam <- rownames(object$p.value)
    STAT <- object$dist

    if (!is.null(object$alternative)) {
        if (object$alternative == "less"){
            H0 <- paste(rnam[grp2], "-", cnam[grp1], ">=", "0")
            PVAL <- paste("Pr(<", STAT, ")", sep="")
        } else if (object$alternative == "greater"){
            H0 <- paste(rnam[grp2], "-", cnam[grp1], "<=", "0")
            PVAL <- paste("Pr(>", STAT, ")", sep="")
        } else {
            H0 <- paste(rnam[grp2], "-", cnam[grp1], "==", "0")
            PVAL <- paste("Pr(>|", STAT, "|)", sep="")
        }
    } else {
        H0 <- paste(rnam[grp2], "-", cnam[grp1], "==", "0")
        PVAL <- paste("Pr(>|", STAT, "|)", sep="")
    }

    STAT2 <- paste0(STAT, " value")
    OK <- !is.na(pval)
    message("\n\tPairwise comparisons using ", object$method, "\n")
    message("data: ", object$data.name)
    if (!is.null(object$alternative)){
        message("alternative hypothesis: ", object$alternative)
    }
    message("P value adjustment method: ", object$p.adjust.method)
    ## Symbols
    symp <- symnum(pval[OK], corr=FALSE,
                   cutpoints = c(0,  .001,.01,.05, .1, 1),
                   symbols = c("***","**","*","."," "))

    xdf <- data.frame(statistic = round(stat[OK], 3),
                      p.value = format.pval(pval[OK]),
                      symp)
    rownames(xdf) <- H0[OK]
    names(xdf) <- c(STAT2, PVAL, "")
    message("H0")
    ##
    print(xdf)
    message("---")
    message("Signif. codes: ", attr(symp, 'legend'))
    invisible(object)
}

#' @title Grouped Summary of an PMCMR Object
#' @description
#' Performes a grouped summary on an PMCMR object.
#' @aliases summaryGroup
#' @param x an object of class \code{"PMCMR"}.
#' @param alpha the selected alpha-level. Defaults to 0.05.
#' @param \dots further arguments. Currently ignored.
#' @return
#' Provides summary statistics for each factor level
#' and a letter symbol, whereas different letters indicate
#' significant differences between factor levels based on the
#' selected level of alpha.
#' @keywords methods
#' @seealso
#' \code{\link{summary.PMCMR}}
#' @importFrom multcompView multcompLetters
#' @importFrom stats quantile
#' @importFrom stats median
#' @export
summaryGroup <- function(x, alpha = 0.05, ...){
    OK <- inherits(x, c("PMCMR"))
    if (!OK)
        stop ("Not an object of class PMCMR")
    if (!is.matrix(x$statistic))
        stop ("Matrix object$statistic not found.")
    pval <- as.numeric(x$p.value)
    stat <- as.numeric(x$statistic)
    grp1 <- as.numeric(c(col(x$p.value)))
    cnam <- colnames(x$p.value)
    grp2 <- as.numeric(c(row(x$p.value)))
    rnam <- rownames(x$p.value)
    STAT <- x$dist
    if (!is.null(x$alternative)) {
        if (x$alternative == "less"){
            H0 <- paste(rnam[grp2], "-", cnam[grp1], ">=", "0")
            PVAL <- paste("Pr(<", STAT, ")", sep="")
        } else if (x$alternative == "greater"){
            H0 <- paste(rnam[grp2], "-", cnam[grp1], "<=", "0")
            PVAL <- paste("Pr(>", STAT, ")", sep="")
        } else {
            H0 <- paste(rnam[grp2], "-", cnam[grp1], "==", "0")
            PVAL <- paste("Pr(>|", STAT, "|)", sep="")
        }
    } else {
        H0 <- paste(rnam[grp2], "-", cnam[grp1], "==", "0")
        PVAL <- paste("Pr(>|", STAT, "|)", sep="")
    }

    ##this makes only sense for simultaneous multiple comparisons
    ## it does not make sense for test with one control
    if(length(x$p.value[1,]) != length(x$p.value[,1])){
        stop("Homogeneous groups can not be determined by
many-to-one comparisons!")
    }
    if (!is.null(x$alternative)){
        if (x$alternative != "two.sided"){
            stop ("'summaryGrouped' is only possible for tests with alternative = 'two.sided'")
        }
    }
    ## make a grouped output
    ppval <- get.pvalues(x)
    ## only for today
    out.mcv <- multcompLetters(ppval, threshold = alpha)

    dat <- x$model
    if( any(c(grepl("Nemenyi", x$method),
              grepl("Dunn", x$method),
              grepl("Conover", x$method),
              grepl("Steel", x$method),
              grepl("Miller", x$method),
              grepl("Eisinga", x$method),
              grepl("Demsar", x$method),
              grepl("Siegel", x$method),
              grepl("Durbin", x$method),
              grepl("Anderson", x$method),
              grepl("BWS", x$method)))){
        xmedian <- tapply(dat$x, dat$g, median)
        xQ25 <- tapply(dat$x, dat$g, function(x) quantile(x, 0.25))
        xQ75 <- tapply(dat$x, dat$g, function(x) quantile(x, 0.75))
        xn <-  tapply(dat$x, dat$g, length)

        xdf <- data.frame(round(xmedian, 3),
                          round(xQ25, 3),
                          round(xQ75, 3),
                          xn,
                          out.mcv$Letters)
        rownames(xdf) <- c(colnames(x$statistic)[1],
                           rownames(x$statistic))
        colnames(xdf) <- c("median", "Q25", "Q75", "n", "Sig. group")

    } else {

        xmean <- tapply(dat$x, dat$g, mean)
        xn <- tapply(dat$x, dat$g, length)
        xsd <- tapply(dat$x, dat$g, sd)

        xdf <- data.frame(round(xmean, 3),
                          round(xsd, 3),
                          xn,
                          out.mcv$Letters)
        rownames(xdf) <- c(colnames(x$statistic)[1],
                           rownames(x$statistic))
        names(xdf) <- c("mean", "sd", "n", "Sig. group")
    }

    message("\n\tPairwise comparisons using ", x$method, "\n")
    message("data: ", x$data.name)
    if (!is.null(x$alternative)){
        message("alternative hypothesis: ", x$alternative)
    }
    message("P value adjustment method: ", x$p.adjust.method)
    message("Different letters indicate significant differences ",
        PVAL, " < ", alpha, "\n", sep="")
    print(xdf)
    invisible(x)
}

#' @title PMCMR Printing
#' @description
#' \code{print.PMCMR} is the \emph{PMCMR} method of the generic
#' \code{\link{print}} function which prints its argument
#' and returns it \emph{invisibly} (via \code{\link{invisible}(x)}).
#' @param x an object used to select a method.
#' @param \ldots further arguments. Currently ignored.
#' @aliases print.PMCMR
#' @method print PMCMR
#' @keywords print
#' @export
print.PMCMR <-
function(x, ...)
{
    message("\n\tPairwise comparisons using ", x$method, "\n")
    message("data: ", x$data.name, "\n")
    pp <- format.pval(x$p.value, 2, na.form="-")
    attributes(pp) <- attributes(x$p.value)
    print(pp, quote=FALSE, ...)
    message("\nP value adjustment method: ", x$p.adjust.method)
    if (!is.null(x$alternative)){
        message("alternative hypothesis: ", x$alternative)
    }
    invisible(x)
}


get.pvalues <-
    function(object, ...)
{
    OK <- inherits(object, c("PMCMR", "pairwise.htest", "trendPMCR"))
    if (!OK)
        stop ("Not an object of class ,", sQuote("PMCMR"), ", ", sQuote("pairwise.htest"), " or ", sQuote("trendPMCMR"))
    if (!is.matrix(object$p.value))
        stop ("Matrix object$p.value not found.")
    pval <- as.numeric(object$p.value)
    grp1 <- as.numeric(c(col(object$p.value)))
    cnam <- colnames(object$p.value)
    grp2 <- as.numeric(c(row(object$p.value)))
    rnam <- rownames(object$p.value)
    H0 <- paste(cnam[grp1],"-",rnam[grp2], sep="")
    OK <- !is.na(pval)
    out <- pval[OK]
    names(out) <- H0[OK]
    out
}
