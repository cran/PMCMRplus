## manyOneUTest.R
## Part of the R package: PMCMR
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

#' @name manyOneUTest
#' @title Multiple Comparisons with One Control (U-test)
#' @description
#' Performs pairwise comparisons of multiple group levels with
#' one control.
#' @details
#' This functions performs Wilcoxon, Mann and Whitney's U-test
#' for a one factorial design where each factor level is tested against
#' one control (\eqn{m = k -1} tests). As the data are re-ranked
#' for each comparison, this test is only suitable for
#' balanced (or almost balanced) experimental designs.
#'
#' For the two-tailed test and \code{p.adjust.method = "single-step"}
#' the multivariate normal distribution is used for controlling
#' Type 1 error and to calculate p-values. Otherwise,
#' the p-values are calculated from the standard normal distribution
#' with any latter p-adjustment as available by \code{\link{p.adjust}}.
#'
#' @template class-PMCMR
#'
#' @references
#' OECD (ed. 2006) \emph{Current approaches in the statistical analysis
#' of ecotoxicity data: A guidance to application}, OECD Series
#' on testing and assessment, No. 54.
#' @seealso
#' \code{\link{wilcox.test}}, \code{\link{pmvnorm}}, \code{\link{Normal}}
#' @concept wilcoxonranks
#' @keywords htest nonparametric
#' @export
manyOneUTest <- function(x, ...) UseMethod("manyOneUTest")

#' @rdname manyOneUTest
#' @method manyOneUTest default
#' @aliases manyOneUTest.default
#' @template one-way-parms
#' @param alternative the alternative hypothesis. Defaults to \code{two.sided}.
#' @param p.adjust.method method for adjusting p values
#' (see \code{\link{p.adjust}})
#' @importFrom stats pnorm
#' @importFrom stats p.adjust
#' @importFrom stats p.adjust.methods
#' @importFrom mvtnorm pmvnorm
#' @importFrom stats complete.cases
#' @export
manyOneUTest.default <-
    function(x, g, alternative = c("two.sided", "greater", "less"),
             p.adjust.method = c("single-step", p.adjust.methods),...){
        ## taken from stats::kruskal.test

        if (is.list(x)) {
            stop("'x' must be a list with at least 2 elements")
            DNAME <- deparse(substitute(x))
            x <- lapply(x, function(u) u <- u[complete.cases(u)])
            k <- length(x)
            l <- sapply(x, "length")
            if (any(l == 0))
                stop("all groups must contain data")
            g <- factor(rep(1 : k, l))
            if (is.null(x$alternative)){
                alternative <- "two.sided"
            } else {
                alternative <- x$alternative
            }
            if(is.null(x$p.adjust.method)){
                p.adjust.method <- "single-step"
            } else {
                p.adjust.method <- x$p.adjust.method
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

        # check arguments
	alternative <- match.arg(alternative)
	p.adjust.method <- match.arg(p.adjust.method)


        k <- nlevels(g)
        n <- tapply(x, g, length)
        glev <- levels(g)

        # Function to get ties for tie adjustment
        getties <- function(x){
            t <- table(x)
            C <- sum((t^3 - t) / 12)
            C
        }

    # function for pairwise comparisons with one control
        compare.stats <-function(i){
            n1 <- n[1]
            n2 <- n[i]
            xraw <- c(x[g==glev[1]], x[g==glev[i]])
            rankx <- rank(xraw)
            lev <- c(g[g==glev[1]], g[g==glev[i]])
            R <- tapply(rankx, lev, sum)
            U <- c(n1* n2 + (n1 * (n1 + 1) / 2),
                   n1 * n2 + (n2 * (n2 + 1) / 2)) - R
            Umn <- min(U)
            S <- n1 + n2
            VAR <- (n1 * n2 / (S * (S - 1))) * ((S^3 - S) / 12 - getties(rankx))
            PSTAT <- (Umn - n1 * n2 / 2) / sqrt(VAR)
            if (R[1] < R[2]){
                PSTAT <- PSTAT * (-1)
            }
            PSTAT
        }

        # compute values
	pstat <- sapply(2:k, function(i) compare.stats(i))

	if (p.adjust.method != "single-step"){
	    # use normal approximation, possibly with p.adjust other than
            # single-step
            if (alternative == "two.sided"){
                pval <- 2 * pnorm(abs(pstat), lower.tail = FALSE)
            }
            else if(alternative == "greater"){
		pval <- pnorm(pstat, lower.tail = FALSE)
            } else {
		pval <- pnorm(pstat)
            }
            padj <- p.adjust(pval, method = p.adjust.method)
	} else {

	 ## use function pmvnorm of package mvtnorm
            m <- k - 1
            # correlation matrix
            if(m < 2){
                stop("for 'p.adjust.method = single-step' at least\n
                      k = 3 groups are required.")
                cr <- NULL
            } else {
                ni <- tapply(x, g, length)
                n0 <- ni[1]
                nn <- ni[2:k]
                cr <- diag(m)

                # create cr for unequal (or equal) sample sizes
                for (i in 1:(m-1)){
                    for (j in (i+1):m){
				 #cr[i,j] <- 0.5
				 #cr[j,i] <- 0.5
                        cr[i,j] <- ((nn[i] * nn[j]) /
                                    ((nn[i] + n0) * (nn[j] + n0)))^(1/2)
                        cr[j,i] <- cr[i, j]
                        cr[j,i] <- cr[i, j]

                    }
                }
            }
            if (alternative == "two.sided"){
                padj <- sapply(pstat, function(x)
                    1 - pmvnorm(lower = -rep(abs(x), m),
                                upper = rep(abs(x), m),
                                corr = cr))
            } else if (alternative == "greater"){
                padj <- sapply(pstat, function(x)
                    1 - pmvnorm(lower = -Inf,
                                upper = rep(x, m),
                                corr = cr))
            } else {
                padj <- sapply(pstat, function(x)
                    1 - pmvnorm(lower = rep(x, m),
                                upper = Inf,
                                corr = cr))
            }
        }

	GRPNAMES <- levels(g)
        PSTAT <- cbind(pstat)
	PVAL <- cbind(padj)
	colnames(PSTAT) <- GRPNAMES[1]
	colnames(PVAL) <- GRPNAMES[1]
	rownames(PSTAT) <- GRPNAMES[2:k]
	rownames(PVAL) <- GRPNAMES[2:k]

	DIST <- "z"

        METHOD <- paste("Wilcoxon, Mann, Whittney U-test\n",
		           "\tfor multiple comparisons with one control",
						sep="")
        MODEL <- data.frame(x, g)
        ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL,
                    statistic = PSTAT, p.adjust.method = p.adjust.method,
                    model = MODEL, dist=DIST, alternative = alternative)
        class(ans) <- "PMCMR"
        ans
}

#' @rdname manyOneUTest
#' @method manyOneUTest formula
#' @aliases manyOneUTest.formula
#' @template one-way-formula
#' @export
manyOneUTest.formula <-
    function(formula, data, subset, na.action,
             alternative = c("two.sided", "greater", "less"),
             p.adjust.method = c("single-step", p.adjust.methods),...)
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
    p.adjust.method <- match.arg(p.adjust.method)
    names(mf) <- NULL
    y <- do.call("manyOneUTest", c(as.list(mf), alternative = alternative,
                                  p.adjust.method = p.adjust.method))
    y$data.name <- DNAME
    y
}
