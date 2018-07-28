## frdManyOneNemenyiTest.R
## Part of the R package: PMCMR
##
## Copyright (C)  2017, 2018 Thorsten Pohlert
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
#' @rdname frdManyOneNemenyiTest
#' @title  Nemenyi's Many-to-One Test
#' for Unreplicated Blocked Data
#'
#' @description
#' Performs Nemenyi's non-parametric many-to-one comparison test
#' for Friedman-type ranked data.
#'
#' @details
#' For many-to-one comparisons (pairwise comparisons with one control)
#' in a two factorial unreplicated complete block design
#' with non-normally distributed residuals, Nemenyi's test can be
#' performed on Friedman-type ranked data.
#'
#' Let there be \eqn{k} groups including the control,
#' then the number of treatment levels is \eqn{m = k - 1}.
#' A total of \eqn{m} pairwise comparisons can be performed between
#' the \eqn{i}-th treatment level and the control.
#' H\eqn{_i: \theta_0 = \theta_i} is tested in the two-tailed case against
#' A\eqn{_i: \theta_0 \ne \theta_i, ~~ (1 \le i \le m)}.
#'
#' The \eqn{p}-values are computed from the multivariate normal distribution.
#'
#' @references
#' Hollander, M., Wolfe, D. A., Chicken, E. (2014),
#' \emph{Nonparametric Statistical Methods}. 3rd ed. New York: Wiley. 2014.
#'
#' Miller Jr., R. G. (1996), \emph{Simultaneous Statistical Inference}.
#'  New York: McGraw-Hill.
#'
#' Nemenyi, P. (1963), \emph{Distribution-free Multiple Comparisons}.
#'  Ph.D. thesis, Princeton University.
#'
#' Siegel, S., Castellan Jr., N. J. (1988), \emph{Nonparametric
#'  Statistics for the Behavioral Sciences}. 2nd ed.
#'  New York: McGraw-Hill.
#'
#' Zarr, J. H. (1999), \emph{Biostatistical Analysis}. 4th ed.
#' Upper Saddle River: Prentice-Hall.
#'
#' @concept FriedmanTest
#' @concept Rank
#' @concept ManyToOne
#' @keywords htest nonparametric
#' @example examples/frdManyOne.R
#' @seealso
#' \code{\link{friedmanTest}}, \code{\link[stats]{friedman.test}},
#' \code{\link{frdManyOneExactTest}}, \code{\link{frdManyOneDemsarTest}}
#' \code{\link[mvtnorm]{pmvnorm}}
#'
#' @template class-PMCMR
#' @export
frdManyOneNemenyiTest <- function(y, ...) UseMethod("frdManyOneNemenyiTest")

#' @rdname frdManyOneNemenyiTest
#' @method frdManyOneNemenyiTest default
#' @aliases frdManyOneNemenyiTest.default
#' @param alternative the alternative hypothesis. Defaults to \code{two.sided}.
#' @template two-way-parms
#' @importFrom mvtnorm pmvnorm
#' @export
frdManyOneNemenyiTest.default <-
    function(y, groups, blocks, alternative =
                                    c("two.sided", "greater", "less"), ...)
{
    if ((is.matrix(y)) | (is.data.frame(y))) {
        ## corrected 4. Jun 2017
        DNAME <- paste(deparse(substitute(y)))
        GRPNAMES <- colnames(y)
        k <- length(GRPNAMES)
        BLOCKNAMES <- rownames(y)
        n <- length(BLOCKNAMES)
        groups <- factor(rep(GRPNAMES, times = n))
        blocks <- factor(rep(BLOCKNAMES, each = k))
        y <- as.vector(t(y))
    }
    else {
        if (any(is.na(groups)) || any(is.na(blocks)))
            stop("NA's are not allowed in groups or blocks")
        if (any(diff(c(length(y), length(groups), length(blocks)))))
            stop("y, groups and blocks must have the same length")
        if (any(table(groups, blocks) != 1))
                stop("Not an unreplicated complete block design")

        DNAME <- paste(deparse(substitute(y)), ",",
                       deparse(substitute(groups)), "and",
                       deparse(substitute(blocks)) )
        groups <- factor(groups)
        blocks <- factor(blocks)
        k <- nlevels(groups)
        n <- nlevels(blocks)
        GRPNAMES <- as.character(groups[1:k])
    }

    ## Check arguments
    alternative <- match.arg(alternative)

    mat <- matrix(y, nrow = n, ncol = k, byrow = TRUE)
    r <- t(apply(mat, 1L, rank))


    METHOD <- c("Nemenyi-Wilcoxon-Wilcox-Miller many-to-one test",
                 " for a two-way balanced complete block design")
    df <- k-1

    # correlation matrix for balanced design
    cr <- matrix(0.5 , ncol=df, nrow=df)
    diag(cr) <- 1.0

    Ri <- colSums(r)
    compNem <- function(j){
        dif <- Ri[j] - Ri[1]
        val <- dif / sqrt(n * k * (k + 1) / 6)
        return(val)
    }
    PSTATv <- rep(NA, df)
    for (j in 2:k) {PSTATv[j-1] <- compNem(j)}

    if (alternative == "two.sided"){
        PADJv <- sapply(PSTATv, function(x)
            1 - pmvnorm(lower = -rep(abs(x), df),
                        upper = rep(abs(x), df),
                        corr = cr))
    } else if (alternative == "greater"){
        PADJv <- sapply(PSTATv, function(x)
            1 - pmvnorm(lower = -Inf,
                        upper = rep(x, df),
                        corr = cr))
    } else {
        PADJv <- sapply(PSTATv, function(x)
            1 - pmvnorm(lower = rep(x, df),
                        upper = Inf,
                        corr = cr))
    }
    LNAME <- GRPNAMES[2:k]
    DIST <- "z"
    p.adjust.method = "single-step"
    ## build matrix
    PSTAT <- matrix(data=PSTATv, nrow = (k-1), ncol = 1,
                    dimnames = list(LNAME, GRPNAMES[1]))
    PVAL <- matrix(data=PADJv, nrow = (k-1), ncol = 1,
                   dimnames = list(LNAME, GRPNAMES[1]))
    MODEL <- data.frame(x = y, g = groups, b = blocks)
    ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL,
                statistic = PSTAT, p.adjust.method = p.adjust.method,
                alternative = alternative, dist = DIST, model = MODEL)
    class(ans) <- "PMCMR"
    ans
}
