## frdAllPairsNemenyiTest.R
## Part of the R package: PMCMRplus
##
## Copyright (C) 2014-2020 Thorsten Pohlert
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

#' @name frdAllPairsNemenyiTest
#' @title Nemenyi's All-Pairs Comparisons Test for Unreplicated Blocked Data
#' @description
#'  Performs Nemenyi's all-pairs comparisons tests of Friedman-type ranked data.
#' @details
#' For all-pairs comparisons in a two factorial unreplicated
#' complete block design
#' with non-normally distributed residuals, Nemenyi's test can be
#' performed on Friedman-type ranked data.
#'
#' A total of \eqn{m = k ( k -1 )/2} hypotheses can be tested.
#' The null hypothesis, H\eqn{_{ij}: \theta_i = \theta_j}, is tested
#' in the two-tailed case against the alternative,
#' A\eqn{_{ij}: \theta_i \ne \theta_j, ~~ i \ne j}.
#'
#' The \eqn{p}-values are computed from the studentized range distribution.
#'
#' @references
#' Demsar, J. (2006) Statistical comparisons of classifiers over multiple
#'  data sets, \emph{Journal of Machine Learning Research} \bold{7}, 1--30.
#'
#' Miller Jr., R. G. (1996) \emph{Simultaneous statistical inference}.
#'  New York: McGraw-Hill.
#'
#' Nemenyi, P. (1963), \emph{Distribution-free Multiple Comparisons}.
#'  Ph.D. thesis, Princeton University.
#'
#' Sachs, L. (1997) \emph{Angewandte Statistik}. Berlin: Springer.
#'
#' @keywords htest nonparametric
#' @concept friedmanranks
#'
#' @seealso
#' \code{\link{friedmanTest}}, \code{\link{friedman.test}},
#' \code{\link{frdAllPairsExactTest}}, \code{\link{frdAllPairsConoverTest}},
#' \code{\link{frdAllPairsMillerTest}}, \code{\link{frdAllPairsSiegelTest}}
#' @template class-PMCMR
#' @example examples/frdAllPairs.R
#' @export
frdAllPairsNemenyiTest <-
    function(y, ...)
        UseMethod("frdAllPairsNemenyiTest")

#' @rdname frdAllPairsNemenyiTest
#' @aliases frdAllPairsNemenyiTest.default
#' @method frdAllPairsNemenyiTest default
#' @template two-way-parms
#' @export
frdAllPairsNemenyiTest.default <-
    function(y, groups, blocks, ...)
    {
        ## 2019-10-16
        ## novel external function
        ans <- frdRanks(y, groups, blocks)
        r <- ans$r
        n <- nrow(r)
        k <- ncol(r)
        GRPNAMES <- colnames(r)

        ## continue
        R.mnsum <- colMeans(r)

        compare.stats <- function(i, j) {
            dif <- abs(R.mnsum[i] - R.mnsum[j])
            qval <- dif / sqrt(k * (k + 1) / (6 * n))
            return(qval)
        }


        PSTAT <- pairwise.table(compare.stats, GRPNAMES,
                                p.adjust.method = "none") * sqrt(2)

        PVAL <- 1 - ptukey(PSTAT, nmeans = k, df = Inf)
        METHOD <-
            c(
                "Nemenyi-Wilcoxon-Wilcox all-pairs test for a two-way",
                " balanced complete block design"
            )
        DIST <- "q"
        p.adjust.method <- "single-step"

        colnames(PSTAT) <- GRPNAMES[1:(k - 1)]
        rownames(PSTAT) <- GRPNAMES[2:k]
        colnames(PVAL) <- GRPNAMES[1:(k - 1)]
        rownames(PVAL) <- GRPNAMES[2:k]

        ans <- list(
            method = METHOD,
            data.name = ans$DNAME,
            p.value = PVAL,
            statistic = PSTAT,
            p.adjust.method = p.adjust.method,
            dist = DIST,
            model = ans$inDF
        )

        class(ans) <- "PMCMR"
        ans
    }

## taken from stats::friedman.test
#' @rdname frdAllPairsNemenyiTest
#' @method frdAllPairsNemenyiTest formula
#' @aliases frdAllPairsNemenyiTest.formula
#' @param formula a formula of the form \code{a ~ b | c} where
#'    \code{a, b} and \code{c} give the data values and
#'    the corresponding groups and blocks, respectively.
#' @param data an optional matrix or data frame (or similar: see
#'  \code{\link{model.frame}}) containing the variables in the
#'  formula \code{formula}.  By default the variables are taken from
#'  \code{environment(formula)}.
#' @param subset an optional vector specifying a
#'  subset of observations to be used.
#' @param na.action a function which indicates what should happen when
#'    the data contain \code{NA}s.  Defaults to \code{getOption("na.action")}.
#' @export
frdAllPairsNemenyiTest.formula <-
    function(formula, data, subset, na.action, ...)
    {
        if(missing(formula))
            stop("formula missing")
        ## <FIXME>
        ## Maybe put this into an internal rewriteTwoWayFormula() when
        ## adding support for strata()
        if((length(formula) != 3L)
           || (length(formula[[3L]]) != 3L)
           || (formula[[3L]][[1L]] != as.name("|"))
           || (length(formula[[3L]][[2L]]) != 1L)
           || (length(formula[[3L]][[3L]]) != 1L))
            stop("incorrect specification for 'formula'")
        formula[[3L]][[1L]] <- as.name("+")
        ## </FIXME>
        m <- match.call(expand.dots = FALSE)
        m$formula <- formula
        if(is.matrix(eval(m$data, parent.frame())))
            m$data <- as.data.frame(data)
        ## need stats:: for non-standard evaluation
        m[[1L]] <- quote(stats::model.frame)
        mf <- eval(m, parent.frame())
        DNAME <- paste(names(mf), collapse = " and ")
        names(mf) <- NULL
        y <- do.call("frdAllPairsNemenyiTest", as.list(mf))
        y$data.name <- DNAME
        y
    }
