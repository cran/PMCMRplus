## kwAllPairsNemenyiTest.R
## Part of the R package: PMCMR
##
## Copyright (C) 2014-2018 Thorsten Pohlert
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

#' @rdname kwAllPairsNemenyiTest
#' @title Nemenyi's All-Pairs Rank Comparison Test
#'
#' @description
#' Performs Nemenyi's non-parametric all-pairs comparison test
#' for Kruskal-type ranked data.
#'
#' @details
#' For all-pairs comparisons in an one-factorial layout
#' with non-normally distributed residuals Nemenyi's non-parametric test
#' can be performed. A total of \eqn{m = k(k-1)/2}
#' hypotheses can be tested. The null hypothesis
#' H\eqn{_{ij}: \theta_i(x) = \theta_j(x)} is tested in the two-tailed test
#' against the alternative
#' A\eqn{_{ij}: \theta_i(x) \ne \theta_j(x), ~~ i \ne j}.
#'
#' Let \eqn{R_{ij}} be the rank of \eqn{X_{ij}},
#' where \eqn{X_{ij}} is jointly ranked
#' from \eqn{\left\{1, 2, \ldots, N \right\}, ~~ N = \sum_{i=1}^k n_i},
#' then the test statistic under the absence of ties is calculated as
#'
#' \deqn{
#' t_{ij} = \frac{\bar{R}_j - \bar{R}_i}
#' {\sigma_R \left(1/n_i + 1/n_j\right)^{1/2}} \qquad \left(i \ne j\right),
#' }{%
#'  SEE PDF
#' }
#'
#' with \eqn{\bar{R}_j, \bar{R}_i} the mean rank of the
#' \eqn{i}-th and \eqn{j}-th group and the expected variance as
#' \deqn{
#'  \sigma_R^2 = N \left(N + 1\right) / 12.
#' }{%
#'  SEE PDF
#' }
#'
#' A pairwise difference is significant, if \eqn{|t_{ij}|/\sqrt{2} > q_{kv}},
#' with \eqn{k} the number of groups and \eqn{v = \infty}
#' the degree of freedom.
#'
#' Sachs(1997) has given a modified approach for
#' Nemenyi's test in the presence of ties for \eqn{N > 6, k > 4}
#' provided that the \code{\link{kruskalTest}} indicates significance:
#' In the presence of ties, the test statistic is
#' corrected according to \eqn{\hat{t}_{ij} = t_{ij} / C}, with
#' \deqn{
#'   C = 1 - \frac{\sum_{i=1}^r t_i^3 - t_i}{N^3 - N}.
#' }{%
#'  SEE PDF
#' }
#'
#' The function provides two different \code{dist}
#' for \eqn{p}-value estimation:
#' \describe{
#'  \item{Tukey}{The \eqn{p}-values are computed from the studentized
#'  range distribution (alias \code{\link[stats]{Tukey}}),
#'  \eqn{\mathrm{Pr} \left\{ t_{ij} \sqrt{2} \ge q_{k\infty\alpha} | mathrm{H} \right\} = \alpha}.}
#'  \item{Chisquare}{The \eqn{p}-values are computed from the
#'  \code{\link[stats]{Chisquare}} distribution with \eqn{v = k - 1} degree
#'  of freedom.}
#' }
#'
#' @references
#' Nemenyi, P. (1963) \emph{Distribution-free Multiple Comparisons}.
#'  Ph.D. thesis, Princeton University.
#'
#' Sachs, L. (1997) \emph{Angewandte Statistik}. Berlin: Springer.
#'
#' Wilcoxon, F., Wilcox, R. A. (1964)
#' \emph{Some rapid approximate statistical procedures}.
#' Pearl River: Lederle Laboratories.
#'
#' @template class-PMCMR
#' @keywords nonparametric
#' @concept kruskalranks
#'
#' @seealso
#' \code{\link[stats]{Tukey}}, \code{\link[stats]{Chisquare}},
#' \code{\link[stats]{p.adjust}}, \code{\link{kruskalTest}},
#' \code{\link{kwAllPairsDunnTest}}, \code{\link{kwAllPairsConoverTest}}
#' @example examples/kwAllPairsMC.R
#' @export
kwAllPairsNemenyiTest <- function(x, ...) UseMethod("kwAllPairsNemenyiTest")

#' @rdname kwAllPairsNemenyiTest
#' @method kwAllPairsNemenyiTest default
#' @aliases kwAllPairsNemenyiTest.default
#' @template one-way-parms
#' @param dist the distribution for determining the p-value.
#'  Defaults to \code{"Tukey"}.
#' @importFrom stats ptukey
#' @importFrom stats pchisq
#' @importFrom stats complete.cases
#' @importFrom stats pairwise.table
#' @export
kwAllPairsNemenyiTest.default <-
function(x, g, dist = c("Tukey","Chisquare"), ...){
        ## taken from stats::kruskal.test
        if (is.list(x)) {
            if (length(x) < 2L)
                stop("'x' must be a list with at least 2 elements")
            DNAME <- deparse(substitute(x))
            x <- lapply(x, function(u) u <- u[complete.cases(u)])
            k <- length(x)
            l <- sapply(x, "length")
            if (any(l == 0))
                stop("all groups must contain data")
            g <- factor(rep(1 : k, l))
            # Get distribution from formula call, defaults to Tukey
            dist <- ifelse(!is.null(x$dist), x$dist, "Tukey")
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
        dist <- match.arg(dist)
        x.rank <- rank(x)
        R.bar <- tapply(x.rank, g, mean,na.rm=T)
        R.n <- tapply(!is.na(x), g, length)
        g.unique <- unique(g)
        k <- length(g.unique)
        n <- sum(R.n)

	if(dist == "Chisquare") {
          METHOD <- "Nemenyi's all-pairs test with chi-square approximation"
          DIST <- "X^2"
         compare.stats <- function(i,j) {
            dif <- abs(R.bar[i] - R.bar[j])
            A <- n * (n+1) / 12
            B <- (1 / R.n[i] + 1 / R.n[j])
            chisqval <- dif^2 / (A * B)
            return(chisqval)
         }

        PSTAT <- pairwise.table(compare.stats,levels(g),
                                p.adjust.method="none" )
          C <- gettiesKruskal(x.rank)
          if (C != 1) {
              warning("Ties are present. Chi-sq was corrected for ties.")
          }
          ## Must be devided by C, same as in stats::kruskal.test
          PSTAT <- PSTAT / C
          PVAL <- 1 - pchisq(PSTAT, df=(k-1))
        } else {
            METHOD <- "Tukey-Kramer-Nemenyi all-pairs test with Tukey-Dist approximation"
            DIST <- "q"
            compare.stats <- function(i,j) {
                dif <- abs(R.bar[i] - R.bar[j])
                qval <- dif / sqrt((n * (n + 1) / 12) *
                                   (1/R.n[i] + 1/R.n[j] ))
                return(qval)
            }
            PSTAT <- pairwise.table(compare.stats,levels(g),
                                    p.adjust.method="none" )*sqrt(2)
            C <- gettiesKruskal(x.rank)
            if (C != 1) warning("Ties are present, p-values are not corrected.")
            PVAL <- 1 - ptukey(PSTAT, nmeans=k, df=Inf)
        }
    MODEL <- data.frame(x, g)
    ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL,
                statistic = PSTAT, p.adjust.method = "single-step",
                model = MODEL, dist = DIST, alternative = "two.sided")
    class(ans) <- "PMCMR"
    ans
}


#' @rdname kwAllPairsNemenyiTest
#' @method kwAllPairsNemenyiTest formula
#' @aliases kwAllPairsNemenyiTest.formula
#' @template one-way-formula
#' @export
kwAllPairsNemenyiTest.formula <-
function(formula, data, subset, na.action, dist = c("Tukey","Chisquare"), ...)
{
    mf <- match.call(expand.dots=FALSE)
    m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    dist <- match.arg(dist)
    mf[[1L]] <- quote(stats::model.frame)

   if(missing(formula) || (length(formula) != 3L))
        stop("'formula' missing or incorrect")
    mf <- eval(mf, parent.frame())
    if(length(mf) > 2L)
       stop("'formula' should be of the form response ~ group")
    DNAME <- paste(names(mf), collapse = " by ")
    names(mf) <- NULL
    dist <- match.arg(dist)
    y <- do.call("kwAllPairsNemenyiTest", c(as.list(mf), dist = dist))
    y$data.name <- DNAME
    y
}
