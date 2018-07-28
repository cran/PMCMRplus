# steelTest.R
# Part of the R package: PMCMR
#
# Copyright (C) 2018 Thorsten Pohlert
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

#' @name steelTest
#' @title Steel's Many-to-One Rank Test
#' @description
#' Performs Steel's non-parametric many-to-one comparison
#' test for Wilcox-type ranked data.
#' @details
#' For many-to-one comparisons (pairwise comparisons with one control)
#' in an one-factorial balanced layout with non-normally distributed
#' residuals Steels's non-parametric single-step test can be performed.
#' Let there be \eqn{k} treatment levels (excluding the control),
#' then \eqn{k} pairwise comparisons can be performed between
#' the \eqn{i}-th treatment level and the control.
#' H\eqn{_i: \theta_0 = \theta_i} is tested in the one-tailed case (less) against
#' A\eqn{_i: \theta_0 > \theta_i, ~~ (1 \le i \le k)}.
#'
#' For each control - treatment level the data are ranked in increasing order.
#' The ranksum \eqn{R_i} for the \eqn{i}-th treatment level is compared
#' to a critical \eqn{R} value and is significantly(\eqn{p = 0.05}) less,
#' if \eqn{R_i \le R}. For the \code{alternative = "greater"} the sign is changed.
#'
#' The function does not return p-values. Instead the critical \eqn{R}-values
#' as given in the tables of USEPA (2002) for \eqn{\alpha = 0.05} (one-sided, less)
#' are looked up according to the balanced sample sizes (\eqn{n}) and the order number of the
#' dose level (\eqn{i}).
#'
#' @note
#' Steel's Many-to-One Rank test is only applicable for balanced designs and
#' directional hypotheses. An error message will occur, if the design is unbalanced.
#' In the current implementation, only one-sided tests on
#' the level of \eqn{\alpha = 0.05} can be performed.
#'
#' @return
#' A list with class \code{"steelTest"} containing the following components:
#' \describe{
#'  \item{method}{a character string indicating what type of test was performed.}
#'  \item{data.name}{a character string giving the name(s) of the data.}
#'  \item{statistic}{lower-triangle matrix of the ranksum for the i-th tretent level}
#'  \item{R.crit}{lower-triangle matrix of critical R-values for \eqn{\alpha = 0.05}.}
#'  \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{model}{a data frame of the input data.}
#' \item{dist}{a string that denotes the test distribution.}
#' }
#' There are print and summary methods available.
#'
#' @section Source:
#' The critical rank sum values were taken from Table E.5 of USEPA (2002).
#'
#'  USEPA (2002) \emph{Short-term Methods for Estimating the
#'  Chronic Toxicity of Effluents and Receiving
#'  Waters to Freshwater Organisms}, 4th edition, EPA-821-R-02-013.
#'
#' @seealso
#' \code{\link[stats]{wilcox.test}}, \code{\link[stats]{pairwise.wilcox.test}},
#' \code{\link{manyOneUTest}},
#' \code{\link{shirleyWilliamsTest}}, \code{\link{kwManyOneDunnTest}},
#' \code{\link{kwManyOneNdwTest}}, \code{\link{kwManyOneConoverTest}},
#' \code{\link{print.steel}}, \code{\link{summary.steel}}
#'
#' @references
#' Steel, R. G. D. (1959) A multiple comparison rank sum test:
#' treatments versus control, \emph{Biometrics} \bold{15}, 560--572.
#'
#' @keywords htest
#' @examples
#' ## Example from Sachs (1997, p. 402)
#' x <- c(106, 114, 116, 127, 145,
#' 110, 125, 143, 148, 151,
#' 136, 139, 149, 160, 174)
#' g <- gl(3,5)
#' levels(g) <- c("0", "I", "II")
#'
#' ## Steel's Test
#' steelTest(x ~ g)
#'
#'
#' ## Example from USEPA (2002):
#' ## Reproduction data from a Ceriodaphnia dubia
#' ## 7-day chronic test to several concentrations
#' ## of effluent. Dose level 50% is excluded.
#' x <- c(20, 26, 26, 23, 24, 27, 26, 23, 27, 24,
#' 13, 15, 14, 13, 23, 26, 0, 25, 26, 27,
#' 18, 22, 13, 13, 23, 22, 20, 22, 23, 22,
#' 14, 22, 20, 23, 20, 23, 25, 24, 25, 21,
#' 9, 0, 9, 7, 6, 10, 12, 14, 9, 13,
#' rep(0,10))
#' g <- gl(6, 10)
#' levels(g) <- c("Control", "3%", "6%", "12%", "25%", "50%")
#'
#' ## NOEC at 3%, LOEC at 6%
#' steelTest(x ~ g, subset = g != "50%", alternative = "less")
#'
#'
#'
#' @export
steelTest <- function(x, ...)
  UseMethod("steelTest")

#' @rdname steelTest
#' @aliases steelTest.default
#' @method steelTest default
#' @template one-way-parms
#' @param alternative the alternative hypothesis. Defaults to \code{greater}
#' @export
steelTest.default <-
  function(x,
           g,
           alternative = c("greater", "less"),
           ...)
  {
    ## taken from stats::steelTest
    if (is.list(x)) {
      if (length(x) < 2L)
        stop("'x' must be a list with at least 2 elements")
      DNAME <- deparse(substitute(x))
      x <- lapply(x, function(u)
        u <- u[complete.cases(u)])
      k <- length(x)
      l <- sapply(x, "length")
      if (any(l == 0))
        stop("all groups must contain data")
      g <- factor(rep(1:k, l))
      alternative <- x$alternative
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
      if (k < 2) {
        stop("all observations are in the same group")
      }
    }

    alternative <- match.arg(alternative)

    ## check for balanced design
    ni <- tapply(x, g, length)
    if (!length(unique((ni))) == 1) {
      stop("Steel's test is only valid for balanced designs.")
    }

    ## Check for number of n
    n <- unique(ni)
    if (n < 4) {
      stop("Steel's test requires a minimum of n = 4 replicates")
    } else if (n > 20) {
      stop("Steel's critical R-values are only available for up to n = 20 replicates")
    }

    ##check for number of treatments k excluding control
    kk <- nlevels(g)
    k <- kk - 1
    if (k > 9){
      stop("Steel's critical R-values are only available for up to k = 9 dose levels")
    } else if (n == 4 & k > 6) {
      stop("Steel's critical R-values for n = 4 are only available for up to k = 6 dose levels")
    }


    xold <- x
    ## Critical values are for less
    if (alternative == "greater") {
      x <- -x
    }

    ## Simply rank sums of treatment levels
    Dose <- levels(g)[2:kk]
    Ctrl <- levels(g)[1]

    R <- sapply(Dose, function(j){
      tmpR <- rank(x[g == Ctrl | g == j])
      sum(tmpR[(n+1):(2*n)])
    } )

    ## Get critical value
    K <- paste(k)
    N <- paste(n)

    Rcrit <- rep(TabCrit$steel.R005[N, K], k)


    ## Create output matrices
    STAT <- cbind(Ctrl = R)
    row.names(STAT) <- Dose
    STATCRIT <- cbind(Ctrl = Rcrit)
    row.names(STATCRIT) <- Dose

    METHOD <- paste("Steel's Many-to-One-Rank Test")
    ans <- list(
      method = METHOD,
      data.name = DNAME,
      R.crit = STATCRIT,
      statistic = STAT,
      alternative = alternative,
      dist = "R",
      model = data.frame(x=xold, g = g)
    )
    class(ans) <- "steel"
    ans
  }

#' @rdname steelTest
#' @aliases steelTest.formula
#' @method steelTest formula
#' @template one-way-formula
#' @export
steelTest.formula <-
  function(formula, data, subset, na.action, alternative = c("greater", "less"), ...)
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
    names(mf) <- NULL
    y <- do.call("steelTest", c(as.list(mf), alternative = alternative))
    y$data.name <- DNAME
    y
  }
