## power.dunnett.test.R
## Part of the R package: PMCMRplus
##
## Copyright (C) 2022 Thorsten Pohlert
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
#' @title Power Calculations for Balanced Dunnett's
#' Many-to-One Comparison Test
#'
#' @description
#' Compute average per-pair power of Dunnetts's multiple comparison
#' test with one control.
#'
#' @param n Number of observations (per group)
#' @param groups Number of groups (including control)
#' @param delta true difference in means
#' @param within.var Within group variance
#' @param sig.level Significance level (Type I error probability)
#'
#' @details
#' The function has implemented the following Eq.
#' to estimate average per-pair power for two-sided tests:
#'
#' \deqn{
#'  1 - \beta = 1 - t( T_{\alpha \rho v}, v, \mathrm{ncp}) +
#'   t(-T_{\alpha \rho v}, v, \mathrm{ncp}),
#' }
#'
#' with \eqn{T_{\alpha \rho v}} the two-sided
#' \eqn{\alpha} quantile of
#' the multivariate t-distribution, with \eqn{v = k (n - 1)}
#' degree of freedom, \eqn{k} the number of groups
#' and correlation matrix \eqn{\rho_{ij} = 0.5 ~ (i \neq j)}.
#'
#' The non-centrality parameter for the non-central student t-distribution
#' is
#'
#' \deqn{
#'  \mathrm{ncp} = |\Delta| / \sqrt{s_{\mathrm{in}}^2 ~ 2 / n }.
#' }
#'
#' @inherit power.tukey.test source
#'
#' @inherit power.tukey.test return
#'
#' @note
#' The results for power are seed depending.
#'
#' @seealso
#' \code{\link[stats]{TDist}} \code{\link[mvtnorm]{qmvt}}
#' \code{\link{powerMCTests}}
#'
#' @keywords htest
#' @examples
#' set.seed(113)
#' power.dunnett.test(n = 9, groups = 5, delta = 30,
#'  within.var = 333.7)
#'
#' ## compare with t-test, bonferroni corrected
#' power.t.test(n = 9, delta = 30, sd = sqrt(333.7),
#' sig.level = 0.05 / 4)
#'
#' \dontrun{
#' ## asymptotic Monte-Carlo power analysis
#'  set.seed(113)
#'  powerMCTests(mu = c(rep(0,4), 30), n = 9,
#'  parms = list(mean = 0, sd = sqrt(333.7)),
#'  test = "dunnettTest", alternative = "two.sided")
#' }
#' @export
power.dunnett.test <- function(n,
                             groups,
                             delta,
                             within.var,
                             sig.level = 0.05) {
  ## balanced design
  df <- groups * (n - 1)
  ncp <- abs(delta) / (sqrt(within.var * 2 / n))

  ## m tests
  m <- groups - 1

  ## new call to qDunnett
  stat <- qDunnett(p = sig.level/ 2,
                   n0 = n,
                   n = rep(n, m))

  ## correlation matrix, balanced design
  # corM <- matrix(0.5, m, m)
  # diag(corM) <- 1
  #
  # #' @importFrom mvtnorm qmvt
  # stat <- qmvt(p = 1-sig.level,
  #              tail = "both.tails", ## two-sided
  #              df = df,
  #              corr = corM)$quantile

  #' @importFrom stats pt
  power <-
    pt(q = -stat,
       df = df,
       ncp = ncp) +
    1 - pt(q = stat,
           df = df,
           ncp = ncp)

  power <- min(1, power)

  out <-
    list(
      n = n,
      groups = groups,
      delta = delta,
      within.var = within.var,
      method = "Dunnett-Test power calculation",
      note = "n is number in *each* group",
      power = power,
      alternative = "two-sided"
    )

  class(out) <- "power.htest"
  out
}
