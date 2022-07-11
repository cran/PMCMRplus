## power.tukey.test.R
## Part of the R package: PMCMR
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

#' @title Power Calculations for Balanced Tukey's
#' Multiple Comparison Test
#'
#' @description
#' Compute average per-pair power of Tukey's test for
#' multiple comparison of means.
#'
#' @param n number of observations (per group)
#' @param groups number of groups
#' @param delta true difference in means
#' @param within.var within group variance
#' @param sig.level significance level (Type I error probability)
#'
#' @details
#' The function has implemented the following Eq.
#' to estimate average per-pair power for two-sided tests:
#'
#' \deqn{
#'  1 - \beta = 1 - t(q_{\alpha v k}/\sqrt{2}, v, \mathrm{ncp}) +
#'   t(-q_{\alpha v k}/\sqrt{2}, v, \mathrm{ncp}),
#' }
#'
#' with \eqn{q_{\alpha v k}} the upper \eqn{\alpha} quantile of
#' the studentised range distribution, with \eqn{v = k (n - 1)}
#' degree of freedom and \eqn{k} the number of groups;
#' and \eqn{t(. ~\mathrm{ncp})}
#' the probability function of the non-central student t-distribution
#' with non-centrality parameter
#'
#' \deqn{
#'  \mathrm{ncp} = |\Delta| / \sqrt{s_{\mathrm{in}}^2 ~ 2 / n }.
#' }
#'
#' @seealso
#' \code{\link[stats]{TDist}} \code{\link[stats]{Tukey}}
#' \code{\link{powerMCTests}}
#'
#' @source
#' The Eqs. were taken from Lecture 5, *Determining Sample Size*,
#' Statistics 514, Fall 2015, Purdue University, IN, USA.
#'
#' @returns
#' Object of class \sQuote{\code{power.htest}},
#' a list of the arguments
#' (including the computed one) augmented with
#' \code{method} and \code{note} elements.
#'
#' @keywords htest
#' @examples
#' power.tukey.test(n = 11, groups = 5, delta = 30,
#'  within.var = 333.7)
#'
#' ## compare with t-test, Bonferroni-correction
#' power.t.test(n = 11, delta = 30, sd = sqrt(333.7),
#' sig.level = 0.05 / 10)
#'
#' \dontrun{
#' powerMCTests(mu = c(rep(0,4), 30), n = 11,
#'  parms = list(mean = 0,sd = sqrt(333.7)),
#'  test = "tukeyTest")
#' }
#' @export
power.tukey.test <- function(n,
                             groups,
                             delta,
                             within.var,
                             sig.level = 0.05) {
  df <- groups * (n - 1)
  ncp <- abs(delta) / (sqrt(within.var * 2 / n))

  ## Tukey distribution
  #' @importFrom stats qtukey
  stat <- qtukey(
    p = sig.level,
    nmeans = groups,
    df = df,
    lower.tail = FALSE
  ) / sqrt(2)

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
      method = "Tukey-Test power calculation",
      note = "n is number in *each* group",
      power = power,
      alternative = "two-sided"
    )

  class(out) <- "power.htest"
  out

}
