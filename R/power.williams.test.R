## power.williams.test.R
## Part of the R package: PMCMRplut
##
## Copyright (C) 2022-2023 Thorsten Pohlert
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
#' @title Power calculations for
#' minimum detectable difference of the Williams' test
#' @description
#' Compute the power of a Williams' test,
#' or determine parameters to obtain a target power.
#'
#' @param n number of observations (per group).
#' @param k number of treatment groups.
#' @param delta clinically meaningful minimal difference
#' (between a treatment group and control).
#' @param sd common standard deviation.
#' @param power power of test (1 minus Type II error probability).
#' @param \ldots further arguments, currently ignored.
#'
#' @details
#' Exactly one of the parameters \code{n} or \code{power}
#' must be passed as \code{NULL}, and that
#' parameter is determined from the others.
#'
#' The function has implemented the following Eq. in order to
#' estimate power (Chow et  al. 2008):
#'
#' \deqn{
#'  1 - \beta = 1 - \Phi \left(T_{K \alpha v} -
#'  |\Delta| / \sigma \sqrt{2/n}\right)
#' }{%
#'  1 - beta = 1 - Phi( TKalpha - [|Delta| / sigma * sqrt(2/n)])
#' }
#'
#' with \eqn{|\Delta|}{|Delta|} the clinically meaningful minimal difference,
#' \eqn{T_{K \alpha v}}{TKalpha} the critical Williams' t-statistic
#' for \eqn{\alpha = 0.05}{alpha = 0.05}, \eqn{v = \infty} degree of freedom
#' and \eqn{\Phi}{Phi} the probability function of the standard normal function.
#'
#' The required sample size (balanced design) is estimated
#' based on the expression as given by the PASS manual, p. 595-2:
#'
#' \deqn{
#' n = 2 \sigma^2 ~ \left(T_{K \alpha v} + z_{\beta} \right)^2 ~ / ~ \Delta^2
#' }{%
#' n = 2 * sigma^2 * (TKaplha + zBeta)^2 / Delta^2
#' }
#'
#' @note
#' The current function calculates power for \code{sig.level = 0.05}
#' significance level (Type I error probability) only (one-sided test).
#'
#' @seealso
#' \code{\link[stats]{optimise}} \code{\link{williamsTest}}
#'
#' @return
#' Object of class \sQuote{\code{power.htest}}, a list of the arguments
#' (including the
#' computed one) augmented with method and note elements.
#'
#' @references
#' Chow, S.-C., Shao, J., Wan, H., 2008,
#' *Sample Size Calculations in Clinical Research*, 2nd ed,
#' Chapman & Hall/CRC: Boca Raton, FL.
#'
#' @examples
#' ## Chow et al. 2008, p. 288 depicts 53 (rounded),
#' ## better use ceiling for rounding
#' power.williams.test(power = 0.8, k = 3, delta = 11, sd = 22)
#' power.williams.test(n = 54, k = 3, delta = 11, sd = 22)
#'
#' ## PASS manual example:
#' ## up-rounded n values are:
#' ## 116, 52, 29, 14, 8 and 5
#' ## according to PASS manual, p. 595-5
#' D <- c(10, 15, 20, 30, 40, 50)
#' y <- sapply(D, function(delta) {
#'  power.williams.test(power = 0.9, k = 4, delta = delta, sd = 25)$n
#'  })
#' ceiling(y)
#'
#' \dontrun{
#'  ## compare with power.t.test
#'  ## and bonferroni correction
#'  power.t.test(power = 0.9, delta = 50, sd = 25,
#'  sig.level = 0.05 / 4, alternative = "one.sided")
#' }
#'
#' @export
power.williams.test <- function(n = NULL,
                                k,
                                delta,
                                sd = 1,
                                power = NULL,
                                ...) {

  sig.level = 0.05
  delta <- abs(delta)

  if (all(is.null(n), is.null(power))) {
    stop("One of n or power must be specified.")
  }

  ## check k number of treatment groups
  if (k < 2 | k > 10) {
    stop(sQuote(k), " must be in the interval, 2 <= k <= 10!")
  }

  if (is.null(n)) {

    ## from PASS Manual p. 595-2
    Beta <- 1 - power
    TKAlpha <- getTkalpha(k)
    #' @importFrom stats qnorm
    n <- 2 * sd^2 * (TKAlpha + qnorm(Beta, lower.tail = FALSE))^2 / delta ^2


  } else {
    # go for power
    power <- pwr.fn(n, delta, sd, k)

  }

  METHOD <-
    "William's test power calculation"

  NOTE <- "n is number in *each* group"

  ## add a print method
  out <- structure(
    list(
      n = n,
      k = k,
      delta = delta,
      sd = sd,
      sig.level = sig.level,
      power = power,
      method = METHOD,
      note = NOTE
    ),
    class = "power.htest"
  )

  out
}


getTkalpha <- function(k) {
  ## Williams 1971,
  ## alpha = 0.05, df = Inf
  Tcrit <- c(1.645,
             1.716,
             1.739,
             1.750,
             1.756,
             1.760,
             1.763,
             1.765,
             1.767,
             1.768)
  Tcrit[k]
}


## function to obtain power
pwr.fn <- function(n, delta, sd, k) {
  ## assume infinte df
  tkalpha <- getTkalpha(k = k)
  z <- tkalpha - (delta / (sd * sqrt(2 / n)))

  #' @importFrom stats pnorm
  power <- 1 - pnorm(z)
  power
}
