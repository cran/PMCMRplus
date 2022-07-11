## frdHouseTest.R
## Part of the R package: PMCMR
##
## Copyright (C)  2022 Thorsten Pohlert
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
#' @rdname frdHouseTest
#' @title House Test
#'
#' @description
#' Performs House nonparametric equivalent of William's test
#' for contrasting increasing dose levels of a treatment in
#' a complete randomized block design.
#'
#'
#' @details
#' House test is a non-parametric step-down trend test for testing several treatment levels
#' with a zero control. Let there be \eqn{k} groups including the control and let
#' the zero dose level be indicated with \eqn{i = 0} and the highest
#' dose level with \eqn{i = m}, then the following \code{m = k - 1} hypotheses are tested:
#'
#' \deqn{
#' \begin{array}{ll}
#' \mathrm{H}_{m}: \theta_0 = \theta_1 = \ldots = \theta_m, & \mathrm{A}_{m} = \theta_0 \le \theta_1 \le \ldots \theta_m, \theta_0 < \theta_m \\
#' \mathrm{H}_{m-1}: \theta_0 = \theta_1 = \ldots = \theta_{m-1}, & \mathrm{A}_{m-1} = \theta_0 \le \theta_1 \le \ldots \theta_{m-1}, \theta_0 < \theta_{m-1} \\
#' \vdots & \vdots \\
#' \mathrm{H}_{1}: \theta_0 = \theta_1, & \mathrm{A}_{1} = \theta_0 < \theta_1\\
#' \end{array}
#' }
#'
#' Let \eqn{Y_{ij} ~ (1 \leq i \leq n, 0 \leq j \leq k)} be a i.i.d. random variable
#' of at least ordinal scale. Further, let \eqn{\bar{R}_0,~\bar{R}_1, \ldots,~\bar{R}_k}
#' be Friedman's average ranks and set \eqn{\bar{R}_0^*, \leq \ldots \leq \bar{R}_k^*}
#' to be its isotonic regression estimators under the order restriction
#' \eqn{\theta_0 \leq \ldots \leq \theta_k}.
#'
#' The statistics is
#' \deqn{
#' T_j = \left(\bar{R}_j^* - \bar{R}_0 \right)~ \left[ \left(V_j - H_j \right)
#' \left(2 / n \right) \right]^{-1/2} \qquad (1 \leq j \leq k),
#' }
#'
#' with
#' \deqn{
#' V_j = \left(j + 1\right) ~ \left(j + 2 \right) / 12
#' }
#'
#' and
#' \deqn{
#' H_j = \left(t^3 - t \right) / \left(12 j n \right),
#' }
#'
#' where \eqn{t} is the number of tied ranks.
#'
#' The critical \eqn{t'_{i,v,\alpha}}-values
#' as given in the tables of Williams (1972) for \eqn{\alpha = 0.05} (one-sided)
#' are looked up according to the degree of freedoms (\eqn{v = \infty}) and the order number of the
#' dose level (\eqn{j}).
#'
#' For the comparison of the first dose level \eqn{(j = 1)} with the control, the critical
#' z-value from the standard normal distribution is used (\code{\link[stats]{Normal}}).
#'
#' @references
#' Chen, Y.-I., 1999. Rank-Based Tests for Dose Finding in
#' Nonmonotonic Dose–Response Settings.
#' *Biometrics* **55**, 1258--1262. \doi{10.1111/j.0006-341X.1999.01258.x}
#'
#' House, D.E., 1986. A Nonparametric Version of Williams’ Test for
#' Randomized Block Design. *Biometrics* **42**, 187--190.
#'
#' @concept friedmanranks
#' @keywords htest nonparametric
#'
#' @example examples/frdManyOne.R
#'
#' @seealso
#' \code{\link{friedmanTest}}, \code{\link[stats]{friedman.test}},
#' \code{\link{frdManyOneExactTest}}, \code{\link{frdManyOneDemsarTest}}
#'
#' @template class-PMCMR
#' @export
frdHouseTest <- function(y, ...) UseMethod("frdHouseTest")

#' @rdname frdHouseTest
#' @method frdHouseTest default
#' @aliases frdHouseTest.default
#' @param alternative the alternative hypothesis. Defaults to \code{greater}.
#' @template two-way-parms
#' @export
frdHouseTest.default <-
    function(y,
             groups,
             blocks,
             alternative = c("greater", "less"),
             ...)
{


    ## Check arguments
    alternative <- match.arg(alternative)
    if (alternative == "less") {
      y <- -y
    }

    ## Friedman-type ranking
    ans <- frdRanks(y)
    Rij <- ans$r

    ## number of levels including 0-control
    k <- ncol(Rij)
    if (k-1 > 10)
      stop("Critical t-values are only available for up to 10 dose levels.")

    ## sample size per group
    n <- nrow(Rij)

    ## mean ranks per group
    Rj <- as.vector(
      colMeans(Rij, na.rm = TRUE)
    )
    names(Rj) <- NULL

    ## ties
    ties <- sapply(seq_len(n),
                   function(i){
                     t <- table(Rij[i,])
                     sum(t - 1)
                   })
    t <- sum(ties)

    ##
    ## call to own pava
    Rjiso <- .Fortran(
      "pava",
      y = as.double(Rj),         # vector
      w = as.double(rep(n, k)),  # vector
      kt = integer(k),
      n = as.integer(k)
    )$y

    ## zero control
    Rj0 <- Rj[1]
    Rjiso <- Rjiso[-1]

    Tj <- sapply(1:(k-1),
                 function(j) {

                    Vj <- (j + 1) * (j + 2) / 12 -
                      (t^3 - t) / (12 * j * n)

                   (Rjiso[j] - Rj0) /
                     sqrt(Vj * 2 / n)
                 }
    )

    ## critical t'-values for df = Inf
    Tkdf <- getTkalpha(1:(k-1))


    ## Create output matrices
    STAT <- cbind(ctr = Tj)
    row.names(STAT) <- sapply(1:(k - 1), function(i)
      paste0("mu", i))
    STATCRIT <- cbind(ctr = Tkdf)
    row.names(STATCRIT) <- row.names(STAT)

    DAT <- ans$inDF
    names(DAT) <- c("y", "groups", "blocks")

    ## re-scale original data
    if (alternative == "less") {
      DAT$y <- -DAT$y
    }

    METHOD <- c("House (Williams) test on Friedman-Ranks")
    parameter <- Inf
    names(parameter) <- "df"
    ans <- list(
      method = METHOD,
      data.name = ans$DNAME,
      crit.value = STATCRIT,
      statistic = STAT,
      parameter = parameter,
      alternative = alternative,
      dist = "t\'",
      model = DAT
    )
    class(ans) <-"osrt"
    return(ans)
}
