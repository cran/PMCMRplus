## chenJanTest.R
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
##
#' @name chenJanTest
#' @title Chen and Jan Many-to-One Comparisons Test
#' @description
#' Performs Chen and Jan nonparametric test for contrasting increasing
#' (decreasing) dose levels of a treatment in a randomized block design.
#'
#' @details
#' Chen's test is a non-parametric step-down trend test for
#' testing several treatment levels with a zero control. Let
#' there be \eqn{k} groups including the control and let
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
#' Let \eqn{Y_{ij1}, Y_{ij2}, \ldots, Y_{ijn_{ij}}}
#' \eqn{(i = 1, 2, \dots, b, j = 0, 1, \ldots, k ~ \mathrm{and} ~ n_{ij} \geq 1)} be
#' a i.i.d. random variable of at least ordinal scale. Further,the zero dose
#' control is indicated with \eqn{j = 0}.
#'
#' The Mann-Whittney statistic is
#'
#' \deqn{
#' T_{ij} = \sum_{u=0}^{j-1} \sum_{s=1}^{n_{ij}}
#' \sum_{r=1}^{n_{iu}} I(Y_{ijs} - Y_{iur}),
#' \qquad i = 1, 2, \ldots, b, ~ j = 1, 2, \ldots, k,
#' }
#'
#' where where the indicator function returns \eqn{I(a) = 1, ~ \mathrm{if}~ a > 0, 0.5 ~ \mathrm{if} a = 0}
#' otherwise \eqn{0}.
#'
#' Let
#' \deqn{
#' N_{ij} = \sum_{s=0}^j n_{is} \qquad  i = 1, 2, \ldots, b, ~ j = 1, 2, \ldots, k,
#' }
#'
#' and
#'
#' \deqn{
#' T_j = \sum_{i=1}^b T_{ij} \qquad j = 1, 2, \ldots, k.
#' }
#'
#' The mean and variance of \eqn{T_j} are
#'
#' \deqn{
#' \mu(T_j) = \sum_{i=1}^b n_{ij} ~ N_{ij-1} / 2 \qquad \mathrm{and}
#' }
#'
#' \deqn{
#'  \sigma(T_j) = \sum_{i=1}^b n_{ij} ~ N_{ij-1} \left[
#'  \left(N_{ij} + 1\right) - \sum_{u=1}^{g_i}
#'  \left(t_u^3 - t_u \right) /
#'  \left\{N_{ij} \left(N_{ij} - 1\right) \right\} \right]/ 2,
#' }
#'
#' with \eqn{g_i} the number of ties in the \eqn{i}th block and
#' \eqn{t_u} the size of the tied group \eqn{u}.
#'
#' The test statistic \eqn{T_j^*} is asymptotically multivariate normal
#' distributed.
#'
#' \deqn{
#' T_j^* = \frac{T_j - \mu(T_j)}{\sigma(T_j)}
#' }
#'
#' If \code{p.adjust.method =  "single-step"} than the p-values
#' are calculated with the probability function of the multivariate
#' normal distribution with \eqn{\Sigma = I_k}. Otherwise
#' the standard normal distribution is used to calculate
#' p-values and any method as available
#' by \code{\link{p.adjust}} or by the step-down procedure as proposed
#' by Chen (1999), if \code{p.adjust.method = "SD1"} can be used
#' to account for \eqn{\alpha}-error inflation.
#'
#' @template class-PMCMR
#'
#' @examples
#' ## Example from Chen and Jan (2002, p. 306)
#' ## MED is at dose level 2 (0.5 ppm SO2)
#' y <- c(0.2, 6.2, 0.3, 0.3, 4.9, 1.8, 3.9, 2, 0.3, 2.5, 5.4, 2.3, 12.7,
#' -0.2, 2.1, 6, 1.8, 3.9, 1.1, 3.8, 2.5, 1.3, -0.8, 13.1, 1.1,
#' 12.8, 18.2, 3.4, 13.5, 4.4, 6.1, 2.8, 4, 10.6, 9, 4.2, 6.7, 35,
#' 9, 12.9, 2, 7.1, 1.5, 10.6)
#' groups <- gl(4,11, labels = c("0", "0.25", "0.5", "1.0"))
#' blocks <- structure(rep(1:11, 4), class = "factor",
#' levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"))
#'
#' summary(chenJanTest(y, groups, blocks, alternative = "greater"))
#' summary(chenJanTest(y, groups, blocks, alternative = "greater", p.adjust = "SD1"))
#'
#' @references
#' Chen, Y.I., Jan, S.L., 2002. Nonparametric Identification of
#' the Minimum Effective Dose for Randomized Block Designs.
#' *Commun Stat-Simul Comput* **31**, 301--312.
#'
#' @seealso
#' \code{\link[stats]{Normal}} \code{\link[mvtnorm]{pmvnorm}}
#' @concept wilcoxonranks
#' @keywords htest nonparametric
#' @export
chenJanTest <- function(y, ...)
  UseMethod("chenJanTest")

#' @rdname chenJanTest
#' @method chenJanTest default
#' @aliases chenJanTest.default
#' @param alternative the alternative hypothesis. Defaults to \code{greater}.
#' @template two-way-parms
#' @param p.adjust.method method for adjusting p values
#' (see \code{\link{p.adjust}})
#' @importFrom stats pnorm
#' @importFrom stats p.adjust
#' @importFrom stats p.adjust.methods
#' @importFrom stats complete.cases
#' @importFrom mvtnorm pmvnorm
#' @export
chenJanTest.default <-
  function(y,
           groups,
           blocks,
           alternative = c("greater", "less"),
           p.adjust.method =  c("single-step", "SD1", p.adjust.methods),
           ...) {
    DNAME <- deparse(substitute(y))
    if (is.matrix(y)) {
      groups <- factor(c(col(y)))
      blocks <- factor(c(row(y)))

      GRPNAMES <- colnames(y)
      ROWNAMES <- rownames(y)
    }
    else {
      if (anyNA(groups) || anyNA(blocks))
        stop("NA's are not allowed in 'groups' or 'blocks'")
      if (any(diff(c(
        length(y), length(groups), length(blocks)
      )) != 0L))
        stop("'y', 'groups' and 'blocks' must have the same length")
      DNAME <- paste0(DNAME,
                      ", ",
                      deparse(substitute(groups)),
                      " and ",
                      deparse(substitute(blocks)))
      # if (any(table(groups, blocks) != 1))
      #   stop("not an unreplicated complete block design")
      groups <- factor(groups)
      blocks <- factor(blocks)
      ## Need to ensure consistent order of observations within
      ## blocks.
      o <- order(groups, blocks)
      y <- y[o]
      groups <- groups[o]
      blocks <- blocks[o]

      GRPNAMES <- levels(groups)
      ROWNAMES <- levels(blocks)
    }


    # check arguments
    alternative <- match.arg(alternative)
    p.adjust.method <- match.arg(p.adjust.method)


    if (alternative == "less") {
      y <- -y
    }


    ## Indicator function
    I <- function(a) {
      if (a > 0) {
        ans <- 1
      } else if (a < 0) {
        ans <- 0
      } else {
        ans <- 1 / 2
      }
      ans
    }


    ## dose levels (1 <= j <= k) )
    k <- nlevels(groups)
    glev <- levels(groups)

    ## blocks
    b <- nlevels(blocks)
    blev <- levels(blocks)

    ## treatment levels (m = k-1)
    m <- k - 1
    # Chen and Jan 2002, p.304
    Tij <- matrix(0, nrow = b, ncol = m)
    nij <- matrix(0, nrow = b, ncol = k)

    for (i in seq_len(b)) {
      ## 0-dose control
      nij[i, 1] <- length(y[blocks == blev[i] & groups == glev[1]])
      for (j in seq_len(m)) {
        Yij <- y[blocks == blev[i] &  groups == glev[j + 1]]
        nij[i, j + 1] <- length(Yij)
        usum <- 0
        for (u in 0:(j - 1)) {
          Yiu <- y[blocks == blev[i] &  groups == glev[u + 1]]
          niu <- length(Yiu)
          ssum <- 0
          for (s in 1:nij[i, j + 1]) {
            rsum <- 0
            for (r in 1:niu) {
              rsum <- rsum + I(Yij[s] - Yiu[r])
            }
            ssum <- ssum + rsum
          }
          usum <- usum + ssum
        }
        Tij[i, j] <- usum
      }
    }

    ## test statistic per treatment
    T <- as.vector(colSums(Tij))

    ##
    Nij <- matrix(0, nrow = b, ncol = k)
    for (i in seq_len(b)) {
      Nij[i, ] <- cumsum(nij[i, ])
    }

    ##mean
    mu <- rep(0, m)
    for (j in seq_len(m)) {
      for (i in seq_len(b)) {
        mu[j] <- mu[j] + nij[i, j + 1] * Nij[i, j] / 2
      }
    }

    ## variance without ties
    sigma2 <- rep(0, m)
    for (j in seq_len(m)) {
      for (i in seq_len(b)) {
        Yb <- b[blocks == blev[i]]
        ties <- sum((table(Yb) - 1) ^ 3 - (table(Yb) - 1))

        sigma2[j] <- sigma2[j] +
          nij[i, j + 1] * Nij[i, j] *
          ((Nij[i, j + 1] + 1) - ties /
             (Nij[i, j + 1] * (Nij[i, j + 1] - 1)))  / 12
      }
    }

    ##
    STATISTIC <- (T - mu) / sqrt(sigma2)

    ## p-values, one-sided
    if (p.adjust.method == "single-step") {
      SIGMA <- diag(m)
      pvalv <- sapply(STATISTIC, function(x)
        1 - pmvnorm(
          lower = -Inf,
          upper = rep(x, m),
          sigma = SIGMA
        ))
      padj <- pvalv

    } else if (p.adjust.method == "SD1") {
      pval <- pnorm(STATISTIC , lower.tail = FALSE)
      ## p-adjustment
      padj <- SD1p(pval)
    } else {
      pval <- pnorm(STATISTIC , lower.tail = FALSE)
      ## p-adjustment
      padj <- p.adjust(pval, method = p.adjust.method)
    }

    ##
    if (alternative == "less") {
      STATISTIC <- -STATISTIC
    }

    ## prepare output
    GRPNAMES <- c("crt", paste0("mu", 1:m))
    PSTAT <- cbind(STATISTIC)
    PVAL <- cbind(padj)
    colnames(PSTAT) <- GRPNAMES[1]
    colnames(PVAL) <- GRPNAMES[1]
    rownames(PSTAT) <- GRPNAMES[2:k]
    rownames(PVAL) <- GRPNAMES[2:k]

    DIST <- "T*"
    METHOD <- paste("Chen-Jan test\n",
                    "\tfor multiple comparisons with one control",
                    sep = "")
    MODEL <- data.frame(y, groups, blocks)

    ans <-
      list(
        method = METHOD,
        data.name = DNAME,
        p.value = PVAL,
        statistic = PSTAT,
        p.adjust.method = p.adjust.method,
        model = MODEL,
        dist = DIST,
        alternative = alternative
      )
    class(ans) <- "PMCMR"
    ans
  }
