# shirleyWilliamsTest.R
# Part of the R package: PMCMR
#
##  Copyright (C) 2017, 2018 Thorsten Pohlert
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

#' @name shirleyWilliamsTest
#' @title Shirley-Williams Test
#' @description
#' Performs Shirley's nonparametric equivalent of William's test
#' for contrasting increasing dose levels of a treatment.
#'
#' @details
#' The Shirley-William test is a non-parametric step-down trend test for testing several treatment levels
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
#' Let \eqn{R_{ij}} be the rank of \eqn{X_{ij}},
#' where \eqn{X_{ij}} is jointly ranked
#' from \eqn{\left\{1, 2, \ldots, N \right\}, ~~ N = \sum_{i=1}^k n_i},
#' then the test statistic is
#'
#' \deqn{
#'   t_{i} = \frac{\max_{1 \le u \le i} \left(\sum_{j=u}^i n_j \bar{R}_j / \sum_{j=u}^i n_j \right) - \bar{R}_0}
#' {\sigma_{R_i} \sqrt{1/n_i + 1/n_0}},
#' }{%
#'  SEE PDF
#' }
#'
#' with expected variance of
#' \deqn{
#' \sigma_{R_i}^2 = N_i \left(N_i + 1 \right) / 12 - T_i,
#' }{%
#'  SEE PDF
#' }
#'
#' where \eqn{N_i = n_0 + n_1 + n_2 + \ldots + n_i} and
#' \eqn{T_i} the ties for the \eqn{i}-th comparison is given by
#'
#' \deqn{
#'  T_i = \sum_{j=1}^i \frac{t_j^3 - t_j}{12 \left(N_i - 1\right)}.
#' }{%
#'  SEE PDF
#' }
#'
#' The procedure starts from the highest dose level (\eqn{m}) to the the lowest dose level (\eqn{1}) and
#' stops at the first non-significant test. The consequent lowest effect dose
#' is the treatment level of the previous test number. This function has
#' included the modifications as recommended by Williams (1986), i.e.
#' the data are re-ranked for each of the \eqn{i}-th comparison.
#'
#' If \code{method = "look-up"} is selected, the function does not return p-values.
#' Instead the critical \eqn{t'_{i,v,\alpha}}-values
#' as given in the tables of Williams (1972) for \eqn{\alpha = 0.05} (one-sided)
#' are looked up according to the degree of freedoms (\eqn{v = \infty}) and the order number of the
#' dose level (\eqn{i}) and (potentially) modified according to the given extrapolation
#' coefficient \eqn{\beta}.
#'
#' Non tabulated values are linearly interpolated with the function
#' \code{\link[stats]{approx}}.
#'
#' For the comparison of the first dose level (i = 1) with the control, the critical
#' z-value from the standard normal distribution is used (\code{\link[stats]{Normal}}).
#'
#' If \code{method = "boot"}, the p-values are estimated through an assymptotic
#' boot-strap method. The p-values for H\eqn{_1}
#' are calculated from the t distribution with infinite degree of freedom.
#'
#' @note
#' For \code{method = "look-up"}, only tests on the level of \eqn{\alpha = 0.05}
#' can be performed for alternative hypotheses less or greater.
#'
#' For \code{method = "boot"} only the alternative \code{"two.sided"} can be calculated.
#' One may increase the number of permutations to e.g. \code{nperm = 10000}
#' in order to get more precise p-values. However, this will be on the expense of
#' computational time.
#'
#' @references
#' Shirley, E., (1977) Nonparametric Equivalent of Williams Test for Contrasting Increasing
#' Dose Levels of a Treatment, \emph{Biometrics} \bold{33}, 386--389.
#'
#' Williams, D. A. (1986) Note on Shirley's nonparametric test for comparing
#' several dose levels with a zero-dose control, \emph{Biometrics} \bold{42}, 183--186.
#'
#' @return
#' Either a list with class \code{"osrt"} or a list with class \code{"PMCMR"}.
#'
#' @template returnOsrt
#' @template class-PMCMR
#'
#' @seealso
#' \code{\link{williamsTest}}
#'
#' @keywords htest nonparametric
#' @concept trendtest
#'
#' @example examples/shirleyEx.R
#'
#'
#' @export
shirleyWilliamsTest <-
  function(x, ...)
    UseMethod("shirleyWilliamsTest")

#' @rdname shirleyWilliamsTest
#' @aliases shirleyWilliamsTest.default
#' @method shirleyWilliamsTest default
#' @template one-way-parms
#' @param alternative the alternative hypothesis. Defaults to \code{two.sided}
#' @param method a character string specifying the test statistic to use.
#' Defaults to \code{"look-up"} that uses published Table values of Williams (1972).
#' @param nperm number of permutations for the asymptotic permutation test.
#' Defaults to \code{1000}. Ignored, if \code{method = "look-up"}.
#' @importFrom stats qnorm approx
#' @export
shirleyWilliamsTest.default <-
  function(x,
           g,
           alternative = c("two.sided", "greater", "less"),
           method = c("look-up", "boot"),
           nperm = 1E4,
           ...) {
    ## taken from stats::kruskal.test

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
      nperm <- x$nperm
      method <- x$method
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

    alternative <- match.arg(alternative)
    method <- match.arg(method)
    if (method == "boot" & alternative != "two.sided") {
      alternative <- "two.sided"
      warning("As method is 'boot', alternative was set to 'two.sided'.")

    } else if (alternative == "two.sided" & method != "boot") {
      warning("As alternative is 'two.sided', method was set to 'boot'.")
      method <- "boot"
    }

    xold <- x

    if (method == "look-up" & alternative == "less") {
      x <- -x
    }

    nj <- tapply(x, g, length)
    k <- nlevels(g)
    kk <- k - 1
    if (kk > 10)
      stop("Critical t-values are only available for up to 10 dose levels.")
    N <- sum(nj)

    ## comparisons
    compfn <- function(x, ix, g, nj) {
      k <- length(nj)
      ti <- rep(NA, k)
      x <- x[ix]
      for (i in 2:k) {
        N <- sum(nj[1:i])
        r <- rank(x[1:N])
        gg <- g[1:N]
        Rj <- tapply(r, gg, mean)
        t <- table(r)
        names(t) <- NULL
        T <- sum((t ^ 3 - t) / (12 * (N - 1)))
        Vi <- N * (N + 1) / 12 - T
        u <- 2:i
        j <- u
        enum <- sapply(j, function(j)
          sum(nj[j:i] * Rj[j:i]))
        denom <- sapply(j, function(j)
          sum(nj[j:i]))

        ti[i] <- (max(enum / denom) - Rj[1]) /
          sqrt(Vi * (1 / nj[i] + 1 / nj[1]))
      }
      return(ti[2:k])
    }

    l <- 1:N
    STATISTIC <- compfn(x, l, g, nj)

    if (method == "boot") {
      ## permutation
      mt <- matrix(NA, ncol = (k - 1), nrow = nperm)
      for (i in 1:nperm) {
        ix <- sample(l)
        mt[i, ] <- compfn(x, ix, g, nj)
      }

      ## pvalues
      PVAL <- sapply(1:(k - 1), function(j) {
        p <- (sum(mt[, j] <= -abs(STATISTIC[j])) +
                (sum(mt[, j] >= abs(STATISTIC[j])))) / nperm
        p
      })

      ## exact p-value from student t distribution
      PVAL[1] <-
        2 * min(0.5, pt(STATISTIC[1], df = Inf, lower.tail = FALSE))

      STAT <- cbind(ctr = STATISTIC)
      row.names(STAT) <- sapply(1:(k - 1), function(i)
        paste0("mu", i))
      P <- cbind(ctr = PVAL)
      row.names(P) <- row.names(STAT)

      DAT <- data.frame(x, g)
      METH <- c("Shirley-Williams test")
      ans <- list(
        statistic = STAT,
        p.value = P,
        data = DAT,
        method = METH,
        data.name = DNAME,
        alternative = "two.sided",
        dist = "t",
        p.adjust.method = "boot"
      )
      class(ans) <- "PMCMR"
      return(ans)

    } else {
      ## Extrapolation function, see Williams, 1972, p. 530
      extrapolFN <- function(Tki, beta, r, c) {
        out <- Tki - 1E-2 * beta * (1 - r / c)
        return(out)
      }

      ## load critical t values (are in sysdata.rda)
      df <- 1E6  # Originally Inf
      c <- nj[1]
      r <- nj[2:k]
      nrows <- nrow(TabCrit$williams.tk005)
      Tkdf <- numeric(kk)
      dft <- as.numeric(rownames(TabCrit$williams.tk005))
      xx <- c(2:6, 8, 10)
      for (i in 2:kk) {
        if (i <= 6 | i == 8 | i == 10) {
          ## here only df needs to be interpolated
          yt <- TabCrit$williams.tk005[, paste0(i)]
          yb <- TabCrit$williams.beta005[, paste0(i)]

        } else {
          yt <- sapply(1:nrows, function(j) {
            approx(x = xx,
                   y = TabCrit$williams.tk005[j,],
                   xout = i)$y
          })
          yb <- sapply(1:nrows, function(j) {
            approx(x = xx,
                   y = TabCrit$williams.beta005[j,],
                   xout = i)$y
          })
        }

        tt <- approx(x = dft, y = yt, xout = df)$y
        tb <- approx(x = dft, y = yb, xout = df)$y

        Tkdf[i] <- extrapolFN(tt, tb, r[i], c)
      }

      ## Critical t_Inf value = z value for i = 1
      Tkdf[1] <- qnorm(0.05, lower.tail = FALSE) ## greater

  #    if (alternative == "less") {
  #      STATISTIC <- -1 * STATISTIC
  #      Tkdf <- -1 * Tkdf
  #    }

      ## Create output matrices
      STAT <- cbind(ctr = STATISTIC)
      row.names(STAT) <- sapply(1:(k - 1), function(i)
        paste0("mu", i))
      STATCRIT <- cbind(ctr = Tkdf)
      row.names(STATCRIT) <- row.names(STAT)

      DAT <- data.frame(xold, g)
      METHOD <- c("Shirley-Williams test")
      parameter <- Inf
      names(parameter) <- "df"
      ans <- list(
        method = METHOD,
        data.name = DNAME,
        crit.value = STATCRIT,
        statistic = STAT,
        parameter = parameter,
        alternative = alternative,
        dist = "t\'",
        model = DAT
      )
      class(ans) <-"osrt"
   #   class(ans) <- "williams"
      return(ans)
    }
  }

#' @rdname shirleyWilliamsTest
#' @aliases shirleyWilliamsTest.formula
#' @method shirleyWilliamsTest formula
#' @template one-way-formula
#' @export
shirleyWilliamsTest.formula <-
  function(formula,
           data,
           subset,
           na.action,
           alternative = c("two.sided", "greater", "less"),
           method = c("look-up", "boot"),
           nperm = 1E4,
           ...)
  {
    mf <- match.call(expand.dots = FALSE)
    m <-
      match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(stats::model.frame)

    if (missing(formula) || (length(formula) != 3L))
      stop("'formula' missing or incorrect")
    mf <- eval(mf, parent.frame())
    if (length(mf) > 2L)
      stop("'formula' should be of the form response ~ group")
    DNAME <- paste(names(mf), collapse = " by ")
    alternative <- match.arg(alternative)
    method <- match.arg(method)
    names(mf) <- NULL
    y <-
      do.call(
        "shirleyWilliamsTest",
        c(
          as.list(mf),
          alternative = alternative,
          method = method,
          nperm = nperm
        )
      )
    y$data.name <- DNAME
    y
  }
