# williamsTest.R
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

#' @name williamsTest
#' @title Williams Trend Test
#' @description
#' Performs Williams' test for contrasting increasing (decreasing) dose levels of a treatment.
#' @details
#' Williams' test is a step-down trend test for testing several treatment levels
#' with a zero control in a one-factorial design with normally distributed
#' errors of homogeneous variance. Let there be \eqn{k} groups including the control and let
#' the zero dose level be indicated with \eqn{i = 0} and the treatment
#' levels indicated as \eqn{1 \le i \le m}, then the following \eqn{m = k - 1} hypotheses are tested:
#'
#' \deqn{
#' \begin{array}{ll}
#' \mathrm{H}_{m}: \bar{x}_0 = m_1 = \ldots = m_m, & \mathrm{A}_{m}: \bar{x}_0 \ge m_1 \ge \ldots m_m, \bar{x}_0 < m_m \\
#' \mathrm{H}_{m-1}: \bar{x}_0 = m_1 = \ldots = m_{m-1}, & \mathrm{A}_{m-1}: \bar{x}_0 \ge m_1 \ge \ldots m_{m-1}, \bar{x}_0 < m_{m-1} \\
#' \vdots & \vdots \\
#' \mathrm{H}_{1}: \bar{x}_0 = m_1, & \mathrm{A}_{1}: \bar{x}_0 < m_1,\\
#' \end{array}
#' }
#'
#' where \eqn{m_i} denotes the isotonic mean of the \eqn{i}th dose level group.
#' The procedure starts from the highest dose level (\eqn{m}) to the the lowest dose level (\eqn{1}) and
#' stops at the first non-significant test. The consequent lowest effect dose
#' is the treatment level of the previous test number.
#'
#' The function does not return p-values. Instead the critical t-values
#' as given in the tables of Williams (1972) for \eqn{\alpha = 0.05} (one-sided)
#' are looked up according to the degree of freedoms (\eqn{v}) and the order number of the
#' dose level (\eqn{i}) and (potentially) modified according to the given extrapolation
#' coefficient \eqn{\beta}.
#'
#' Non tabulated values are linearly interpolated as recommended by Williams (1972).
#' The function \code{\link[stats]{approx}} is used.
#'
#' For the comparison of the first dose level (i = 1) with the control, the critical t-value
#' from the Student t distribution is used (\code{\link[stats]{TDist}}).
#'
#' @note
#' In the current implementation, only tests on the level of \eqn{\alpha = 0.05}
#' can be performed. The included extrapolation function assumes either
#' a balanced design, or designs, where the number of replicates in the control excdeeds the number of replicates
#' in the treatment levels. A warning message appears, if the following
#' condition is not met, \eqn{1 \le n_0 / n_i \le 6} for \eqn{1 \le i \le m}.
#'
#' @return
#' A list with class \code{"williamsTest"} containing the following components:
#' \describe{
#'  \item{method}{a character string indicating what type of test was performed.}
#'  \item{data.name}{a character string giving the name(s) of the data.}
#'  \item{statistic}{lower-triangle matrix of the estimated
#' quantiles of the pairwise test statistics.}
#'  \item{t.value}{lower-triangle matrix of the critical t\'-values for \eqn{\alpha = 0.05}.}
#'  \item{df.residual}{the degree of freedom}
#'  \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{model}{a data frame of the input data.}
#' \item{dist}{a string that denotes the test distribution.}
#' }
#' There are print and summary methods available.
#'
#' @section Source:
#' The source code for the application of the pool adjacent violators
#' theorem to calculate the isotonic means
#' was taken from the file \code{"pava.f"}, which is included in the
#' package \pkg{Iso}:
#'
#'  Rolf Turner (2015). Iso: Functions to Perform Isotonic Regression. R
#'  package version 0.0-17. \url{https://CRAN.R-project.org/package=Iso}.
#'
#' The file \code{pava.f} is a Ratfor modification of Algorithm AS 206.1:
#'
#' Bril, G., Dykstra, R., Pillers, C., Robertson, T. (1984)
#' Statistical Algorithms: Algorithm AS 206: Isotonic
#' Regression in Two Independent Variables, \emph{Appl. Statist.},
#' 34, 352--357.
#'
#'  The Algorith AS 206 is available from StatLib
#' \url{http://lib.stat.cmu.edu/apstat/}. The Royal Statistical Society
#' holds the copyright to these routines,
#' but has given its permission for their distribution provided that
#' no fee is charged.
#'
#' @seealso
#' \code{\link[stats]{TDist}}, \code{\link[stats]{approx}}, \code{\link{print.williams}},
#' \code{\link{summary.williams}}
#'
#' @useDynLib 'PMCMRplus', .registration = TRUE, .fixes = "F_"
#' @references
#' Williams, D. A. (1971) A test for differences between treatment means
#'   when several dose levels are compared with a zero dose control,
#'   \emph{Biometrics} \bold{27}, 103--117.
#'
#' Williams, D. A. (1972) The comparison of several dose levels with a zero
#'   dose control, \emph{Biometrics} \bold{28}, 519--531.
#' @keywords htest
#' @importFrom stats qt approx var
#' @examples
#' ## Example from Sachs (1997, p. 402)
#' x <- c(106, 114, 116, 127, 145,
#' 110, 125, 143, 148, 151,
#' 136, 139, 149, 160, 174)
#' g <- gl(3,5)
#' levels(g) <- c("0", "I", "II")
#'
#' ## Williams Test
#' williamsTest(x ~ g)
#'
#' @export
williamsTest <- function(x, ...)
  UseMethod("williamsTest")

#' @rdname williamsTest
#' @aliases williamsTest.default
#' @method williamsTest default
#' @template one-way-parms
#' @param alternative the alternative hypothesis. Defaults to \code{greater}
#' @export
williamsTest.default <-
  function(x,
           g,
           alternative = c("greater", "less"),
           ...)
  {
    ## taken from stats::williamsTest
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
    xold <- x
    if (alternative == "less") {
      x <- -x
    }

    xi <- tapply(x, g, mean, na.rm = T)
    ni <- tapply(x, g, length)


    k <- nlevels(g)
    kk <- k - 1
    if (kk > 10) stop("Critical t-values are only available for up to 10 dose levels.")
    N <- sum(ni)

    ## pooled within-group variance
    df <- N - k
    s2i <- tapply(x, g, var)
    s2in <- 1 / df * sum(s2i * (ni - 1))

    ## call to own pava
    xiiso <- .Fortran(
      "pava",
      y = as.double(xi),
      w = as.double(ni),
      kt = integer(k),
      n = as.integer(k)
    )$y


    mui <- rep(NA, k)

## This is for alternative greater
##    if (alternative == "greater") {
      for (i in 1:k) {
        v <- k
        tmp <- rep(NA, length(1:i))
        for (u in 1:i) {
          j <- u
          tmp01 <-
            sapply(i:k, function(v)
              sum(ni[j:v] * xiiso[j:v]) /
                sum(ni[j:v]))
          tmp[u] <- min(tmp01)
        }
        mui[i] <- max(tmp, na.rm = TRUE)
      }

      if (alternative == "greater") {
      Tk <- sapply(2:k, function(i) {
        (mui[i] - xi[1]) / sqrt((s2in / ni[i] + s2in / ni[1]))
      })
      } else {
        Tk <- sapply(2:k, function(i) {
          (xi[1] - mui[i] ) / sqrt((s2in / ni[i] + s2in / ni[1]))
        })
      }

    ## Extrapolation function, see Williams, 1972, p. 530
    extrapolFN <- function(Tki, beta, r, c) {
      out <- Tki - 1E-2 * beta * (1 - r / c)
      return(out)
    }

    ## load critical t values (are in sysdata.rda)
    c <- ni[1]
    r <- ni[2:k]

    ## check ratio o
    o <- c / r
    for (i in 1:kk) {
      if (o[i] < 1 | o[i] > 6) {
        warning(paste0("Ratio n0 / n", i, " is ", o[i], " and outside the range.\n
                       Test results may not be accurate."))
      }
    }


    nrows <- nrow(TabCrit$williams.tk005)
    Tkdf <- numeric(kk)
    dft <- as.numeric(rownames(TabCrit$williams.tk005))
    xx <- c(2:6,8,10)
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

    ## Critical t-value from Students t-distribution for i = 1
    Tkdf[1] <- qt(0.05, df = df, lower.tail = FALSE)

    ## Check alternative
    Tkdf <- sapply(Tkdf, function(i) ifelse(alternative == "greater", i, -i))

    ## Create output matrices
    STAT <- cbind(ctr = Tk)
    row.names(STAT) <- sapply(1:(k - 1), function(i)
      paste0("mu", i))
    STATCRIT <- cbind(ctr = Tkdf)
    row.names(STATCRIT) <- row.names(STAT)

    METHOD <- paste("Williams trend test")
    ans <- list(
      method = METHOD,
      data.name = DNAME,
      t.value = STATCRIT,
      statistic = STAT,
      df.residual = df,
      alternative = alternative,
      dist = "t\'",
      model = data.frame(x=xold, g = g)
    )
    class(ans) <- "williams"
    ans
  }

#' @rdname williamsTest
#' @aliases williamsTest.formula
#' @method williamsTest formula
#' @template one-way-formula
#' @export
williamsTest.formula <-
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
    y <- do.call("williamsTest", c(as.list(mf), alternative = alternative))
    y$data.name <- DNAME
    y
  }
