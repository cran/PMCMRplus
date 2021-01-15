#  Part of the R package PMCMRplus
#
#  Copyright (C) 2021 T. Pohlert
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
# F(x, y) =

#' @name mrrTest
#' @title Madhava Rao-Raghunath Test for Testing Treatment vs. Control
#'
#' @description
#' The function has implemented the nonparametric test of
#' Madhava Rao and Raghunath (2016) for testing paired two-samples
#' for symmetry. The null hypothesis \eqn{H: F(x,y) = F(y,x)}
#' is tested against the alternative \eqn{A: F(x,y) \ne F(y,x)}{A: F(x,y) != F(y,x)}.
#'
#' @param x numeric vector of data values. Non-finite (e.g., infinite or missing) values will be omitted.
#' @param y an optional numeric vector of data values:
#'          as with x non-finite values will be omitted.
#' @param m numeric, optional integer number, whereas \eqn{n = k m} needs to be
#'          full filled.
#'
#' @details
#' Let \eqn{X_i} and \eqn{Y_i, ~ i \le n} denote
#' continuous variables that were observed
#' on the same \eqn{i}th test item (e.g. patient)
#' with \eqn{i = 1, \ldots n}. Let
#' \deqn{
#'   U_i = X_i + Y_i \qquad V_i = X_i - Y_i
#' }{%
#'  SEE PDF
#' }
#'
#' Let \eqn{U_{(i)}} be the \eqn{i}th order statistic,
#' \eqn{U_{(1)} \le U_{(2)} \le \ldots U_{(n)}} and \eqn{k} the
#' number of clusters, with the condition:
#'
#' \deqn{
#'  n = k ~ m.
#' }{%
#'  SEE PDF
#' }
#'
#' Further, let the divider denote \eqn{d_0 = -\infty}, \eqn{d_k = \infty}, and else
#' \deqn{
#'  d_j = \frac{ U_{(jm)} +  U_{(jm+1)} }{2}, ~ 1 \le j \le k -1
#' }{%
#'  SEE PDF
#' }
#'
#' The two counts are
#' \deqn{
#'  n_j^{+} = \left\{
#'            \begin{array}{lr}
#'            1 & \mathrm{if}~ d_{j-1} < u_i < d_j, v_i > 0 \\
#'            0 &
#'            \end{array}
#'            \right.
#' }{%
#'  SEE PDF
#' }
#'
#' and
#' \deqn{
#'  n_j^{-} = \left\{
#'            \begin{array}{lr}
#'            1 & \mathrm{if}~ d_{j-1} < u_i < d_j, v_i \le 0 \\
#'            0 &
#'            \end{array}
#'            \right.
#' }{%
#'  SEE PDF
#' }
#'
#' The test statistic is
#' \deqn{
#'  M = \sum_{j = 1}^k \frac{\left(n_j^{+} - n_j^{-}\right)^2}
#'                          {m}
#' }{%
#'  SEE PDF
#' }
#'
#' The exact p-values for \eqn{5 \le n \le 30} are taken from an
#' internal look-up table. The exact p-values were taken
#' from Table 7, Appendix B  of Madhava Rao and Raghunath (2016).
#'
#' If \code{m = NULL} the function uses \eqn{n = m} for
#' all prime numbers, otherwise it tries to find an value for
#' m in such a way, that for \eqn{k = n / m} all variables
#' are integer.
#'
#' @note
#' The function returns an error code if a value for \code{m}
#' is provided  that does not lead to an integer of the ratio
#' \eqn{k = n /m}.
#'
#' The function also returns an error code, if a tabulated
#' value for given \eqn{n}, \eqn{m} and calculated \eqn{M}
#' can not be found in the look-up table.
#'
#' @template class-htest
#'
#' @references
#' Madhava Rao, K.S., Ragunath, M. (2016) A Simple Nonparametric Test
#' for Testing Treatment Versus Control. \emph{J Stat Adv Theory Appl} \bold{16},
#' 133–162. \doi{10.18642/jsata_7100121717}
#'
#' @examples
#' ## Madhava Rao and Raghunath (2016), p. 151
#' ## Inulin clearance of living donors
#' ## and recipients of their kidneys
#' x <- c(61.4, 63.3, 63.7, 80.0, 77.3, 84.0, 105.0)
#' y <- c(70.8, 89.2, 65.8, 67.1, 87.3, 85.1, 88.1)
#' mrrTest(x, y)
#'
#' ## formula method
#' ## Student's Sleep Data
#' mrrTest(extra ~ group, data = sleep)
#'
#' @importFrom gmp isprime
#' @export
mrrTest <- function(x, ...)
  UseMethod("mrrTest")

#' @rdname mrrTest
#' @method mrrTest default
#' @aliases mrrTest.default
#' @export
mrrTest.default <-
  function(x, y = NULL, m = NULL, ...)
  {
    paired <- TRUE

    if (!is.numeric(x))
      stop("'x' must be numeric")
    if (!is.null(y)) {
      if (!is.numeric(y))
        stop("'y' must be numeric")
      DNAME <- paste(deparse(substitute(x)), "and",
                     deparse(substitute(y)))
      if (paired) {
        if (length(x) != length(y))
          stop("'x' and 'y' must have the same length")
        OK <- complete.cases(x, y)
        x <- x[OK]
        y <- y[OK]
      }

    } else {
      DNAME <- deparse(substitute(x))
      if (paired)
        stop("'y' is missing for paired test")
      x <- x[is.finite(x)]
    }

    n <- length(x)

    if (n < 5L) {
      stop("not enough (finite) 'x' observations")
    } else if (n > 30) {
      stop("current version has only exact p for n=30")
    }

    ## aux function
    is.ratiointeger <- function(x, y)
      ifelse(x %% y > 0, FALSE, TRUE)

    ## check for optional m
    if (is.null(m)) {
      ## check for prime number
      if (isprime(n) == 2) {
        m <- n
      } else {
        ## determine m
        if (n < 25) {
          div <- c(3:2)
        } else {
          div <- c(5, 3, 2)
        }

        for (i in seq_along(div)) {
          if (is.ratiointeger(n, div[i])) {
            m <- n / div[i]
            break
          }
        }
      }
    } else {
      ## check whether ratio is integer
      if (!is.ratiointeger(n, m)) {
        stop("Select m in such a way that n = k m")
      }
    }


    ## continue with k
    k <- n / m

    ## estimate statistic
    U <- x + y
    V <- x - y

    Us <- sort(U)
    dj <- numeric(k + 1)
    dj[1] <- -Inf
    dj[k + 1] <- Inf
    if (k >= 2) {
      for (jj in 2:k) {
        j <- jj - 1
        dj[jj] <- (Us[j * m] + Us[j * m + 1]) / 2
      }
    }
    ## count
    njp <- rep(0, k)
    njm <- rep(0, k)

    for (jj in 2:(k + 1)) {
      j <- jj - 1
      for (i in 1:n) {
        if (dj[jj - 1] < U[i] & U[i] < dj[jj]) {
          if (V[i] > 0) {
            njp[j] <- njp[j] + 1
          } else {
            njm[j] <- njm[j] + 1
          }
        }
      }
    }

    STATISTIC <- sum((njp - njm) ^ 2 / m)

    ## binomial distribution, exact probabilities
    ## see Madhava Rao and Raghunath (2016), p. 141
    ## values taken from Table 7, Appendix B
    #  Njp <- matrix(rep(seq(m, 0), k), nrow = m + 1, ncol=k)
    #  Njm <- matrix(rep(rev(Njp), k), nrow = m + 1, ncol = k)
    #  M <- rowSums((Njp - Njm)^2) #/ m

    Tab <- TabCrit$mrr.exact

    Tab$n <- as.character(Tab$n)
    Tab$m <- as.character(Tab$m)

    ok1 <- Tab$n %in% as.character(n)
    ok2 <- Tab$m %in% as.character(m)

    if (!any(ok1 & ok2)) {
      stop("Could not find M value for n = ", n , " and m =", m)
    }

    Tab <- Tab[ok1 & ok2, ]

    ok <- sapply(Tab$M, function(b)
      identical(all.equal(b , STATISTIC), TRUE))

    if (!any(ok)) {
      stop("Could not find p-value for M = ", STATISTIC, ", n = ", n, ", m = ", m)
    }

    PVAL <- Tab$p[ok]
    METHOD <-  "Madhava Rao-Raghunath Exact Symmetry Test"

    RVAL <- list(
      statistic = c(M = STATISTIC),
      parameter = c(n = n, m = m),
      p.value = as.numeric(PVAL),
      method = METHOD,
      data.name = DNAME,
      alternative = "two.sided",
      null.value = c("F(x,y)" = "F(y,x)")
    )
    class(RVAL) <- "htest"
    RVAL
  }

#' @rdname mrrTest
#' @aliases mrrTest.formula
#' @method mrrTest formula
#' @template one-way-formula
#' @importFrom stats model.frame
#' @export
mrrTest.formula <-
  function(formula, data, subset, na.action, ...)
  {
    if (missing(formula)
        || (length(formula) != 3L)
        || (length(attr(terms(formula[-2L]), "term.labels")) != 1L))
      stop("'formula' missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame())))
      m$data <- as.data.frame(data)
    ## need stats:: for non-standard evaluation
    m[[1L]] <- quote(stats::model.frame)
    m$... <- NULL
    mf <- eval(m, parent.frame())
    DNAME <- paste(names(mf), collapse = " by ")
    names(mf) <- NULL
    response <- attr(attr(mf, "terms"), "response")
    g <- factor(mf[[-response]])
    if (nlevels(g) != 2L)
      stop("grouping factor must have exactly 2 levels")
    DATA <- setNames(split(mf[[response]], g),
                     c("x", "y"))
    y <- do.call("mrrTest", c(DATA, list(...)))
    y$data.name <- DNAME
    y
  }
