## chenTest.R
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
#' @name chenTest
#' @title Chen's Many-to-One Comparisons Test
#' @description
#' Performs Chen's nonparametric test for contrasting increasing
#' (decreasing) dose levels of a treatment.
#'
#' @details
#' Chen's test is a non-parametric step-down trend test for
#' testing several treatment levels with a zero control.
#' Let \eqn{X_{0j}} denote a variable with the \eqn{j}-th
#' realization of the control group (\eqn{1 \le j \le n_0})
#' and \eqn{X_{ij}} the \eqn{j}-the realization
#' in the \eqn{i}-th treatment group (\eqn{1 \le i \le k}).
#' The variables are i.i.d. of a least ordinal scale with
#' \eqn{F(x) = F(x_0) = F(x_i), ~ (1 \le i \le k)}.
#' A total of \eqn{m = k} hypotheses can be tested:
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
#' The statistics \eqn{T_i} are based on a Wilcoxon-type ranking:
#'
#' \deqn{
#' T_i = \sum_{j=0}^{i=1} \sum_{u=1}^{n_i} \sum_{v=1}^{n_j} I(x_{iu} - x_{jv}), \qquad (1 \leq i \leq k),
#' }
#'
#' where the indicator function returns \eqn{I(a) = 1, ~ \mathrm{if}~ a > 0, 0.5 ~ \mathrm{if} a = 0}
#' otherwise \eqn{0}.
#'
#' The expected \eqn{i}th mean is
#' \deqn{
#' \mu(T_i) = n_i N_{i-1} / 2,
#' }
#'
#' with \eqn{N_j = \sum_{j =0}^i n_j} and the \eqn{i}th variance:
#'
#' \deqn{
#' \sigma^2(T_i) = n_i N_{i-1} / 12 ~ \left\{N_i + 1 -
#' \sum_{j=1}^g t_j \left(t_j^2 - 1 \right) /
#' \left[N_i \left( N_i - 1 \right)\right]\right\}.
#' }
#'
#' The test statistic \eqn{T_i^*} is asymptotically standard normal
#'
#' \deqn{
#'  T_i^* = \frac{T_i - \mu(T_i)}
#'  {\sqrt{\sigma^2(T_i)}}, \qquad (1 \leq i \leq k).
#' }
#'
#' The p-values are calculated from the standard normal distribution.
#' The p-values can be adjusted with any method as available
#' by \code{\link{p.adjust}} or by the step-down procedure as proposed
#' by Chen (1999), if \code{p.adjust.method = "SD1"}.
#'
#' @inherit cuzickTest note
#' @template class-PMCMR
#'
#' @examples
#' ## Chen, 1999, p. 1237,
#' ## Minimum effective dose (MED)
#' ## is at 2nd dose level
#' df <- data.frame(x = c(23, 22, 14,
#' 27, 23, 21,
#' 28, 37, 35,
#' 41, 37, 43,
#' 28, 21, 30,
#' 16, 19, 13),
#' g = gl(6, 3))
#' levels(df$g) <- 0:5
#' ans <- chenTest(x ~ g, data = df, alternative = "greater",
#'                 p.adjust.method = "SD1")
#' summary(ans)
#'
#' @references
#' Chen, Y.-I., 1999, Nonparametric Identification of the
#' Minimum Effective Dose. *Biometrics* **55**, 1236--1240.
#' \doi{10.1111/j.0006-341X.1999.01236.x}
#'
#' @seealso
#' \code{\link{wilcox.test}}, \code{\link{Normal}}
#' @concept wilcoxonranks
#' @keywords htest nonparametric
#' @export
chenTest <- function(x, ...)
  UseMethod("chenTest")

#' @rdname chenTest
#' @method chenTest default
#' @aliases chenTest.default
#' @template one-way-parms
#' @param alternative the alternative hypothesis. Defaults to \code{two.sided}.
#' @param p.adjust.method method for adjusting p values
#' (see \code{\link{p.adjust}})
#' @importFrom stats pnorm
#' @importFrom stats p.adjust
#' @importFrom stats p.adjust.methods
#' @importFrom stats complete.cases
#' @export
chenTest.default <-
  function(x,
           g,
           alternative = c("greater", "less"),
           p.adjust.method =  c("SD1", p.adjust.methods),
           ...) {
    ## taken from stats::kruskal.test

    if (is.list(x)) {
      stop("'x' must be a list with at least 2 elements")
      DNAME <- deparse(substitute(x))
      x <- lapply(x, function(u)
        u <- u[complete.cases(u)])
      k <- length(x)
      l <- sapply(x, "length")
      if (any(l == 0))
        stop("all groups must contain data")
      g <- factor(rep(1:k, l))
      if (is.null(x$alternative)) {
        alternative <- "two.sided"
      } else {
        alternative <- x$alternative
      }
      if (is.null(x$p.adjust.method)) {
        p.adjust.method <- "single-step"
      } else {
        p.adjust.method <- x$p.adjust.method
      }
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

    # check arguments
    alternative <- match.arg(alternative)
    p.adjust.method <- match.arg(p.adjust.method)


    if (alternative == "less") {
      x <- -x
    }

    k <- nlevels(g)
    n <- tapply(x, g, length)
    glev <- levels(g)

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

    # Chen, 1999, p.1237
    m <- k - 1
    T <- rep(0, m)
    for (i in 1:m) {
      for (j in 0:(i - 1)) {
        Yi <- x[g == glev[i + 1]]
        Yj <- x[g == glev[j + 1]]
        for (u in 1:n[i + 1]) {
          for (v in 1:n[j + 1]) {
            T[i] <- T[i] +
              sum(I(Yi[u] - Yj[v]))
          }
        }
      }
    }

    ## mean
    N <- cumsum(n)
    mu <- rep(0, m)
    for (i in 1:m) {
      mu[i] <- n[i + 1] * N[i] / 2
    }

    ## variance
    sigma2 <- rep(0, m)
    for (i in 1:m) {
      tmp <- x[g %in% glev[1:(i + 1)]]
      C <- gettiesFW(tmp)
      sigma2[i] <- n[i + 1] * N[i] *
        (N[i + 1] + 1 - C / (N[i + 1] * (N[i + 1] - 1))) /
        12
    }

    ##
    STATISTIC <- (T - mu) / sqrt(sigma2)

    ## p-values, one-sided
    pval <- pnorm(STATISTIC , lower.tail = FALSE)

    ## p-adjustment
    if (p.adjust.method == "SD1") {
      padj <- SD1p(pval)
    } else {
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
    METHOD <- paste("Chen's test\n",
                    "\tfor multiple comparisons with one control",
                    sep = "")
    MODEL <- data.frame(x, g)

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

#' @rdname chenTest
#' @method chenTest formula
#' @aliases chenTest.formula
#' @template one-way-formula
#' @export
chenTest.formula <-
  function(formula,
           data,
           subset,
           na.action,
           alternative = c("greater", "less"),
           p.adjust.method = c("SD1", p.adjust.methods),
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
    p.adjust.method <- match.arg(p.adjust.method)
    names(mf) <- NULL
    y <-
      do.call(
        "chenTest",
        c(
          as.list(mf),
          alternative = alternative,
          p.adjust.method = p.adjust.method
        )
      )
    y$data.name <- DNAME
    y
  }
