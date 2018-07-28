## welchManyOneTTest.R
## Part of the R package: PMCMRplus
##
## Copyright (C) 2017, 2018 Thorsten Pohlert
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

#' @name welchManyOneTTest
#' @title Welchs's Many-To-One Comparison Test
#' @template class-PMCMR
#'
#' @description
#' Performs Welchs's t-test for multiple comparisons with one control.
#'
#' @details
#' For many-to-one comparisons in an one-factorial layout
#' with normally distributed residuals and unequal variances
#' Welch's t-test can be used. A total of \eqn{m = k-1}
#' hypotheses can be tested. The null hypothesis
#' H\eqn{_{i}: \mu_0(x) = \mu_i(x)} is tested in the two-tailed test
#' against the alternative
#' A\eqn{_{i}: \mu_0(x) \ne \mu_i(x), ~~ 1 \le i \le k-1}.
#'
#' This function is basically a wrapper function for
#' \code{\link[stats]{t.test}(..., var.equal = FALSE)}. The p-values for the test
#' are calculated from the t distribution
#' and can be adusted with any method that is implemented in
#' \code{\link[stats]{p.adjust.methods}}.
#'
#' @references
#'  Welch, B. L. (1947) The generalization of "Student's" problem
#'  when several different population variances are involved,
#'  \emph{Biometrika} \bold{34}, 28--35.
#'
#'  Welch, B. L. (1951) On the comparison of several mean values:
#'  An alternative approach, \emph{Biometrika} \bold{38}, 330--336.
#'
#' @keywords htest
#' @concept ManyToOneComparisons
#' @examples
#' set.seed(245)
#' mn <- rep(c(1, 2^(1:4)), each=5)
#' sd <- rep(1:5, each=5)
#' x <- mn + rnorm(25, sd = sd)
#' g <- factor(rep(1:5, each=5))
#'
#' fit <- aov(x ~ g)
#' shapiro.test(residuals(fit))
#' bartlett.test(x ~ g) # var1 != varN
#' anova(fit)
#' summary(welchManyOneTTest(x, g, alternative = "greater", p.adjust="holm"))
#'
#' @seealso
#' \code{\link[stats]{pairwise.t.test}}, \code{\link[stats]{t.test}},
#' \code{\link[stats]{p.adjust}}, \code{\link{tamhaneDunnettTest}}
#'
#' @importFrom stats t.test
#' @importFrom stats p.adjust.methods p.adjust
#' @export
welchManyOneTTest <- function(x, ...)
  UseMethod("welchManyOneTTest")

#' @rdname welchManyOneTTest
#' @method welchManyOneTTest default
#' @aliases welchManyOneTTest.default
#' @template one-way-parms
#' @param alternative the alternative hypothesis.
#' Defaults to \code{two.sided}.
#' @param p.adjust.method  method for adjusting p values
#'    (see \code{\link{p.adjust}}).
#' @export
welchManyOneTTest.default <-
  function(x,
           g,
           alternative = c("two.sided", "greater", "less"),
           p.adjust.method = p.adjust.methods,
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
      ##
      if (is.null(x$alternative)){
        alternative <- "two.sided"
      } else {
        alternative <- x$alternative
      }
      if(is.null(x$p.adjust.method)){
        p.adjust.method <- "holm"
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

    alternative <- match.arg(alternative)
    p.adjust.method <- match.arg(p.adjust.method)

    ## prepare factors
    kk <- k - 1
    levNames <- levels(g)

    ## prepare output
    statistic <- numeric(kk)
    p.value <- numeric(kk)

    ## Control is x0
    x0 <- x[g == levNames[1]]

    for (i in 1:kk) {
      out <- t.test(
        y = x0,
        x = x[g == levNames[i + 1]],
        alternative = alternative,
        var.equal = FALSE
      )

      statistic[i] <- out$statistic
      p.value[i] <- out$p.value
    }

    p.value <- p.adjust(p.value,
                        method = p.adjust.method)

    METHOD <- "Welch's t-test"

    STAT <- cbind(statistic)
    colnames(STAT) <- levNames[1]
    rownames(STAT) <- levNames[2:k]
    PVAL <- cbind(p.value)
    colnames(PVAL) <- colnames(STAT)
    rownames(PVAL) <- rownames(STAT)

    MODEL <- data.frame(x, g)
    DIST <- "t"
    ans <- list(
      method = METHOD,
      data.name = DNAME,
      p.value = PVAL,
      statistic = STAT,
      p.adjust.method = p.adjust.method,
      model = MODEL,
      dist = DIST,
      alternative = alternative
    )
    class(ans) <- "PMCMR"
    ans
  }

#' @rdname welchManyOneTTest
#' @method welchManyOneTTest formula
#' @aliases welchManyOneTTest.formula
#' @template one-way-formula
#' @export
welchManyOneTTest.formula <-
  function(formula, data, subset, na.action,
           alternative = c("two.sided", "greater", "less"),
           p.adjust.method = p.adjust.methods,...)
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
    names(mf) <- NULL
    alternative <- match.arg(alternative)
    p.adjust.method <- match.arg(p.adjust.method)
    y <- do.call("welchManyOneTTest", c(as.list(mf),
                                        alternative = alternative,
                                        p.adjust.method = p.adjust.method))
    y$data.name <- DNAME
    y
  }
