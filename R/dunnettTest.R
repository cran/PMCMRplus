## dunnettTest.R
## Part of the R package: PMCMR
##
## Copyright (C) 2017-2020 Thorsten Pohlert
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
##
## Requires package mvtnorm
#' @name dunnettTest
#' @title Dunnett's Many-to-One Comparisons Test
#'
#' @description
#' Performs Dunnett's multiple comparisons test with one control.
#' @details
#' For many-to-one comparisons in an one-factorial layout
#' with normally distributed residuals Dunnett's test
#' can be used.
#' Let \eqn{X_{0j}} denote a continuous random variable
#' with the \eqn{j}-the realization of the control group
#' (\eqn{1 \le j \le n_0}) and \eqn{X_{ij}} the \eqn{j}-the realization
#' in the \eqn{i}-th treatment group (\eqn{1 \le i \le k}).
#' Furthermore, the total sample size is \eqn{N = n_0 + \sum_{i=1}^k n_i}.
#' A total of \eqn{m = k} hypotheses can be tested: The null hypothesis is
#' H\eqn{_{i}: \mu_i = \mu_0} is tested against the alternative
#' A\eqn{_{i}: \mu_i \ne \mu_0} (two-tailed). Dunnett's test
#' statistics are given by
#'
#' \deqn{
#'  t_{i} \frac{\bar{X}_i - \bar{X_0}}
#'  {s_{\mathrm{in}} \left(1/n_0 + 1/n_i\right)^{1/2}}, ~~
#'  (1 \le i \le k)
#' }{%
#'  SEE PDF
#' }
#'
#' with \eqn{s^2_{\mathrm{in}}} the within-group ANOVA variance.
#' The null hypothesis is rejected if
#' \eqn{|t_{ij}| > |T_{kv\rho\alpha}|} (two-tailed),
#' with \eqn{v = N - k} degree of freedom and \eqn{rho} the correlation:
#'
#' \deqn{
#'  \rho_{ij} = \sqrt{\frac{n_i n_j}
#'              {\left(n_i + n_0\right) \left(n_j+ n_0\right)}} ~~
#'              (i \ne j)
#' .}{%
#'  SEE PDF
#' }
#'
#' The p-values are computed with the function \code{\link{pDunnett}}
#' that is a wrapper to the the multivariate-t distribution as implemented in the function
#' \code{\link[mvtnorm]{pmvt}}.
#'
#' @template class-PMCMR
#'
#' @references
#' Dunnett, C. W. (1955) A multiple comparison procedure for comparing several
#'   treatments with a control. \emph{Journal of the American Statistical Association}
#'   \bold{50}, 1096–1121.
#'
#'  OECD (ed. 2006) \emph{Current approaches in the statistical analysis
#'    of ecotoxicity data: A guidance to application - Annexes}. OECD Series
#'    on testing and assessment, No. 54.
#' @seealso
#' \code{\link[mvtnorm]{pmvt}} \code{\link{pDunnett}}
#' @examples
#' fit <- aov(Y ~ DOSE, data = trout)
#' shapiro.test(residuals(fit))
#' bartlett.test(Y ~ DOSE, data = trout)
#'
#' ## works with fitted object of class aov
#' summary(dunnettTest(fit, alternative = "less"))
#'
#' @keywords htest
#' @concept parametric
#' @importFrom stats var
#' @importFrom stats complete.cases
#' @export
dunnettTest <- function(x, ...) UseMethod("dunnettTest")

#' @rdname dunnettTest
#' @aliases dunnettTest.default
#' @method dunnettTest default
#' @template one-way-parms-aov
#' @param alternative the alternative hypothesis. Defaults to \code{two.sided}.
#' @export
dunnettTest.default <-
function(x, g, alternative = c("two.sided", "greater", "less"), ...){
        ## taken from stats::kruskal.test
    if (is.list(x)) {
        if (length(x) < 2L)
            stop("'x' must be a list with at least 2 elements")
        DNAME <- deparse(substitute(x))
        x <- lapply(x, function(u) u <- u[complete.cases(u)])
        k <- length(x)
        l <- sapply(x, "length")
        if (any(l == 0))
            stop("all groups must contain data")
        g <- factor(rep(1 : k, l))
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
        if (k < 2)
            stop("all observations are in the same group")
    }
    alternative <- match.arg(alternative)

        # Parametric
    x.mean <- tapply(x, g, mean, na.rm = TRUE)
    x.var <- tapply(x, g, var, na.rm = TRUE)
    n <- tapply(!is.na(x), g, length)

    x.var.pool <- sum((n - 1) * x.var) / sum(n - 1)
    x.sd.pool <- sqrt(x.var.pool)

    g.unique <- unique(g)
    k <- length(g.unique)

    METHOD <- paste("Dunnett's-test for multiple","
                         comparisons with one control", sep="\t")

        # control is x.mean[1]
    compare.stats <- function(j) {
        numer <- x.mean[j] - x.mean[1]
        denom <- x.sd.pool * sqrt(1 / n[j] + 1 / n[1])
        STATISTIC <- numer / denom
        return(STATISTIC)
    }

	# correlation matrix
    n0 <- n[1]
    nn <- n[2:k]
    kk <- k - 1

    df <- length(x) - k
    STATISTIC <- rep(NA, k - 1)

	# Get statistic values
    for (j in 2:k) {
        STATISTIC[j-1] <- compare.stats(j)
    }

  # get p-values from novel pDunnett
    if (alternative == "two.sided") {
      PVAL <- pmin(1, 2 * pDunnett(q = abs(STATISTIC),
                               n0 = n0,
                               n = nn,
                               lower.tail = FALSE))
    } else if (alternative == "greater") {
      PVAL <- pDunnett(q = STATISTIC,
                                   n0 = n0,
                                   n = nn,
                                   lower.tail = FALSE)
    } else {
      PVAL <- pDunnett(q = STATISTIC,
                       n0 = n0,
                       n = nn,
                       lower.tail = TRUE)

    }


        # Names
    LNAME <- levels(g)[2:k]
    PARMS <- c(kk, df)
    names(PARMS) <- c("k", "df")

        # build matrix
    PSTAT <- matrix(data=STATISTIC, nrow = (k-1), ncol = 1,
                    dimnames = list(LNAME, levels(g)[1]))
    PVAL <- matrix(data=PVAL, nrow = (k-1), ncol = 1,
                   dimnames = list(LNAME, levels(g)[1]))
    MODEL <- data.frame(x, g)
    DIST <- "t"
    ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL,
                statistic = PSTAT, p.adjust.method = "single-step",
                alternative = alternative, parameter = PARMS,
                model = MODEL, dist = DIST)
    class(ans) <- "PMCMR"
    ans
}

#' @rdname dunnettTest
#' @aliases dunnettTest.formula
#' @method dunnettTest formula
#' @template one-way-formula
#' @export
dunnettTest.formula <-
function(formula, data, subset, na.action, alternative = c("two.sided", "greater", "less"), ...)
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
    y <- do.call("dunnettTest", c(as.list(mf), alternative = alternative))
    y$data.name <- DNAME
    y
}

##
#' @rdname dunnettTest
#' @aliases dunnettTest.aov
#' @method dunnettTest aov
# @param obj A fitted model object, usually an \link[stats]{aov} fit.
#' @export
dunnettTest.aov <- function(x, alternative = c("two.sided", "greater", "less"), ...) {
    model <- x$model
    DNAME <- paste(names(model), collapse = " by ")
    names(model) <- c("x", "g")
    alternative <- match.arg(alternative)
    parms <- c(as.list(model), list(alternative = alternative))
    y <- do.call("dunnettTest", parms)
    y$data.name <- DNAME
    y
}


