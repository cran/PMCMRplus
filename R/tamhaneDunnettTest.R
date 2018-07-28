## tamhaneDunnettTest.R
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

#' @name tamhaneDunnettTest
#' @title Tamhane-Dunnett Many-to-One Comparison Test
#' @description
#' Performs Tamhane-Dunnett's multiple comparisons test with one control.
#' For many-to-one comparisons in an one-factorial layout
#' with normally distributed residuals and unequal variances
#' Tamhane-Dunnett's test can be used. A total of \eqn{m = k-1}
#' hypotheses can be tested. The null hypothesis
#' H\eqn{_{i}: \mu_0(x) = \mu_i(x)} is tested in the two-tailed test
#' against the alternative
#' A\eqn{_{i}: \mu_0(x) \ne \mu_i(x), ~~ 1 \le i \le k-1}.
#'
#' The p-values for the test are calculated from the multivariate t distribution
#' as implemented in the function \code{\link[mvtnorm]{pmvt}}.
#'
#' @template class-PMCMR
#' @references
#'  OECD (ed. 2006) \emph{Current approaches in the statistical analysis
#'    of ecotoxicity data: A guidance to application - Annexes}. OECD Series
#'    on testing and assessment, No. 54.
#' @seealso
#' \code{\link[mvtnorm]{pmvt}}, \code{\link{welchManyOneTTest}}
#' @examples
#' set.seed(245)
#' mn <- c(1, 2, 2^2, 2^3, 2^4)
#' x <- rep(mn, each=5) + rnorm(25)
#' g <- factor(rep(1:5, each=5))
#'
#' fit <- aov(x ~ g - 1)
#' shapiro.test(residuals(fit))
#' bartlett.test(x ~ g - 1)
#' anova(fit)
#' summary(tamhaneDunnettTest(x, g, alternative = "greater"))
#'
#' @keywords htest
#' @concept ManyOneComparisons
#' @export
tamhaneDunnettTest <- function(x, ...) UseMethod("tamhaneDunnettTest")

#' @rdname tamhaneDunnettTest
#' @aliases tamhaneDunnettTest.default
#' @method tamhaneDunnettTest default
#' @template one-way-parms
#' @param alternative the alternative hypothesis.
#' Defaults to \code{"two.sided"}.
#' @importFrom mvtnorm pmvt
#' @importFrom stats var
#' @importFrom stats complete.cases
#' @export
tamhaneDunnettTest.default <-
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
            alternative <- x$alternative
            g <- factor(rep(1 : k, l))
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
        x.n <- tapply(!is.na(x), g, length)
        g.unique <- unique(g)
        k <- length(g.unique)
        n <- sum(x.n)

        METHOD <- paste("Tamhane-Dunnett's-test for multiple","
                         comparisons with one control", sep="\t")

        # control is x.mean[1]
        compare.stats <- function(j) {
            numer <- x.mean[j] - x.mean[1]
            denom <- x.var[j] / x.n[j] + x.var[1] / x.n[1]
            STATISTIC <- numer / sqrt(denom)
            return(STATISTIC)
        }
        STATISTIC <- rep(NA, k-1)

        df <- rep(NA, k - 1)
        for (i in 2:k){
            df[i-1] <- x.n[1] + x.n[i] - 2
        }

        for (j in 2:k) {STATISTIC[j-1] <- compare.stats(j)}

        cr <- diag(k-1)
        m <- k-1
        PVAL <- rep(NA, m)




        if (alternative == "greater") {
        # use studentized maximum distribution
        # aka one-sided multivariate t distribution with corr = 0
            for (i in 1:m){
                lo <- -Inf
                up <- rep(STATISTIC[i], m)
                PVAL[i] <- 1 - pmvt(lower = lo,
                                    upper = up,
                                    df = df[i],
                                    corr = cr)
           }
        #    PVAL <- sapply(list(STATISTIC, df),
        #                   function(x) 1 - pmvt(lower=-Inf,
        #                                        upper=rep(x$STATISTIC, m),
        #                                        df = x$df,
        #                                        corr = cr))


        } else if (alternative == "less"){
            for (i in 1:m){
                lo <- rep(STATISTIC[i], m)
                up <- Inf
                PVAL[i] <- 1 - pmvt(lower = lo,
                                    upper = up,
                                    df = df[i],
                                    corr = cr)
            }
          #  PVAL <- sapply(list(STATISTIC, df),
          #                 function(x) 1 - pmvt(lower = rep(x$STATISTIC, m),
          #                                      upper = Inf,
          #                                      df = x$df,
          #                                      corr = cr))
        } else {

           # use Studentized Maximum Modulus Distribution
           # equals Two-sided Multivariate t districution
           # use mvtnorm
           # critical values

            for (i in 1:m){
                lo <- -rep(abs(STATISTIC[i]), m)
                up <- rep(abs(STATISTIC[i]), m)
                PVAL[i] <- 1 - pmvt(lower = lo,
                                    upper = up,
                                    df = df[i],
                                    corr = cr)
            }

         #   PVAL <- sapply(list(STATISTIC, df),
         #                  function(x) 1 - pmvt(lower = -rep(abs(STATISTIC), m),
#                                                upper = rep(abs(STATISTIC), m),
#                                                df = x$df,
#                                                corr = cr))
        }
        LNAME <- levels(g)[2:k]

        # build matrix
        PSTAT <- matrix(data=STATISTIC, nrow = (k-1), ncol = 1,
                        dimnames = list(LNAME, levels(g)[1]))
        PVAL <- matrix(data=PVAL, nrow = (k-1), ncol = 1,
                        dimnames = list(LNAME, levels(g)[1]))
        MODEL <- data.frame(x, g)
        DIST <- "t"
        ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL,
                    statistic = PSTAT, p.adjust.method = "single-step",
                    alternative = alternative, df = df, model =MODEL,
                    dist = DIST)
        class(ans) <- "PMCMR"
        ans
}

#' @rdname tamhaneDunnettTest
#' @aliases tamhaneDunnettTest.formula
#' @method tamhaneDunnettTest formula
#' @template one-way-formula
#' @export
tamhaneDunnettTest.formula <-
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
    y <- do.call("tamhaneDunnettTest", c(as.list(mf), alternative = alternative))
    y$data.name <- DNAME
    y
}
