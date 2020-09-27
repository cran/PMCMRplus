## siegelTukeyTest.R
## Part of the R package: PMCMR
##
## Copyright (C) 2018 Thorsten Pohlert
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
#' @name siegelTukeyTest
#' @title Siegel-Tukey Rank Dispersion Test
#' @description Performs Siegel-Tukey non-parametric
#' rank dispersion test.
#' @details
#' Let \eqn{x} and \eqn{y} denote two identically and independently
#' distributed variables of at least ordinal scale.
#'  Further, let
#'  \eqn{\theta}{theta}, and \eqn{\lambda}{lambda} denote
#'  location and scale parameter of the common, but unknown distribution.
#'  Then for the two-tailed case, the null hypothesis
#'  H: \eqn{\lambda_x / \lambda_y = 1 | \theta_x = \theta_y}{%
#'  lambdaX / lambdaY = 1 | thetaX = thetaY} is
#'  tested against the alternative,
#'  A: \eqn{\lambda_x / \lambda_y \ne 1}{lambdaX / lambdaY != 1}.
#'
#'  The data are combinedly ranked according to Siegel-Tukey.
#'  The ranking is done by alternate extremes (rank 1 is lowest,
#'  2 and 3 are the two highest, 4 and 5 are the two next lowest, etc.).
#'  If no ties are present, the p-values are computed from
#'  the Wilcoxon distribution (see \code{\link{Wilcoxon}}).
#'  In the case of ties, a tie correction is done according
#'  to Sachs (1997) and approximate p-values are computed
#'  from the standard normal distribution (see \code{\link{Normal}}).
#'
#'  If both medians differ, one can correct for medians to
#' increase the specificity of the test.
#'
#' @source
#' The algorithm for the Siegel-Tukey ranks was
#' taken from the code of Daniel Malter. See also the
#' blog from Tal Galili (02/2010, \url{https://www.r-statistics.com/2010/02/siegel-tukey-a-non-parametric-test-for-equality-in-variability-r-code/},
#' accessed 2018-08-05).
#'
#' @references
#' Sachs, L. (1997), \emph{Angewandte Statistik}. Berlin: Springer.
#'
#' Siegel, S., Tukey, J. W. (1960), A nonparametric sum of ranks
#' procedure for relative spread in unpaired samples,
#' \emph{Journal of the American Statistical Association} \bold{55}, 429--455.
#' @examples
#' ## Sachs, 1997, p. 376
#' A <- c(10.1, 7.3, 12.6, 2.4, 6.1, 8.5, 8.8, 9.4, 10.1, 9.8)
#' B <- c(15.3, 3.6, 16.5, 2.9, 3.3, 4.2, 4.9, 7.3, 11.7, 13.7)
#' siegelTukeyTest(A, B)
#'
#' ## from example var.test
#' x <- rnorm(50, mean = 0, sd = 2)
#' y <- rnorm(30, mean = 1, sd = 1)
#' siegelTukeyTest(x, y, median.corr = TRUE)
#'
#' ## directional hypothesis
#' A <- c(33, 62, 84, 85, 88, 93, 97)
#' B <- c(4, 16, 48, 51, 66, 98)
#' siegelTukeyTest(A, B, alternative = "greater")
#'
#' @keywords htest
#' @keywords nonparametric
#' @template class-htest
#' @export
siegelTukeyTest <- function(x, ...) UseMethod("siegelTukeyTest")

#' @rdname siegelTukeyTest
#' @aliases siegelTukeyTest.default
#' @method siegelTukeyTest default
#' @param x,y numeric vectors of data values.
#' @param alternative a character string specifying the
#' alternative hypothesis, must be one of \code{"two.sided"} (default),
#' \code{"greater"} or \code{"less"}.
#' You can specify just the initial letter.
#' @param median.corr logical indicator, whether median correction
#' should be performed prior testing. Defaults to \code{FALSE}.
#' @importFrom stats pwilcox pnorm median
#' @export
siegelTukeyTest.default <- function(x, y, alternative = c("two.sided", "greater", "less"),
                                    median.corr = FALSE, ...) {

    alternative <- match.arg(alternative)
    DNAME <- paste0(deparse(substitute(x)), " and ",
                    deparse(substitute(y)))

    ## correct for median differences
    if (median.corr) {
        med.x <- median(x)
        med.y <- median(y)
        x <- x - (med.x - med.y) / 2
        y <- y - (med.y - med.x) / 2
    }
    xy <- c(x, y)

    g <- as.factor(c(rep(1, length(x)), rep(2, length(y))))
    ord <- order(xy)
    gord <- g[ord]
    xyord <- xy[ord]
    N <- length(xy)

    ## rank sequence
    ## based on code from Daniel Malter
    a <- rep(seq(ceiling(N / 4)), each=2)
    b <- rep(c(0, 1), ceiling(N)/4)
    suppressWarnings(
        rk.up <- c(1, (a * 4 + b))[1:ceiling(N / 2)]
    )
    suppressWarnings(
        rk.down <- rev(c(a * 4 + b - 2)[1:floor(N / 2)])
    )

    r <- c(rk.up, rk.down)

    ## rank sums and sample sizes
    R <- tapply(r, gord, sum)
    n <- tapply(r, gord, length)
    names(R) <- c("R1", "R2")
    names(n) <- NULL

    ## check for ties
    ties <- sum(table(xy) - 1) > 0

    if (!ties) {
        ## pWilcox
        ## valid for samples without ties
        W <- R[1]  - n[1] * (n[1] + 1) / 2

        p.value <- suppressWarnings(
            switch(alternative,
                   "two.sided" = {
                       if (W > (n[1] * n[2]) / 2) {
                           p <- 2 * min(0.5, pwilcox(W - 1, n[1], n[2], lower.tail = FALSE))
                       } else {
                           p <- 2 * min(0.5, pwilcox(W, n[1], n[2]))
                       }
                       p
                   },
                   "less" = pwilcox(W-1, n[1], n[2], lower.tail = FALSE),

                   "greater" = pwilcox(W, n[1], n[2])
                   )
        )
        statistic <- W
        names(statistic) <- "W"
        parameters <- NULL
        ##names(parameters) <- c("m","n")

    } else {

        warning("Ties are present. Compute approximate p-values.")

        ## approximate for ties
        uniqVals <- unique(xyord)
        rcor <- tapply(r, xyord, mean)
        corected <- data.frame(rcor, uniqVals)
        uncor <- data.frame(xyord, gord, r)
        dat <- merge(x = uncor, y = corected,
                     by.x = "xyord",
                     by.y = "uniqVals")

        tmp <- as.numeric(names(which(table(dat$xyord) > 1)))
        f <- subset(dat, subset = xyord %in% tmp)

        S1 <- sum(f$r^2)
        S2 <- sum(f$rcor^2)

        ## statistic (Sachs, 1997, p. 375)
        a <- ifelse(2 * R[1] > n[1] * (n[1] + n[2] + 1), -1, 1)
        enum <- 2 * R[1] - n[1] * (n[1] + n[2] + 1) + a
        denom <- sqrt(n[1] * (n[1] + n[2] + 1) * (n[2] / 3) -
                      4 * (n[1] * n[2] /
                           ((n[1] + n[2]) * (n[1] + n[2] - 1))) *
                      (S1 - S2))
        statistic <- enum / denom
        names(statistic) <- "z"
        parameters <- NULL

        ## two-sided test
        p.value <- switch(alternative,

                          "two.sided" = 2 * min(0.5, pnorm(abs(statistic), lower.tail = FALSE)),

                          "less" = pnorm(statistic, lower.tail = FALSE),

                          "greater" = pnorm(statistic)
                          )
    }

    method <- "Siegel-Tukey rank dispersion test"

    ans <- list(p.value = p.value,
                statistic = statistic,
                method = method,
                alternative = alternative,
                parameters = parameters,
                null.value = c("ratio of scales" = 1),
                data.name = DNAME)
    class(ans) <- "htest"
    ans
}

#' @rdname siegelTukeyTest
#' @method siegelTukeyTest formula
#' @aliases siegelTukeyTest.formula
#' @template one-way-formula
#' @importFrom stats setNames terms
#' @export
siegelTukeyTest.formula <-
function(formula, data, subset, na.action, ...)
{
    if(missing(formula)
       || (length(formula) != 3L)
       || (length(attr(terms(formula[-2L]), "term.labels")) != 1L))
        stop("'formula' missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if(is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)
    ## need stats:: for non-standard evaluation
    m[[1L]] <- quote(stats::model.frame)
    m$... <- NULL
    mf <- eval(m, parent.frame())
    DNAME <- paste(names(mf), collapse = " by ")
    names(mf) <- NULL
    response <- attr(attr(mf, "terms"), "response")
    g <- factor(mf[[-response]])
    if(nlevels(g) != 2L)
        stop("grouping factor must have exactly 2 levels")
    DATA <- setNames(split(mf[[response]], g), c("x", "y"))
    y <- do.call("siegelTukeyTest", c(DATA, list(...)))
    y$data.name <- DNAME
    y
}
