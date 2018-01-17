# shirleyWilliamsTest.R
# Part of the R package: PMCMR
#
##  Copyright (C) 2017 Thorsten Pohlert
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
#' The procedure starts from the highest dose level (\eqn{m}) to the the lowest dose level (\eqn{1}) and
#' stops at the first non-significant test. The consequent lowest effect dose
#' is the treatment level of the previous test number.
#'
#' The p-values are estimated through an assymptotic boot-strap method. The p-values for H\eqn{_1}
#' are calculated from the t distribution with infinite degree of freedom. This function has
#' included the modifications as recommended by Williams (1986).
#'
#' @note
#' One may increase the number of permutations to e.g. \code{nperm = 10000}
#' in order to get more precise p-values. However, this will be on the expense of
#' computational time.
#' 
#' @references
#' Shirley, E., (1977), Nonparametric Equivalent of Williams Test for Contrasting Increasing
#' Dose Levels of a Treatment. \emph{Biometrics}, 33, 386--389.
#' 
#' Williams, D.A. (1986), Note on Shirley's nonparametric test for comparing
#' several dose levels with a zero-dose control. \emph{Biometrics} 42, 183--186.
#' 
#' @template class-PMCMR
#' @seealso
#' \code{\link{sample}}
#'
#' @keywords htest nonparametric
#' @concept StepDownTrendTest
#' @concept DoseFinding
#' 
#' @examples
#' ## Example from Sachs (1997, p. 402)
#' x <- c(106, 114, 116, 127, 145,
#' 110, 125, 143, 148, 151,
#' 136, 139, 149, 160, 174)
#' g <- gl(3,5)
#' levels(g) <- c("0", "I", "II")
#'
#' ## Shirley-Williams Test
#' shirleyWilliamsTest(x ~ g)
#' 
#' @export
shirleyWilliamsTest <- function(x, ...) UseMethod("shirleyWilliamsTest")

#' @rdname shirleyWilliamsTest
#' @aliases shirleyWilliamsTest.default
#' @method shirleyWilliamsTest default
#' @template one-way-parms
#' @param nperm number of permutations for the assymptotic permutation test.
#' Defaults to \code{1000}.
#' @importFrom stats pt
#' @export
shirleyWilliamsTest.default <-
function(x, g, nperm=1000, ...){
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

    nj <- tapply(x, g, length)
    k <- nlevels(g)

    ## comparisons
    compfn <- function(x, ix, g, nj){
        k <- length(nj)
        ti <- rep(NA, k)
        x <- x[ix]
        for (i in 2:k){
            N <- sum(nj[1:i])
            r <- rank(x[1:N])
            gg <- g[1:N]
            Rj <- tapply(r, gg, mean)
            t <- table(r)
            names(t) <- NULL
            T <- sum((t^3 - t) / (12 * (N - 1)))
            Vi <- N * (N + 1) / 12 - T
            u <- 2:i
            j <- u
            enum <- sapply(j, function(j) sum(nj[j:i] * Rj[j:i]))
            denom <- sapply(j, function(j) sum(nj[j:i]))
        
            ti[i] <- (max(enum / denom) - Rj[1]) /
                sqrt(Vi * (1/nj[i] + 1/nj[1]))
            }
        return(ti[2:k])
    }

    N <- sum(nj)
    l <- 1:N
    STATISTIC <- compfn(x, l, g, nj)

    ## permutation
    mt <- matrix(NA, ncol=(k-1), nrow=nperm)
    for(i in 1:nperm){
        ix <- sample(l)
        mt[i,] <- compfn(x, ix, g, nj)
    }

    ## pvalues
    PVAL <- sapply(1:(k-1), function(j) {
        p <- (sum(mt[,j] <= -abs(STATISTIC[j])) +
              (sum(mt[,j] >= abs(STATISTIC[j]))))/nperm
        p}
        )

    ## exact p-value from student t distribution
    PVAL[1] <- 2 * min(0.5, pt(STATISTIC[1], df=Inf, lower.tail=FALSE))
     
    STAT <- cbind(ctr = STATISTIC)
    row.names(STAT) <- sapply(1:(k-1), function(i) paste0("mu",i))
    P <- cbind(ctr = PVAL)
    row.names(P) <- row.names(STAT)
    
    DAT <- data.frame(x, g)
    METH <- c("Shirley-Williams test")
    ans <- list(statistic = STAT, 
                p.value = P,
                data = DAT,
                method = METH,
                data.name = DNAME,
                alternative = "two.sided",
                dist = "t",
                p.adjust.method = "boot")
    class(ans) <- "PMCMR"
    ans
}

#' @rdname shirleyWilliamsTest
#' @aliases shirleyWilliamsTest.formula
#' @method shirleyWilliamsTest formula
#' @template one-way-formula
#' @export
shirleyWilliamsTest.formula <-
    function(formula, data, subset, na.action, nperm=1000,...)
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
    names(mf) <- NULL
    y <- do.call("shirleyWilliamsTest", c(as.list(mf)))
    y$data.name <- DNAME
    y
}
