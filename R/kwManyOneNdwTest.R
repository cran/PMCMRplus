## kwManyOneNdwTest.R
## Part of the R package: PMCMR
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

##
## basically the same as dunn.test.control but no tie correction

#' @rdname kwManyOneNdwTest
#' @title Nemenyi-Damico-Wolfe Many-to-One Rank Comparison Test
#' @description
#' Performs Nemenyi-Damico-Wolfe non-parametric many-to-one comparison
#' test for Kruskal-type ranked data.
#'
#' @details
#' For many-to-one comparisons (pairwise comparisons with one control)
#' in an one-factorial layout with non-normally distributed
#' residuals the Nemenyi-Damico-Wolfe non-parametric test can be performed.
#' Let there be \eqn{k} groups including the control,
#' then the number of treatment levels is \eqn{m = k - 1}.
#' Then \eqn{m} pairwise comparisons can be performed between
#' the \eqn{i}-th treatment level and the control.
#' H\eqn{_i: \theta_0 = \theta_i} is tested in the two-tailed case against
#' A\eqn{_i: \theta_0 \ne \theta_i, ~~ (1 \le i \le m)}.
#'
#' If \code{p.adjust.method == "single-step"} is selected,
#' the \eqn{p}-values will be computed
#' from the multivariate normal distribution. Otherwise,
#' the \eqn{p}-values are computed from the standard normal distribution using
#' any of the \eqn{p}-adjustment methods as included in \code{\link{p.adjust}}.
#'
#' @note
#' This function is essentially the same as \code{\link{kwManyOneDunnTest}}, but
#' there is no tie correction included. Therefore, the implementation of
#' Dunn's test is superior, when ties are present.
#'
#' @references
#' Damico, J. A., Wolfe, D. A. (1989) Extended tables of the exact distribution of
#' a rank statistic for treatments versus control multiple comparisons in one-way
#' layout designs, \emph{Communications in Statistics - Theory and Methods} \bold{18},
#' 3327--3353.
#'
#' Nemenyi, P. (1963) \emph{Distribution-free Multiple Comparisons},
#' Ph.D. thesis, Princeton University.
#'
#' @template class-PMCMR
#' @concept ManyToOneComparisons
#' @concept Kruskal
#' @concept Rank
#' @keywords htest nonparametric
#' @seealso
#' \code{\link{pmvt}}, \code{\link{TDist}}, \code{\link{kruskalTest}},
#' \code{\link{kwManyOneDunnTest}}, \code{\link{kwManyOneConoverTest}}
#' @example examples/kwManyOneMC.R
#' @export
kwManyOneNdwTest <- function(x, ...) UseMethod("kwManyOneNdwTest")

#' @rdname kwManyOneNdwTest
#' @method kwManyOneNdwTest default
#' @aliases kwManyOneNdwTest.default
#' @template one-way-parms
#' @param alternative the alternative hypothesis. Defaults to \code{two.sided}.
#' @param p.adjust.method  method for adjusting p values
#'    (see \code{\link{p.adjust}}).
#' @importFrom stats pnorm
#' @importFrom mvtnorm pmvnorm
#' @export
kwManyOneNdwTest.default <-
    function(x, g, alternative = c("two.sided", "greater", "less"),
             p.adjust.method = c("single-step", p.adjust.methods), ...){
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
            ##
            if (is.null(x$alternative)){
                alternative <- "two.sided"
            } else {
                alternative <- x$alternative
            }
            if(is.null(x$p.adjust.method)){
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

        # Check arguments
        p.adjust.method <- match.arg(p.adjust.method)
        alternative <- match.arg(alternative)

        # Preparation
        x.rank <- rank(x)
        R.bar <- tapply(x.rank, g, mean, na.rm=T)
        R.n <- tapply(!is.na(x), g, length)
        k <- nlevels(g)
        n <- sum(R.n)

	# tie function
        getties <- function(x){
            n <- length(x)
            t <- table(x)
            C <- 1 - (sum(t^3 - t)) / (n^3 - n)
            C <- min(c(1,C))
            return(C)
        }

        C <- getties(x)
        if(C < 1){
            warning("Ties are present. p values are not corrected.")
        }

        compare.stats <- function(i) {
		# Control is in first element
            dif <- abs(R.bar[i] - R.bar[1])
            qval <- dif / sqrt((n * (n + 1) / 12) * (1/R.n[i] + 1/R.n[1] ))
            return(qval)
        }

        pstat <- as.vector(sapply(2:k, function(i) compare.stats(i)))

        if (p.adjust.method != "single-step") {
            if (alternative == "two.sided")
            {
                pval <- 2 * pnorm(abs(pstat), lower.tail = FALSE)
            } else if (alternative == "greater"){
                pval <- pnorm(pstat,
                              lower.tail = FALSE)
            } else {
                pval <- pnorm(pstat)
            }
            pvalv <- p.adjust(pval, method=p.adjust.method)
        } else {
            ## use function pmvnorm of package mvtnorm
            m <- k - 1
            df <- n - k
            # correlation matrix
            ni <- tapply(x, g, length)
            n0 <- ni[1]
            nn <- ni[2:k]
            cr <- diag(m)
            # valid also for unequal sample sizes
            for (i in 1:(m-1)){
                for (j in (i+1):m){
                    cr[i,j] <- ((nn[i] * nn[j]) /
                                ((nn[i] + n0) * (nn[j] + n0)))^(1/2)
                    cr[j,i] <- cr[i, j]
                    cr[j,i] <- cr[i, j]
                }
            }
            if (alternative == "two.sided"){
                pvalv <- sapply(pstat, function(x)
                    1 - pmvnorm(lower = -rep(abs(x), m),
                                upper = rep(abs(x), m),
                                corr = cr))
            } else if (alternative == "greater"){
                pvalv <- sapply(pstat, function(x)
                    1 - pmvnorm(lower = -Inf,
                                upper = rep(x, m),
                                corr = cr))
            } else {
                pvalv <- sapply(pstat, function(x)
                    1 - pmvnorm(lower = rep(x, m),
                                upper = Inf,
                                corr = cr))
            }
        }

        METHOD <- "Nemenyi-Damico-Wolfe many-to-one test"
        #grpn <- levels(g)
        LNAME <- levels(g)[2:k]
        ## PVAL <- cbind(pvalv)
        ## PSTAT <- cbind(pstat)
        PSTAT <- matrix(data=pstat, nrow = (k-1), ncol = 1,
                    dimnames = list(LNAME, levels(g)[1]))
        PVAL <- matrix(data=pvalv, nrow = (k-1), ncol = 1,
                   dimnames = list(LNAME, levels(g)[1]))
        MODEL <- data.frame(x = x, g = g)
        ## colnames(PVAL) <- grpn[1]
        ## colnames(PSTAT) <- grpn[1]
        ## rownames(PVAL) <- grpn[2:k]
        ## rownames(PSTAT) <- grpn[2:k]
        ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL,
                    statistic = PSTAT, p.adjust.method = p.adjust.method,
                    dist = "z", model = MODEL, alternative = alternative)
        class(ans) <- "PMCMR"
        ans
}

#' @rdname kwManyOneNdwTest
#' @method kwManyOneNdwTest formula
#' @aliases kwManyOneNdwTest.formula
#' @template one-way-formula
#' @export
kwManyOneNdwTest.formula <-
    function(formula, data, subset, na.action,
             alternative = c("two.sided", "greater", "less"),
             p.adjust.method = c("single-step", p.adjust.methods), ...)
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
    alternative <- match.arg(alternative)
    p.adjust.method <- match.arg(p.adjust.method)
    y <- do.call("kwManyOneNdwTest", c(as.list(mf),
                               alternative = alternative,
                               p.adjust.method = p.adjust.method))
    y$data.name <- DNAME
    y
}
