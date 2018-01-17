## adKSampleTest.R
## 
## Copyright (C) 2017 Thorsten Pohlert
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
#' @name adKSampleTest
#' @title Anderson-Darling k-Sample Test
#' 
#' @description
#' Performs Anderson-Darling k-sample test.
#' 
#' @details
#' The null hypothesis, H\eqn{_0: F_1 = F_2 = \ldots = F_k}
#' is tested against the alternative,
#' H\eqn{_\mathrm{A}: F_i \ne F_j ~~(i \ne j)}, with at least
#' one unequality beeing strict.
#'
#' This function only evaluates version 1 of the k-sample Anderson-Darling
#' test (i.e. Eq. 6) of Scholz and Stephens (1987).
#' The p-values are estimated with the extended empirical function
#' as implemented in \code{\link[kSamples]{ad.pval}} of
#' the package \pkg{kSamples}.
#' 
#' @template class-htest
#'
#' @seealso
#' \code{\link{adAllPairsTest}}, \code{\link{adManyOneTest}},
#' \code{\link[kSamples]{ad.pval}}.
#' @examples
#' ## Hollander & Wolfe (1973), 116.
#' ## Mucociliary efficiency from the rate of removal of dust in normal
#' ## subjects, subjects with obstructive airway disease, and subjects
#' ## with asbestosis.
#' x <- c(2.9, 3.0, 2.5, 2.6, 3.2) # normal subjects
#' y <- c(3.8, 2.7, 4.0, 2.4)      # with obstructive airway disease
#' z <- c(2.8, 3.4, 3.7, 2.2, 2.0) # with asbestosis
#'
#' datf <- data.frame(
#'   g = g <- c(rep("ns", length(x)), rep("oad",
#'       length(y)), rep("a", length(z))),
#'   x = x <- c(x, y, z))
#'
#' adKSampleTest(x ~ g, datf)
#' 
#' @references
#' Scholz, F.W., Stephens, M.A. (1987) K-Sample Anderson-Darling Tests.
#' \emph{Journal of the American Statistical Association}, 82, 9118--924.
#' @export adKSampleTest
adKSampleTest <- function(x, ...) UseMethod("adKSampleTest")

#' @rdname adKSampleTest
#' @method adKSampleTest default
#' @aliases adKSampleTest.default
#' @template one-way-parms
#' @importFrom stats complete.cases
#' @importFrom kSamples ad.pval
#' @export
adKSampleTest.default <-
    function(x, g, ...)
{
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
  ##      dist <- x$dist
  ##      replicates <- x$replicates
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

##    dist <- match.arg(dist)
    ix <- order(as.character(g))
    g <- g[ix]
    x <- x[ix]
    
    ## prepare
    n <- tapply(x, g, length)
    N <- sum(n)

#    ## Function ADstat without ties
#    ADstat <- function(x, g){      
#        o <- rank(x)
#        ## local variables
#        M <- matrix(NA, ncol=N, nrow=k)
#        Z <- 1:N
#        lev <- levels(g)
#        for (i in 1:k){
#            M[i,] <- sapply(Z, function(z)
#                length(o[g == lev[i] & o <= z]))
#        }
#        
#        tmp1 <- 0
#        j <- 1:(N-1)
#        for (i in 1:k){
#            tmp2 <- sum( (N * M[i,j] - j * n[i])^2 / (j * (N - j)))            # 
#            tmp1 <- tmp1 + 1/n[i] * tmp2
#        }
#        AsqkN <- 1/N * tmp1
#        names(AsqkN) <- NULL
#        return(AsqkN)
#    }

    ## Function ADstatV1 with ties
    ADstatV1 <- function(x, g){
        lev <- levels(g)
        Zstar <- sort(x)
        Zstar <- unique(Zstar)
        L <- length(Zstar)
        f <- matrix(0, ncol=L, nrow=k)
        for (i in 1:k){
            tmp <- x[g == lev[i]]
            f[i,] <- sapply(Zstar, function(Z){
                return(length(tmp[tmp == Z]))
                }
                )
        }

        l <- sapply(1:L, function(j) sum(f[,j]))
        
        tmp <- rep(NA,k)
        for (i in 1:k){
            tm <- 0
            for (j in 1:(L-1)){
                Bj <- sum(l[1:j])
                Mij <- sum(f[i, 1:j])
                tm <- tm + (l[j] / N) *
                    ((N * Mij - n[i] * Bj)^2 /
                    (Bj * (N - Bj)))
            }
            tmp[i] <- tm
        } 
        AsqkN <- sum(1/n * tmp)
        return(AsqkN)
    }
    
    ## Check for ties
    ## ties <- sum(table(x) - 1)

 #   if (ties == 0){
#        AsqkN <- ADstat(x, g)
#    } else {
    AsqkN <- ADstatV1(x, g)
#    }
    
    ## Calculate varAsqkN
    H <- sum(1/n)
    h <- sum(1/(1:(N-1)))
    G <- 0
    for (i in 1:(N-2)){
        G <- G + sum( 1 / ((N - i) * (i+1):(N-1)))
    }
    
    ## coefficients
    a <- (4 * G - 6) * (k - 1) + (10 - 6 * G) * H
    b <- (2 * G - 4) * k^2 + 8 * h * k +
        (2 * G - 14 * h - 4) * H - 8 * h + 4 * G - 6
    c <- (6 * h + 2 * G - 2) *  k^2 +
        (4 * h - 4 * G + 6) * k +
        (2 * h - 6) * H + 4 * h
    d <- (2 * h + 6) * k^2 - 4 * h * k

    varAsqkN <- (a * N^3 + b * N^2 + c * N + d) /
        ((N - 1) * (N - 2) * (N- 3))
    
    m <- k - 1
    TkN <- (AsqkN - m) / sqrt(varAsqkN)

    ## p-value from package kSample
    PVAL <- ad.pval(tx=TkN, m=m, version=1)
    
    ## Interpolation coefficients according to
    ## Scholz and Stephens, Tab 2, p. 920
 ###   alpha <- c(0.25, 0.1, 0.05, 0.025, 0.01)
 ##   b0 <- c(0.675, 1.281, 1.645, 1.960, 2.326)
 ##   b1 <- c(-0.245, 0.250, 0.678, 1.149, 1.822)
 ##   b2 <- c(-0.105, -0.305, -0.362, -0.391, - 0.396)
##
 ##   tm <- b0 + b1 / sqrt(m) + b2 / m
 ##   mod <- lm(log(alpha) ~ tm)
##
##    if (dist != "asymptotic"){
##        PVAL <- exp(predict(mod, newdata=data.frame(tm=TkN)))
##        PVAL <- min(1, PVAL)
##    } else {
##        AA <- replicate(replicates, ADstatV1(rnorm(N), g))
##        t.AA <- (AA - m) / sqrt(varAsqkN)
##     ##   gt.TkN <- sum(ifelse(t.AA > TkN, 1, 0))
##        gt.TkN <- length(t.AA[t.AA > TkN])
##        PVAL <- gt.TkN / replicates
##    }
##    
    METHOD <- paste("Anderson-Darling k-sample test")
	
    ans <- list(method = METHOD,
                data.name = DNAME,
                p.value = PVAL,
                statistic = c(TkN = TkN),
                parameter = c(m = m),
                estimate = c(A2kN = AsqkN, "sigmaN" = sqrt(varAsqkN)))
    class(ans) <- "htest"
    ans
}

#' @rdname adKSampleTest
#' @method adKSampleTest formula
#' @aliases adKSampleTest.formula
#' @template one-way-formula
#' @export
adKSampleTest.formula <-
function(formula, data, subset, na.action, ...)
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
    y <- do.call("adKSampleTest", as.list(mf))
    y$data.name <- DNAME
    y
}
