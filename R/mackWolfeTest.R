## mackWolfeTest.R
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

#' @rdname mackWolfeTest
#' @title Mack-Wolfe Test for Umbrella Alternatives
#' @description
#' Performs Mack-Wolfe non-parametric test for umbrella alternatives.
#'
#' @details
#' In dose-finding studies one may assume an increasing treatment
#' effect with increasing dose level. However, the test
#' subject may actually succumb to toxic effects at high doses,
#' which leads to decresing treatment effects.
#'
#' The scope of the Mack-Wolfe Test is to test for umbrella alternatives
#' for either a known or unknown point \eqn{p} (i.e. dose-level),
#' where the peak (umbrella point) is present.
#'
#' H\eqn{_i: \theta_0 = \theta_i = \ldots = \theta_k} is tested
#' against the alternative A\eqn{_i: \theta_1 \le \ldots \theta_p \ge
#' \theta_k} for some \eqn{p}, with at least one strict inequality.
#'
#' If \code{p = NULL} (peak unknown), the upper-tail \eqn{p}-value is computed
#' via an asymptotic bootstrap permutation test.
#'
#' If an integer value for \code{p} is given (peak known), the
#' upper-tail \eqn{p}-value is computed from the standard normal
#' distribution (\code{\link{pnorm}}).
#'
#' @note
#' One may increase the number of permutations to e.g. \code{nperm = 10000}
#' in order to get more precise p-values. However, this will be on
#' the expense of computational time.
#'
#' @references
#' Chen, I. Y. (1991) Notes on the Mack-Wolfe and Chen-Wolfe
#' Tests for Umbrella Alternatives, \emph{Biom. J.} \bold{33}, 281--290.
#'
#' Mack, G. A., Wolfe, D. A. (1981) K-sample rank tests for
#' umbrella alternatives, \emph{J. Amer. Statist. Assoc.} \bold{76}, 175--181.
#'
#' @template class-htest
#' @examples
#' ## Example from Table 6.10 of Hollander and Wolfe (1999).
#' ## Plates with Salmonella bacteria of strain TA98 were exposed to
#' ## various doses of Acid Red 114 (in mu g / ml).
#' ## The data are the numbers of visible revertant colonies on 12 plates.
#' ## Assume a peak at D333 (i.e. p = 3).
#' x <- c(22, 23, 35, 60, 59, 54, 98, 78, 50, 60, 82, 59, 22, 44,
#'   33, 23, 21, 25)
#' g <- as.ordered(rep(c(0, 100, 333, 1000, 3333, 10000), each=3))
#' plot(x ~ g)
#' mackWolfeTest(x=x, g=g, p=3)
#'
#' @concept umbrellatest
#' @concept kruskalranks
#' @keywords htest nonparametric
#' @seealso
#' \code{\link{pnorm}}, \code{\link{sample}}.
#' @export
mackWolfeTest <- function(x, ...) UseMethod("mackWolfeTest")

#' @rdname mackWolfeTest
#' @aliases mackWolfeTest.default
#' @method mackWolfeTest default
#' @template one-way-parms
#' @param p the a-priori known peak as an ordinal number of the treatment
#' group including the zero dose level, i.e. \eqn{p = \{1, \ldots, k\}}.
#' Defaults to \code{NULL}.
#' @param nperm number of permutations for the assymptotic permutation test.
#' Defaults to \code{1000}.
## @importFrom mvtnorm pmvnorm
#' @importFrom stats complete.cases
#' @importFrom stats pnorm
#' @export
mackWolfeTest.default <-
    function(x, g, p = NULL, nperm=1000,...)
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
        p <- x$p
        nperm <- x$nperm
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

    ## check p
    if (!is.null(p)){
        if (p > k){
            stop("Selected 'p' > Nr. of groups: ", p, " > ", k)
        } else if  (p < 1){
            stop("Selected 'p' < 1: ", p)
        }
    }

    ## rank the data
    Rij <- rank(x)
    n <- tapply(x, g, length)

    ## U statistic
    Ustat.fn <- function(Rij, g, k){
        lev <- levels(g)
        U <- diag(k)

        .fn <- function(Ri, Rj){
            tmp <- sum(
                unlist(
                    sapply(Ri, function(i)
                        length(Rj[Rj > i])
                        )
                )
            )
            return(tmp)
        }
        for(i in 2:k){
            for(j in 1:(i-1)){
                U[i,j] <- .fn(Rij[g==lev[i]], Rij[g==lev[j]])
                U[j,i] <- .fn(Rij[g==lev[j]], Rij[g==lev[i]])
            }
        }
        return(U)
    }

    ## Calculate Ap
    Ap.fn <- function(p, U){
        tmp1 <- 0
        if (p > 1){
            for(i in 1:(p-1)){
                for (j in (i+1):p){
                    tmp1 <- tmp1 + U[i,j]
                }
            }
        }
        tmp2 <- 0
        if (p < k){
            for (i in p:(k-1)){
                for(j in (i+1):k){
                    tmp2 <- tmp2 + U[j,i]
                }
            }
        }
        return(tmp1 + tmp2)
    }

    ## N1 and N2 function
    N1.fn <- function(p, n) sum(n[1:p])
    N2.fn <- function(p, n) sum(n[p:k])

    ## E0(At)
    meanAt <- function(p, n){
        N1 <- N1.fn(p, n)
        N2 <- N2.fn(p, n)
        return((N1^2 + N2^2 - sum(n^2) - n[p]^2)/4)
    }

    ## Var0(At)
    varAt <- function(p, n){
        N1 <- N1.fn(p, n)
        N2 <- N2.fn(p, n)
        N <- sum(n)
        var0 <- (2 * (N1^3 + N2^3) +
                 3 * (N1^2 + N2^2) -
                 sum(n^2 * (2 * n + 3)) - n[p]^2 *
                 (2 * n[p] + 3) + 12 * n[p] *
                 N1 * N2 - 12 * n[p]^2 * N) / 72
        return(var0)
    }


    if (!is.null(p)){

        ## This is Mack-Wolfe Test for known p
        ## check for ties
        if(sum(table(x) - 1) > 0){
            warning(c("Ties are present. Exact variances can not be computed.",
                      " p-values are over-estimated."))
        }
        U <- Ustat.fn(Rij, g, k)
        EST <- c("Ap" = Ap.fn(p, U))
        MEAN <- meanAt(p, n)
        SD <- sqrt(varAt(p,n))
       ## ESTIM <- levels(g)[p]
        METHOD <- c("Mack-Wolfe test for umbrella alternatives",
                    "with known peak")
        STAT <- (EST - MEAN)/ SD
        names(STAT) <- NULL
        names(STAT) <- "z"
        PVAL <- pnorm(STAT, lower.tail=FALSE)

    } else {
        ## This is Chen's modified version for unknown p
        U <- Ustat.fn(Rij, g, k)
        Ap <- sapply(1:k, function(p) Ap.fn(p, U))
        mean0 <- sapply(1:k, function(p) meanAt(p, n))
        var0 <- sapply(1:k, function(p) varAt(p, n))
        Astar <- (Ap - mean0) / sqrt(var0)
        STAT <- c("Ap*" = max(Astar))
        p <- which(Astar == STAT)
        EST <- NULL
        ## asymptotic bootstrap permutation
        mt <- sapply(1:nperm, function(i){
            ix <- sample(Rij)
            Uix <- Ustat.fn(ix, g, k)
            Apix <- sapply(1:k, function(p) Ap.fn(p, Uix))
            Astarix <- (Apix - mean0) / sqrt(var0)
            max(Astarix)
            })

        PVAL <- length(mt[mt > STAT]) / nperm
        METHOD <- c("Mack-Wolfe test for umbrella alternatives",
                    "with unknown peak")
    }
        ## Covariance Matrix Cov(As,At)
        ## cov0 <- diag(k-1)
       ## Ni <- sapply(2:k, function(i) sum(n[1:i]))
      ##  N <- sum(n)
      ##  for(s in (2:(k-2))){
      ##      for(t in ((s+1):(k-1))){
      ##
      ##          tmp0 <- 0
      ##          for (i in (2:s)){
      ##              tmp0 <- tmp0 + (n[i] * Ni[i-1] + (Ni[i] + 1))
      ##          }
      ##
      ##          tmp1 <- 0
      ##          for(j in(t:(k-1))){
      ##              tmp1 <- tmp1 + n[j] * (N - Ni[j]) * (N - Ni[j-1] + 1)
      ##          }
      ##          tmp2 <- 0
      ##          for (i in (s:t)){
      ##              n[i] * Ni[i-1] * (Ni[i] + 1) + Ni[s-1] * (Ni[t] - Ni[s-1]) * (N + 1)
      ##          }
      ##
      ##          cov0[(s-1),(t-1)] <- (tmp0 + tmp1 - tmp2) / 12 * sqrt(var0[s] * var0[t])
      ##          cov0[(t-1),(s-1)] <- cov0[(s-1),(t-1)]
      ##      }
      ##  }
      ##
      ##  PVAL <-   1 - pmvnorm(lower = -rep(abs(STAT), k-1),
      ##                        upper = rep(abs(STAT), k-1),
      ##                        sigma = cov0)

  ##  }

    ans <- list(method = METHOD,
                data.name = DNAME,
                p.value = PVAL,
                statistic = STAT,
                alternative = paste0("\ntheta_1 <= ... <= theta_p >= ... >= theta_k, p = ",p),
                estimates = EST)
    class(ans) <- "htest"
    ans
}

#' @rdname mackWolfeTest
#' @aliases mackWolfeTest.formula
#' @method mackWolfeTest formula
#' @template one-way-formula
#' @export
mackWolfeTest.formula <-
    function(formula, data, subset, na.action, p = NULL, nperm=1000, ...)
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
    y <- do.call("mackWolfeTest", c(as.list(mf), p = p, nperm=nperm))
    y$data.name <- DNAME
    y
}

