## powerOneWayTests.R
## Part of the R package: PMCMR
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


#' @name powerOneWayTests
#' @title Power Simulation for One-Factorial Single Hypothesis Tests
#' @description
#' Performs power simulation for one-factorial
#' single hypothesis tests.
#'
#' @param mu numeric vector of group means.
#' @param n number of replicates per group. If \code{n} is a scalar, then
#' a balanced design is assumed. Otherwise, \code{n} must be a vector of same
#' length as \code{mu}.
#' @param errfn the error function. Defaults to \code{"Normal"}.
#' @param parms a list that denotes the arguments for the error function.
#' Defaults to \code{list(mean=0, sd=1)}.
#' @param test the test for which the power analysis is
#' to be performed. Defaults to \code{"kwManyOneConoverTest"}.
#' @param alternative the alternative hypothesis. Defaults to \code{"two.sided"},
#' ignored if the selected error function does not use this argument.
#' @param var.equal a logical variable indicating whether to treat the variances
#'          in the samples as equal.  \code{"TRUE"}, then a simple F test for
#'          the equality of means in a one-way analysis of variance is
#'          performed.  If \code{"FALSE"}, an approximate method of Welch (1951)
#'          is used, which generalizes the commonly known 2-sample Welch
#'          test to the case of arbitrarily many samples. Defaults to \code{"TRUE"}; only relevant,
#' if \code{test = "oneway.test"}, otherwise ignored.
#' @param dist the test distribution. Only relevant for
#' \code{\link{kruskalTest}}. Defaults's to \code{NULL}.
#' @param alpha the nominal level of Type I Error.
#' @param FWER logical, indicates whether the family-wise error should be computed.
#' Defaults to \code{TRUE}.
#' @param replicates the number of Monte Carlo replicates or runs. Defaults to \code{1000}.
#' @param p the a-priori known peak as an ordinal number of the treatment
#' group including the zero dose level, i.e. \eqn{p = \{1, \ldots, k\}}.
#' Defaults to \code{NULL}. Only relevant, if \code{"mackWolfeTest"} is selected.
#' @details
#' The linear model of a one-way ANOVA can be written as: 
#' 
#' \deqn{
#' X_{ij} = \mu_i + \epsilon_{ij}
#' }
#'
#' For each Monte Carlo run, the function simulates \eqn{\epsilon_{ij}} based on the given error function and
#' the corresponding parameters. Then the specified test is performed.
#' Finally, Type I and Type II error rates are calculated.
#'
#' @return
#' An object with class \code{powerOneWayPMCMR}.
#'
#' @examples
#' \dontrun{
#' set.seed(12)
#' mu <- c(0, 0, 1, 2)
#' n <- c(5, 4, 5, 5)
#' parms <- list(mean=0, sd=1)
#' powerOneWayTests(mu, n, parms, test = "cuzickTest",
#' alternative = "two.sided", replicates = 1E4)
#'
#' ## Compare power estimation for
#' ## one-way ANOVA with balanced design
#' ## as given by functions
#' ## power.anova.test, pwr.anova.test
#' ## and powerOneWayTest
#' 
#' groupmeans <- c(120, 130, 140, 150)
#' SEsq <- 500  # within-variance
#' n <- 10
#' k <- length(groupmeans)
#' df <- n * k - k
#' SSQ.E <- SEsq * df
#' SSQ.A <- n * var(groupmeans) * (k - 1)
#' sd.errfn <- sqrt(SSQ.E / (n * k - 1))
#' R2 <- c("R-squared" = SSQ.A / (SSQ.A + SSQ.E))
#' cohensf <- sqrt(R2 / (1 - R2))
#' names(cohensf) <- "Cohens f"
#'
#' ## R stats power function 
#' power.anova.test(groups = k,
#'                  between.var = var(groupmeans),
#'                  within.var = SEsq,
#'                  n = n)
#'
#' ## pwr power function
#' pwr.anova.test(k = k, n = n, f = cohensf, sig.level=0.05)
#'
#' ## this Monte-Carlo based estimation  
#' set.seed(200)
#' powerOneWayTests(mu = groupmeans,
#'                  n = n,
#'                  parms = list(mean=0, sd=sd.errfn),
#'                  test = "oneway.test",
#'                  var.equal = TRUE,
#'                  replicates = 5E3)
#'
#' ## Compare with effect sizes
#' R2
#' cohensf
#'
#' }
#' 
#' @seealso
#' \code{\link{powerMCTests}},
#' \code{\link[pwr]{pwr.anova.test}}, 
#' \code{\link[stats]{power.anova.test}}
#' @importFrom stats runif
#' @importFrom stats oneway.test
#' @importFrom stats as.formula
#' @concept TestPower
#' @keywords misc
#' @export
powerOneWayTests <- function(mu,
                         n = 10,
                         errfn = c("Normal",
                                   "Lognormal",
                                   "Exponential",
                                   "Chisquare",
                                   "TDist",
                                   "Cauchy",
                                   "Weibull"),
                         parms = list(mean=0, sd = 1),                         
                         test = c("kruskalTest",
                                  "leTest",
                                  "vanWaerdenTest",
                                  "normalScoresTest",
                                  "spearmanTest",
                                  "cuzickTest",
                                  "jonckheereTest",
                                  "johnsonTest",
                                  "oneway.test",
                                  "adKSampleTest",
                                  "bwsKSampleTest",
                                  "bwsTrendTest",
                                  "mackWolfeTest"),
                         alternative = c("two.sided", "greater", "less"),
                         var.equal = TRUE,
                         dist = NULL, 
                         alpha = 0.05,
                         FWER = TRUE,
                         replicates=1000,
                         p = NULL)
{

    test <- match.arg(test)
    alternative <- match.arg(alternative)
    errfn <- match.arg(errfn)
    
    if (!is.vector(mu)){
        stop("'mu' must be a vector of type numeric")
    }
    k <- length(mu)
    if (k < 3){
        stop("'mu' must be a vector with at least 3 elements")
    } else if (k > 14){
        stop("'mu' must be a vector with less than 15 elements")
    } else if (replicates > 1e5){
        stop("'replicates' must be less than 1e5")
    }
    if (is.vector(n)){
        if (k == length(n)){
            ni <- n
	} else {
	    ni <- rep(n, k)
	}
    }

    ## group vector
    g <- unlist(sapply(1:k, function(i) rep(i, ni[i])))
    g <- as.factor(g)

    ## Sample size
    N <- sum(ni)

    ##
    loc <- as.vector(unlist(sapply(1:k, function(i) rep(mu[i], ni[i]))))
    if(FWER){
        loc0 <- rep(0, N)
    }

    ## create uniform random numbers
    PUNIF <- matrix(data=runif(N * replicates),
                    nrow=replicates, ncol=N, byrow=TRUE)

    if (errfn == "Normal"){
        FN <- "qnorm"
        if (is.vector(parms$sd) & length(parms$sd) ==k) {
            sd <- as.vector(unlist(sapply(1:k, function(i)
                rep(parms$sd[i], ni[i]))))
        } else if (is.numeric(parms$sd)) {
            sd = rep(parms$sd, N)
        } else {
            sd = rep(1, N)
        }
        fnargs <- list(mean = rep(parms$mean, N), sd = sd)
    } else if (errfn == "Lognormal"){
        FN <- "qlnorm"
        if (is.vector(parms$sdlog) & length(parms$sdlog) ==k) {
            sdlog <- as.vector(unlist(sapply(1:k, function(i)
                rep(parms$sdlog[i], ni[i]))))
        } else if (is.numeric(parms$sdlog)) {
            sdlog = rep(parms$sdlog, N)
        } else {
            sdlog = rep(1, N)
        }
        fnargs <- list(meanlog = rep(parms$meanlog, N), sdlog = sdlog)
    } else if (errfn == "Exponential"){
        FN <- "qexp"
        if (is.vector(parms$rate) & length(parms$rate) ==k) {
            rate <- as.vector(unlist(sapply(1:k, function(i)
                rep(parms$rate[i], ni[i]))))
        } else if (is.numeric(parms$rate)) {
            rate = rep(parms$rate, N)
        } else {
            rate = rep(1, N)
        }
        fnargs <- list(rate=rate)
    } else if (errfn == "Chisquare"){
        FN <- "qchisq"
        if (is.vector(parms$df) & length(parms$df) ==k) {
            df <- as.vector(unlist(sapply(1:k, function(i)
                rep(parms$df[i], ni[i]))))
        } else if (is.numeric(parms$df)) {
            df = rep(parms$df, N)
        } else {
            df = rep(k-1, N)
        }
        fnargs <- list(df = df)
    } else if (errfn == "TDist"){
        FN <- "qt"
        if (is.vector(parms$df) & length(parms$df) ==k) {
            df <- as.vector(unlist(sapply(1:k, function(i)
                rep(parms$df[i], ni[i]))))
        } else if (is.numeric(parms$df)) {
            df = rep(parms$df, N)
        } else {
            df = rep(k, N)
        }
        fnargs <- list(df = df)
    } else if (errfn == "Cauchy"){
        FN <- "qcauchy"
        if (is.vector(parms$scale) & length(parms$scale) ==k) {
            scale <- as.vector(unlist(sapply(1:k, function(i)
                rep(parms$scale[i], ni[i]))))
        } else if (is.numeric(parms$scale)) {
            scale = rep(parms$scale, N)
        } else {
            scale = rep(1, N)
        }
        fnargs <- list(location = parms$location, scale = scale)
    } else if (errfn == "Weibull"){
        FN <- "qweibull"
        if (is.vector(parms$scale) & length(parms$scale) ==k) {
            scale <- as.vector(unlist(sapply(1:k, function(i)
                rep(parms$scale[i], ni[i]))))
        } else if (is.numeric(parms$scale)) {
            scale = rep(parms$scale, N)
        } else {
            scale = rep(1, N)
        }
        fnargs <- list(shape = parms$shape, scale = scale)
    } 

    m <- 1
    ## simulation
    if (FWER) {
        mat0 <- matrix(mat0 <- sapply(1:replicates,
                                      function(j){
                                          fnargs$p <- PUNIF[j,]
                                          X <- loc0 + do.call(FN, fnargs)
                                          if(test == "oneway.test"){
                                              ans <- do.call(test,
                                                             args=list(
                                                                 as.formula( X ~ g),
                                                                 var.equal = var.equal)
                                                             )$p.value
                                          } else {
                                              ans <- do.call(test,
                                                             args=list(
                                                                 as.formula( X ~ g),
                                                                 alternative = alternative,
                                                                 dist=dist,
                                                                 p = p)
                                                             )$p.value
                                          }
                                          ans
                                      }
                                      )
                      ,nrow=replicates, ncol=m, byrow=TRUE)
    } else {
        mat0 <-NULL
    }
        
    mat <- matrix(mat <- sapply(1:replicates,
                                function(j){
                                    fnargs$p <- PUNIF[j,]
                                    X <- loc + do.call(FN, fnargs)

                                    if(test == "oneway.test"){
                                        ans <- do.call(test,
                                                       args=list(
                                                           as.formula( X ~ g),
                                                           var.equal = var.equal)
                                                       )$p.value
                                    } else {
                                        ans <- do.call(test,
                                                       args=list(
                                                           as.formula( X ~ g),
                                                           alternative = alternative,
                                                           dist=dist,
                                                           p = p)
                                                       )$p.value
                                    }
                                    ans
                                }
                                )
                 ,nrow=replicates, ncol=m, byrow=TRUE)
    
    ## If mat0 -> family wise error rate = Type I Error
    ## If mat -> Power = 1 - Type II Error 
    fwerfn <- function(inp)
    {
        ok <- sapply(1:replicates, function(i) any(inp[i,] < alpha))
        tmp2 <- 1:replicates
        V <- length(tmp2[ok])
        return(V/replicates)
    }


    power <- fwerfn(mat)
    if(FWER) {
        fwer <- fwerfn(mat0)
    } else {
        fwer <- NA
    }

    ans <- list(mu = mu,
                parms = parms,
                n = n,
                errfn = errfn,
                test = test,
                alternative = alternative,
                alpha = alpha,
                replicates = replicates,
                fwer = fwer,
                power = power,
                p0 = mat0,
                p = mat)
    class(ans) <- "powerOneWayPMCMR"
    ans
}
