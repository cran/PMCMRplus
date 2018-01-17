## powerMCTests.R
## Part of the R package: PMCMRplus
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

#' @name powerMCTests
#' @title Power Simulation for One-Factorial All-Pairs and Many-To-One Comparison Tests
#' @description
#' Performs power simulation for one-factorial all-pairs and Many-To-One comparison tests.
#'
#' @param mu numeric vector of group means.
#' @param n number of replicates per group. If \code{n} is a scalar, then
#' a balanced design is assumed. Otherwise, \code{n} must be a vector of same
#' length as \code{mu}.
#' @param errfn the error function. Defaults to \code{"Normal"}.
#' @param parms a list that denotes the arguments for the error function.
#' Defaults to \code{list(mean=0, sd=1)}.
#' @param test the multiple comparison test for which the power analysis is
#' to be performed. Defaults to \code{"kwManyOneConoverTest"}.
#' @param alternative the alternative hypothesis. Defaults to \code{"two.sided"},
#' ignored if the selected error function does not use this argument.
#' @param p.adjust.method method for adjusting p values (see \code{\link{p.adjust}}).
#' @param alpha the nominal level of Type I Error.
#' @param FWER logical, indicates whether the family-wise error should be computed.
#' Defaults to \code{TRUE}.
#' @param replicates the number of Monte Carlo replicates or runs. Defaults to \code{1000}.
#'
#' @details
#' The linear model of a one-way ANOVA can be written as: 
#' 
#' \deqn{
#' X_{ij} = \mu_i + \epsilon_{ij}
#' }
#'
#' For each Monte Carlo run, the function simulates \eqn{\epsilon_{ij}} based on the given error function and
#' the corresponding parameters. Then the specified all-pairs
#' or many-to-one comparison test is performed.
#' Finally, several effect sizes (Cohen's f ans R-squared),
#' error rates (per comparison error rate,
#' false discovery rate and familywise error rate)
#' and test powers (any-pair power, average per-pair power
#' and all-pairs power) are calculated.
#' 
#' @return
#' An object with class \code{powerPMCMR}.
#'
#' @examples
#' \dontrun{
#' mu <- c(0, 0, 1, 2)
#' n <- c(5, 4, 5, 5)
#' set.seed(100)
#' powerMCTests(mu, n, errfn="Normal",
#'  parms=list(mean=0, sd=1),
#'  test="dunnettTest", replicates=1E4)
#'
#' powerMCTests(mu, n, errfn="Normal",
#'  parms=list(mean=0, sd=1),
#'  test="kwManyOneDunnTest", p.adjust.method = "bonferroni",
#'  replicates=1E4)
#' 
#' }
#'
#' @importFrom stats pairwise.t.test
#' @importFrom stats runif
#' @importFrom stats qnorm
#' @importFrom stats qlnorm
#' @importFrom stats qexp
#' @importFrom stats qchisq
#' @importFrom stats qt
#' @importFrom stats qcauchy
#' @importFrom stats qweibull
#' @concept TestPower
#' @keywords misc
#' @export
powerMCTests <- function(mu,
                         n = 10,
                         errfn = c("Normal",
                                   "Lognormal",
                                   "Exponential",
                                   "Chisquare",
                                   "TDist",
                                   "Cauchy",
                                   "Weibull"),
                         parms = list(mean=0, sd = 1),                         
                         test = c("kwManyOneConoverTest",
                                  "kwManyOneDunnTest",
                                  "kwManyOneNdwTest",
                                  "vanWaerdenManyOneTest",
                                  "normalScoresManyOneTest",
                                  "dunnettTest",
                                  "tamhaneDunnettTest",
                                  "ManyOneUTest",
                                  "kwAllPairsNemenyiTest",
                                  "kwAllPairsDunnTest",
                                  "kwAllPairsConoverTest",
                                  "normalScoresAllPairsTest",
                                  "vanWaerdenAllPairsTest",
                                  "dscfAllPairsTest",
                                  "gamesHowellTest",
                                  "lsdTest",
                                  "scheffeTest", 
                                  "tamhaneT2Test",
                                  "tukeyTest",
                                  "dunnettT3Test",
                                  "pairwise.t.test",
                                  "pairwise.wilcox.test",
                                  "adManyOneTest",
                                  "adAllPairsTest",
                                  "bwsManyOneTest",
                                  "bwsAllPairsTest"),
                         alternative = c("two.sided", "greater", "less"),
                         p.adjust.method= c("single-step", p.adjust.methods), 
                         alpha = 0.05,
                         FWER = TRUE,
                         replicates=1000)
{

    test <- match.arg(test)
    alternative <- match.arg(alternative)
    p.adjust.method <- match.arg(p.adjust.method)
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

    ## check tests
    ## ManyOne Tests
    ManyOne <- c("kwManyOneConoverTest",
              "kwManyOneDunnTest",
              "kwManyOneNdwTest",
              "vanWaerdenManyOneTest",
              "normalScoresManyOneTest",
              "dunnettTest",
              "tamhaneDunnettTest",
              "ManyOneUTest",
              "adManyOneTest",
              "bwsManyOneTest")

    ## AllPairs Tests
    AllPairs <- c("kwAllPairsNemenyiTest",
              "kwAllPairsDunnTest",
              "kwAllPairsConoverTest",
              "normalScoresAllPairsTest",
              "vanWaerdenAllPairsTest",
              "dscfAllPairsTest",
              "gamesHowellTest",
              "lsdTest",
              "scheffeTest", 
              "tamhaneT2Test",
              "tukeyTest",
              "dunnettT3Test",
              "pairwise.t.test",
              "pairwise.wilcox.test",
              "adAllPairsTest",
              "bwsAllPairsTest")
    
    mt1 <- any(sapply(1:length(ManyOne), function(i) (ManyOne[i] == test)))
    mtm <- any(sapply(1:length(AllPairs), function(i) (AllPairs[i] == test)))
    if (mt1){
        m <- k - 1
    } else if (mtm) {
        m <- k * (k - 1) / 2
    } else {
        stop("Could not find 'test':", test)
    }
    
    H <- logical(m)
    mm <- 0
    if (mt1){
        for (j in 2:k){
            mm <- mm + 1
            H[mm] <- (mu[1] == mu[j])
        }
    } else if(mtm) {
        for (i in 1:(k-1)){
            for (j in (i+1):k){
                mm <- mm + 1
                H[mm] <- (mu[i] == mu[j])
            }
        }
    }
    
    ## nr of Hypothesis with H0 true
    tmp <- 1:m
    m0 <- length(tmp[H])
	
    ## group vector
    g <- unlist(sapply(1:k, function(i) rep(i, ni[i])))
    g <- as.factor(g)

    N <- sum(ni)
    
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

    if (FWER) {
        mat0 <- matrix(mat0 <- sapply(1:replicates,
                                      function(j){
                                          fnargs$p <- PUNIF[j,]
                                          X <- loc0 + do.call(FN, fnargs)
                                          ans <- get.pvalues(
                                              do.call(test,
                                                      args=list(
                                                          x=X,
                                                          g=g,
                                                          alternative = alternative,
                                                          p.adjust.method=p.adjust.method,
                                                          pool.sd = TRUE)
                                                      )
                                              )
                                          ans
                                      }
                                      )
                      ,nrow=replicates, ncol=m, byrow=TRUE)
    } else {
        mat0 <- NULL
    }
    
    mat <- matrix(mat <- sapply(1:replicates,
                                function(j){
                                    fnargs$p <- PUNIF[j,]
                                    X <- loc + do.call(FN, fnargs)
                                    ans <- get.pvalues(
                                        do.call(test,
                                                args=list(
                                                    x=X,
                                                    g=g,
                                                    alternative = alternative,
                                                    p.adjust.method=p.adjust.method,
                                                    pool.sd = TRUE)
                                                )
                                    )
                                    ans
                                }
                                )
                 ,nrow=replicates, ncol=m, byrow=TRUE)
    
    ## returns a vector with counted Pr(H0|H) <= alpha
    indicatorFn <- function(inp, type=c("any", "all", "sum"), alpha)
    {
        if (!is.matrix(inp)){
            inp <- cbind(inp)
        }
        type <- match.arg(type)

        ##
        r <- replicates
        if (type == "any"){
            out <- sapply(1:r, function(i)
                ifelse(any(inp[i,] <= alpha), 1, 0)
                )

        } else if (type == "all") {
            out <- sapply(1:r, function(i)
                ifelse(all(inp[i,] <= alpha), 1, 0)
                )
        } else {
            out <- sapply(1:r, function(i)
                sum( ifelse(inp[i,] <= alpha, 1, 0))
                )
        }
        return(out)
    }
    
    ## Power
    m1 <- m - m0
    if (m1 == 0){
        avep = 0
        allp <- 0
        anyp <- 0
    } else {
        tmp <-  indicatorFn(inp=mat[,!H], type="sum", alpha=alpha)
        avep <- sum(tmp) / length(tmp) / m1

        tmp <- indicatorFn(inp=mat[,!H], type="any", alpha=alpha)
        anyp <- sum(tmp) / length(tmp)

        tmp <- indicatorFn(inp=mat[,!H], type="all", alpha=alpha)
        allp <- sum(tmp) / length(tmp)
    }

    ## FWER for mat0
    vv <- indicatorFn(inp=mat0, type="any", alpha=alpha)
    fwer <- sum(vv) / replicates
    
    ## Type 1 Errors
    if (m0 > 0){
        vv <- indicatorFn(inp=mat[,H], type="sum", alpha=alpha)
        rr <- indicatorFn(inp=mat, type="sum", alpha=alpha)

        ## false discovery proportion
        E.Q <- sum(vv[rr > 0] / rr[rr > 0]) / length(rr[rr> 0])

        ## Probability of rejecting at least one hypothesis 
        P.R <- length(rr[rr > 0]) / length(rr)

        ## false discovery rate
        fdr <- E.Q * P.R

        ## 
        E.V <- sum(vv) / length(vv)
        
        ## pair-wise comparison error
        pcer <- E.V / m
        
    }   else {
        fdr <- 0
        E.Q <- 0
        E.V <- 0
        pcer <- 0
        P.R <- 0
    }

    ans <- list(mu = mu,
                parms = parms,
                n = n,
                m0 = m0,
                m = m,
                errfn = errfn,
                test = test,
                p.adjust.method = p.adjust.method,
                alternative = alternative,
                alpha = alpha,
                replicates = replicates,
                pcer = pcer,
                fdr = fdr,
                EQ = E.Q,
                PR = P.R,
                fwer = fwer,
                anypair = anyp,
                avepair = avep,
                allpair = allp,
                p0 = mat0,
                p = mat)
    class(ans) <- "powerPMCMR"
    ans
}
