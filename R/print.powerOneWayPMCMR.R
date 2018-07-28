## print.powerOneWayPMCMR.R
##
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
##

#' @name print.powerOneWayPMCMR
#' @title PowerOneWayPMCMR Printing
#' @description
#' \code{print.powerOneWayPMCMR} is the
#' \emph{powerOneWayPMCMR} method of the generic
#' \code{\link{print}} function which prints its argument
#' and returns it \emph{invisibly} (via \code{\link{invisible}(x)}).
#' @param x an object used to select a method.
#' @param \ldots further arguments. Currently ignored.
#' @aliases print.powerOneWayPMCMR
#' @method print powerOneWayPMCMR
#' @keywords print
#' @export
print.powerOneWayPMCMR <- function(x, ...)
{
    METH <- paste0("\n\tPower-simulation for ", x$test,
                   "\n\t with alternative = ", x$alternative)
    message(METH)
    message(paste0("\nInput summary:"))
    message(paste0("Tested linear model: X_ij = mu_i + epsilon_ij"))
    message(paste0("with epsilon_ij ~ ", x$errfn, "\n"))
    k <- length(x$mu)
    cohensf <- NULL
    if(x$errfn == "Normal"){
        if (length(x$parms$sd) == k){
            sdE <- x$parms$sd
        } else if (length(x$parms$sd) == 1){
            sdE <- rep(x$parms$sd, k)
            MU <- mean(x$mu)
            if(length(x$n) == k){
                SSA <- sum(x$n * (x$mu - MU)^2)
                SSE <- x$parms$sd^2 * (sum(x$n) - 1)
            } else if (length(x$n) == 1){
                SSA <- sum(x$n * (x$mu - MU)^2)
                SSE <- x$parms$sd^2 * (k * x$n - 1)
            }
            R2 <- SSA / (SSA + SSE)
            cohensf <- sqrt(R2 / (1 - R2))
        } else if (all.equal(x$parms$sd)) {
            sdE <- rep(unique(x$parms$sd),k)
            MU <- mean(x$mu)
            if(length(x$n) == k){
                SSA <- sum(x$n * (x$mu - MU)^2)
                SSE <- x$parms$sd^2 * (sum(x$n) - 1)
            } else if (length(x$n) == 1){
                SSA <- sum(x$n * (x$mu - MU)^2)
                SSE <- x$parms$sd^2 * (k * x$n - 1)
            }
            R2 <- SSA / (SSA + SSE)
            cohensf <- sqrt(R2 / (1 - R2))
        } else { ## arbitrary error variances
            g <- unlist(sapply(1:k, function(i) rep(i, x$n[i])))
            g <- as.factor(g)
            sdE <- tapply(x$parms$sd, g, unique)
        }

        ans <- data.frame(mu = x$mu,
                          n = x$n,
                          meanE = rep(unique(x$parms$mean), k),
                          sdE = sdE)
        print(ans)

    } else if (x$errfn == "Lognormal"){

        if (length(x$parms$sdlog) == k){
            sdlogE = x$parms$sdlog
        } else if (length(x$parms$sdlog) == 1) {
            sdlogE <- rep(unique(x$parms$sdlog),k)
        } else if (all.equal(x$parms$sdlog)) {
            sdlogE <- rep(unique(x$parms$sdlog),k)
        }
        ans <- data.frame(mu = x$mu,
                          n = x$n,
                          meanlogE = rep(unique(x$parms$meanlog), k),
                          sdlogE = sdlogE)
        print(ans)
    } else if (x$errfn == "Cauchy"){

        if (length(x$parms$scale) == k){
            scaleE = x$parms$scale
        } else if (length(x$parms$scale) == 1) {
            scaleE <- rep(unique(x$parms$scale),k)
        } else if (all.equal(x$parms$scale)) {
            scaleE <- rep(unique(x$parms$scale),k)
        }

        ans <- data.frame(mu = x$mu,
                          n = x$n,
                          locationE = rep(unique(x$parms$location), k),
                          scaleE = scaleE)
        print(ans)
    }  else if (x$errfn == "Weibull"){

        if (length(x$parms$scale) == k){
            scaleE = x$parms$scale
        } else if (length(x$parms$scale) == 1) {
            scaleE <- rep(unique(x$parms$scale),k)
        } else if (all.equal(x$parms$scale)) {
            scaleE <- rep(unique(x$parms$scale),k)
        }
        ans <- data.frame(mu = x$mu,
                          n = x$n,
                          shapeE = rep(unique(x$parms$shape), k),
                          scaleE = scaleE)
        print(ans)
    }  else if (x$errfn == "Exponential"){
        ans <- data.frame(mu = x$mu,
                          n = x$n,
                          rateE = rep(unique(x$parms$rate), k))
        print(ans)
    }  else if (x$errfn == "TDist" | x$errfn == "Chisquare"){
        ans <- data.frame(mu = x$mu,
                          n = x$n,
                          df = rep(unique(x$parms$df), k))
        print(ans)
    }
    message("\n")
    inpDF <- data.frame(Value = (c(formatC(x$alpha, digits=2,format="f"),
                                   formatC(c(x$replicates),
                                           format="d"))))
    rownames(inpDF) <- c("Nominal alpha:",
                         "Monte-Carlo replications:")
    print(inpDF)
    if(!is.null(cohensf)){
        message("\nEffect sizes")
        efDF <- data.frame(Value = c(formatC(cohensf, digits=3,format="f"),
                                     formatC(R2, digits=3,format="f")))
        rownames(efDF) <- c("Cohen's f", "R-squared")
        print(efDF)
    }
    message(paste0("\nResults for empirical Type I and Type II Error"))
    t1err <- data.frame(a = format.pval(c(x$fwer, 1 - x$fwer)),
                        b = format.pval(c(x$power, 1 - x$power)))
    rownames(t1err) <- c("Reject H0", "Accept H0")
    names(t1err) <- c("H0 True", "HA True")
    print(t1err)
}
