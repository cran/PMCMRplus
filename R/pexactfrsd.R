#
# Thorsten Pohlert
# v.01 modified on 22. March 2017
# ---------------------------------------------------------------------------------- #
#                                                                                    #
# Exact p-values for pairwise comparison of Friedman rank sums, with application to  #
# comparing classifiers                                                              #
# Eisinga, Heskes, Pelzer & Te Grotenhuis                                            # 
# BMC Bioinformatics, January 2 2017                                                 #
#                                                                                    #
# ---------------------------------------------------------------------------------- # 
#
# pexactfrsd
#
# Description 
#
# The function pexactfrsd computes the exact p-value of the absolute value of rank 
# sum difference d of a pair of Friedman rank sums, in an analysis of k treatments 
# and n blocks. It also offers the possibility to compute the mass point probability, 
# the number of compositions, and the cumulative number of compositions of the 
# absolute value of rank sum difference d, by specifying an optional argument.
#
# Usage
#
# pexactfrsd(d,k,n,option)
#
# Arguments
#
# d       absolute value of rank sum difference. 
# k       number of treatments.
# n       number of blocks.
# option  character string indicating the desired statistic: "pvalue", "probability",
#         "no of compositions", or "cumulative no of compositions". If the character 
#         string is not provided in the function call, the function returns the exact
#         p-value (default).  
#
# Values
#
# The potential range of rank sum difference d is 0 to n(k-1) inclusive. The number 
# of treatments k should be at least 2, and the number of blocks n at least 1. 
# Depending on the option specified, the function computes the following:
#
# "pvalue"                         returns P(D>=d;k,n)
# "probability"                    returns P(D=d;k,n)
# "no of compositions"             returns W(D=d;k,n)
# "cumulative no of compositions"  returns W(D>=d;k,n).
#
# Details
#
# The function pexactfrsd is an implementation of the algorithm provided in Eisinga 
# et al. [1]. The function requires the R package Rmpfr [2] to be installed. In the 
# script below the maximal precision, in bits, is set at (precBits =) 2048, which  
# is sufficient even for rather large values of n (100, say). 
#
# Note
#
# It is important to note that the results pertain to the absolute value of d. The 
# rank sum difference distribution with positive and negative d values is symmetric 
# around 0. Hence the probability mass to the left of d=0 may be folded over, 
# producing a discrete distribution of non-negative d, ranging from 0 to n(k-1).
# The probability P(D=d;k,n) of all d except d=0 in the distribution of non-negative 
# d is doubled relative to the probability in the symmetric distribution, so that 
# they sum to unity. The same doubling goes for the (cumulative) number of 
# compositions. Consequently, the number of compositions of d=n(k-1), for example, 
# equals 2 and not 1, as is the case for the symmetric discrete distribution with 
# support d=[-n(k-1),n(k-1)].  
#
# Examples
#
# Example 1 following the R function with argument values d=k=n=100 returns a 
# p-value of 0.8085251. The result is obtained in about 0.5 secs, using an Intel 
# Core i7-3520M CPU @ 2.9Ghz. Example 2 calculates for k=n=10 the exact p-value of 
# d values ranging from 1 to n(k-1), and subsequently plots and prints out the 
# results. Example 3 generates the statistics presented in the example application 
# in Additional file 4 of Eisinga et al. [1]. Example 4 computes the p-value and the 
# mid p-value for k=n=5 in Table 4 of Eisinga et al. [1].                   
#
# References
#
# [1] Eisinga, Rob, Tom Heskes, Ben Pelzer, Manfred Te Grotenhuis (2017). Exact 
#     p-values for pairwise comparison of Friedman rank sums, with application 
#     to comparing classifiers, BMC Bioinformatics 
# [2] Maechler, Martin. 2015. Rmpfr: R MPFR Multiple Precision Floating-Point 
#     Reliable, Version 0.6-0, December 4 2015,  
#     https://cran.r-project.org/web/packages/Rmpfr/index.html.
#


# tp: library Rmpfr is going to be attached via a namespace.
#library(Rmpfr)


pexactfrsd <- function(d,k,n,option) {
 if (any(n < 1))                     stop("n out-of-bounds: min = 1")
 if (any(k < 2))                     stop("k out-of-bounds: min = 2")
 if (any(d < 0) || any(d > n*(k-1))) stop("d out-of-bounds: min,max = 0,n(k-1)")
 if (missing(option)) {option = "pvalue"}
 result <- 0
 for (h in 0:n) {
    sum1 <- chooseZ(n,h)/mpfr((pow.bigz(k,h) * pow.bigz(1-k,n)), precBits = 2048)
    sum2 <- 0
    for (s in 0:h) {
        if (any(k*s-d-h >= 0)) {
          if (option == "pvalue" || option == "cumulative no of compositions"){
             sum2 <- sum2 + (-1)^s * chooseZ(2*h,h+s) * chooseZ(k*s-d+h,k*s-d-h)}
          if (option == "probability" || option == "no of compositions"){
             sum2 <- sum2 + (-1)^s * chooseZ(2*h,h+s) * chooseZ(k*s-d+h-1,k*s-d-h)}
        }
    }
    result <- result + sum1 * sum2
 }
if (any(d == 0) & option== "pvalue") return(1)
if (any(d != 0) & option== "pvalue") return(as.numeric(2*result))
if (any(d == 0) & option== "probability") return(as.numeric(result))
if (any(d != 0) & option== "probability") return(as.numeric(2*result))
if (any(d != 0) & 
   (option== "no of compositions" || option== "cumulative no of compositions") )            
   return(round(2*result*mpfr(pow.bigz(k*(k-1),n), precBits = 2048)))
if (any(d == 0) & option== "no of compositions")            
   return(round(result*mpfr(pow.bigz(k*(k-1),n), precBits = 2048)))
if (any(d == 0) & option== "cumulative no of compositions") 
   return(round(mpfr(pow.bigz(k*(k-1),n), precBits = 2048)))
}


