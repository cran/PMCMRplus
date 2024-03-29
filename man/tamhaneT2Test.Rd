% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tamhaneT2Test.R
\name{tamhaneT2Test}
\alias{tamhaneT2Test}
\alias{tamhaneT2Test.default}
\alias{tamhaneT2Test.formula}
\alias{tamhaneT2Test.aov}
\title{Tamhane's T2 Test}
\usage{
tamhaneT2Test(x, ...)

\method{tamhaneT2Test}{default}(x, g, welch = TRUE, ...)

\method{tamhaneT2Test}{formula}(formula, data, subset, na.action, welch = TRUE, ...)

\method{tamhaneT2Test}{aov}(x, welch = TRUE, ...)
}
\arguments{
\item{x}{a numeric vector of data values, a list of numeric data
vectors or a fitted model object, usually an \link[stats]{aov} fit.}

\item{\dots}{further arguments to be passed to or from methods.}

\item{g}{a vector or factor object giving the group for the
corresponding elements of \code{"x"}.
Ignored with a warning if \code{"x"} is a list.}

\item{welch}{indicates, whether Welch's approximate solution for
calculating the degree of freedom shall be used or, as usually,
\eqn{df = N - 2}. Defaults to \code{TRUE}.}

\item{formula}{a formula of the form \code{response ~ group} where
\code{response} gives the data values and \code{group} a vector or
factor of the corresponding groups.}

\item{data}{an optional matrix or data frame (or similar: see
\code{\link{model.frame}}) containing the variables in the
formula \code{formula}.  By default the variables are taken from
\code{environment(formula)}.}

\item{subset}{an optional vector specifying a
subset of observations to be used.}

\item{na.action}{a function which indicates what should happen when
the data contain \code{NA}s.  Defaults to \code{getOption("na.action")}.}
}
\value{
A list with class \code{"PMCMR"} containing the following components:
\describe{
\item{method}{a character string indicating what type of test was performed.}
\item{data.name}{a character string giving the name(s) of the data.}
\item{statistic}{lower-triangle matrix of the estimated
quantiles of the pairwise test statistics.}
\item{p.value}{lower-triangle matrix of the p-values for the pairwise tests.}
\item{alternative}{a character string describing the alternative hypothesis.}
\item{p.adjust.method}{a character string describing the method for p-value
adjustment.}
\item{model}{a data frame of the input data.}
\item{dist}{a string that denotes the test distribution.}
}
}
\description{
Performs Tamhane's T2 (or T2') all-pairs comparison test for normally distributed
data with unequal variances.
}
\details{
For all-pairs comparisons in an one-factorial layout
with normally distributed residuals but unequal groups variances
the T2 test (or T2' test) of Tamhane can be performed.
Let \eqn{X_{ij}} denote a continuous random variable
with the \eqn{j}-the realization (\eqn{1 \le j \le n_i})
in the \eqn{i}-th group (\eqn{1 \le i \le k}). Furthermore, the total
sample size is \eqn{N = \sum_{i=1}^k n_i}. A total of \eqn{m = k(k-1)/2}
hypotheses can be tested: The null hypothesis is
H\eqn{_{ij}: \mu_i = \mu_j ~~ (i \ne j)} is tested against the alternative
A\eqn{_{ij}: \mu_i \ne \mu_j} (two-tailed). Tamhane T2 all-pairs
test statistics are given by

\deqn{
 t_{ij} \frac{\bar{X}_i - \bar{X_j}}
 {\left( s^2_j / n_j + s^2_i / n_i \right)^{1/2}}, ~~
 (i \ne j)
}{%
 SEE PDF
}

with \eqn{s^2_i} the variance of the \eqn{i}-th group.
The null hypothesis is rejected (two-tailed) if

\deqn{
 \mathrm{Pr} \left\{ |t_{ij}| \ge t_{v_{ij}\alpha'/2} | \mathrm{H} \right\}_{ij} =
 \alpha.
}{%
 SEE PDF
}

T2 test uses Welch's approximate solution for
calculating the degree of freedom.

\deqn{
 v_{ij} = \frac{\left( s^2_i / n_i + s^2_j / n_j \right)^2}
 {s^4_i / n^2_i \left(n_i - 1\right) + s^4_j / n^2_j \left(n_j - 1\right)}.
}{%
 SEE PDF
}

T2' test applies the following approximation for the degree of freedom
\deqn{
 v_{ij} = n_i + n_j - 2
}{%
 SEE PDF
}

The p-values are computed from the \code{\link[stats]{TDist}}-distribution
and adjusted according to Dunn-Sidak.
\deqn{
 p'_{ij} = \min \left\{1, ~ (1 - (1 - p_{ij})^m)\right\}
}{%
 SEE PDF
}
}
\note{
T2 test is basically an all-pairs pairwise-t-test. Similar results
can be obtained with \code{pairwise.t.test(..., var.equal=FALSE, p.adjust.mehod = FALSE)}.

A warning message appears
in the modified T2' test, if none of in Tamhane (1979) given conditions
for nearly balanced
sample sizes and nearly balanced standard errors is true.

Thanks to Sirio Bolaños for his kind suggestion for adding T2' test
into this function.
}
\examples{
fit <- aov(weight ~ feed, chickwts)
shapiro.test(residuals(fit))
bartlett.test(weight ~ feed, chickwts) # var1 = varN
anova(fit)

## also works with fitted objects of class aov
res <- tamhaneT2Test(fit)
summary(res)
summaryGroup(res)
res

## compare with pairwise.t.test
WT <- pairwise.t.test(chickwts$weight,
                      chickwts$feed,
                      pool.sd = FALSE,
                      p.adjust.method = "none")
p.adj.sidak <- function(p, m) sapply(p, function(p) min(1, 1 - (1 - p)^m))
p.raw <- as.vector(WT$p.value)
m <- length(p.raw[!is.na(p.raw)])
PADJ <- matrix(ans <- p.adj.sidak(p.raw, m),
               nrow = 5, ncol = 5)
colnames(PADJ) <- colnames(WT$p.value)
rownames(PADJ) <- rownames(WT$p.value)
PADJ

## same without Welch's approximate solution
summary(T2b <- tamhaneT2Test(fit, welch = FALSE))

}
\references{
Tamhane, A. C. (1979) A Comparison of Procedures for Multiple Comparisons
of Means with Unequal Variances, \emph{Journal of the American
Statistical Association} \bold{74}, 471--480.
}
\seealso{
\code{\link{dunnettT3Test}} \code{\link{uryWigginsHochbergTest}}
}
\concept{parametric}
\keyword{htest}
