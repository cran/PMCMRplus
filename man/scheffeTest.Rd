% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scheffeTest.R
\name{scheffeTest}
\alias{scheffeTest}
\alias{scheffeTest.default}
\alias{scheffeTest.formula}
\alias{scheffeTest.aov}
\title{Scheffe's Test}
\usage{
scheffeTest(x, ...)

\method{scheffeTest}{default}(x, g, ...)

\method{scheffeTest}{formula}(formula, data, subset, na.action, ...)

\method{scheffeTest}{aov}(x, ...)
}
\arguments{
\item{x}{a numeric vector of data values, a list of numeric data
vectors or a fitted model object, usually an \link[stats]{aov} fit.}

\item{\dots}{further arguments to be passed to or from methods.}

\item{g}{a vector or factor object giving the group for the
corresponding elements of \code{"x"}.
Ignored with a warning if \code{"x"} is a list.}

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
Performs Scheffe's all-pairs comparisons test for normally distributed
data with equal group variances.
}
\details{
For all-pairs comparisons in an one-factorial layout
with normally distributed residuals and equal variances
Scheffe's test can be performed.
Let \eqn{X_{ij}} denote a continuous random variable
with the \eqn{j}-the realization (\eqn{1 \le j \le n_i})
in the \eqn{i}-th group (\eqn{1 \le i \le k}). Furthermore, the total
sample size is \eqn{N = \sum_{i=1}^k n_i}. A total of \eqn{m = k(k-1)/2}
hypotheses can be tested: The null hypothesis is
H\eqn{_{ij}: \mu_i = \mu_j ~~ (i \ne j)} is tested against the alternative
A\eqn{_{ij}: \mu_i \ne \mu_j} (two-tailed). Scheffe's all-pairs test
statistics are given by

\deqn{
 t_{ij} \frac{\bar{X}_i - \bar{X_j}}
 {s_{\mathrm{in}} \left(1/n_j + 1/n_i\right)^{1/2}}, ~~
 (i \ne j)
}{%
 SEE PDF
}

with \eqn{s^2_{\mathrm{in}}} the within-group ANOVA variance.
The null hypothesis is rejected if \eqn{t^2_{ij} > F_{v_{1}v_{2}\alpha}},
with \eqn{v_1 = k - 1, ~ v_2 = N - k} degree of freedom. The p-values
are computed from the \code{\link[stats]{FDist}} distribution.
}
\examples{
fit <- aov(weight ~ feed, chickwts)
shapiro.test(residuals(fit))
bartlett.test(weight ~ feed, chickwts)
anova(fit)

## also works with fitted objects of class aov
res <- scheffeTest(fit)
summary(res)
summaryGroup(res)
}
\references{
Bortz, J. (1993) \emph{Statistik für Sozialwissenschaftler}. 4. Aufl.,
Berlin: Springer.

Sachs, L. (1997) \emph{Angewandte Statistik}, New York: Springer.

Scheffe, H. (1953) A Method for Judging all Contrasts in the Analysis
of Variance, \emph{Biometrika} \bold{40}, 87--110.
}
\seealso{
\code{\link{FDist}}, \code{\link{tukeyTest}}
}
\concept{parametric}
\keyword{htest}
