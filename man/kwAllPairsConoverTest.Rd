% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kwAllPairsConoverTest.R
\name{kwAllPairsConoverTest}
\alias{kwAllPairsConoverTest}
\alias{kwAllPairsConoverTest.default}
\alias{kwAllPairsConoverTest.formula}
\title{Conover's All-Pairs Rank Comparison Test}
\usage{
kwAllPairsConoverTest(x, ...)

\method{kwAllPairsConoverTest}{default}(
  x,
  g,
  p.adjust.method = c("single-step", p.adjust.methods),
  ...
)

\method{kwAllPairsConoverTest}{formula}(
  formula,
  data,
  subset,
  na.action,
  p.adjust.method = c("single-step", p.adjust.methods),
  ...
)
}
\arguments{
\item{x}{a numeric vector of data values, or a list of numeric data
vectors.}

\item{\dots}{further arguments to be passed to or from methods.}

\item{g}{a vector or factor object giving the group for the
corresponding elements of \code{"x"}.
Ignored with a warning if \code{"x"} is a list.}

\item{p.adjust.method}{method for adjusting p values
(see \code{\link{p.adjust}}).}

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
Performs Conover's non-parametric all-pairs comparison test
for Kruskal-type ranked data.
}
\details{
For all-pairs comparisons in an one-factorial layout
with non-normally distributed residuals Conover's non-parametric test
can be performed. A total of \eqn{m = k(k-1)/2}
hypotheses can be tested. The null hypothesis
H\eqn{_{ij}: \mu_i(x) = \mu_j(x)} is tested in the two-tailed test
against the alternative
A\eqn{_{ij}: \mu_i(x) \ne \mu_j(x), ~~ i \ne j}.

If \code{p.adjust.method == "single-step"} the p-values are computed
from the studentized range distribution. Otherwise,
the p-values are computed from the t-distribution using
any of the p-adjustment methods as included in \code{\link{p.adjust}}.
}
\examples{
## Data set InsectSprays
## Global test
kruskalTest(count ~ spray, data = InsectSprays)

## Conover's all-pairs comparison test
## single-step means Tukey's p-adjustment
ans <- kwAllPairsConoverTest(count ~ spray, data = InsectSprays,
                             p.adjust.method = "single-step")
summary(ans)

## Dunn's all-pairs comparison test
ans <- kwAllPairsDunnTest(count ~ spray, data = InsectSprays,
                             p.adjust.method = "bonferroni")
summary(ans)

## Nemenyi's all-pairs comparison test
ans <- kwAllPairsNemenyiTest(count ~ spray, data = InsectSprays)
summary(ans)
}
\references{
Conover, W. J, Iman,  R. L. (1979) \emph{On multiple-comparisons
procedures}, Tech. Rep. LA-7677-MS, Los Alamos Scientific Laboratory.
}
\seealso{
\code{\link[stats]{Tukey}}, \code{\link[stats]{TDist}},
\code{\link[stats]{p.adjust}}, \code{\link{kruskalTest}},
\code{\link{kwAllPairsDunnTest}}, \code{\link{kwAllPairsNemenyiTest}}
}
\concept{kruskalranks}
\keyword{htest}
\keyword{nonparametric}
