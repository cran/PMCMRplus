% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/power.dunnett.test.R
\name{power.dunnett.test}
\alias{power.dunnett.test}
\title{Power Calculations for Balanced Dunnett's
Many-to-One Comparison Test}
\source{
The Eqs. were taken from Lecture 5, \emph{Determining Sample Size},
Statistics 514, Fall 2015, Purdue University, IN, USA.
}
\usage{
power.dunnett.test(n, groups, delta, within.var, sig.level = 0.05)
}
\arguments{
\item{n}{Number of observations (per group)}

\item{groups}{Number of groups (including control)}

\item{delta}{true difference in means}

\item{within.var}{Within group variance}

\item{sig.level}{Significance level (Type I error probability)}
}
\value{
Object of class \sQuote{\code{power.htest}},
a list of the arguments
(including the computed one) augmented with
\code{method} and \code{note} elements.
}
\description{
Compute average per-pair power of Dunnetts's multiple comparison
test with one control.
}
\details{
The function has implemented the following Eq.
to estimate average per-pair power for two-sided tests:

\deqn{
 1 - \beta = 1 - t( T_{\alpha \rho v}, v, \mathrm{ncp}) +
  t(-T_{\alpha \rho v}, v, \mathrm{ncp}),
}

with \eqn{T_{\alpha \rho v}} the two-sided
\eqn{\alpha} quantile of
the multivariate t-distribution, with \eqn{v = k (n - 1)}
degree of freedom, \eqn{k} the number of groups
and correlation matrix \eqn{\rho_{ij} = 0.5 ~ (i \neq j)}.

The non-centrality parameter for the non-central student t-distribution
is

\deqn{
 \mathrm{ncp} = |\Delta| / \sqrt{s_{\mathrm{in}}^2 ~ 2 / n }.
}
}
\note{
The results for power are seed depending.
}
\examples{
set.seed(113)
power.dunnett.test(n = 9, groups = 5, delta = 30,
 within.var = 333.7)

## compare with t-test, bonferroni corrected
power.t.test(n = 9, delta = 30, sd = sqrt(333.7),
sig.level = 0.05 / 4)

\dontrun{
## asymptotic Monte-Carlo power analysis
 set.seed(113)
 powerMCTests(mu = c(rep(0,4), 30), n = 9,
 parms = list(mean = 0, sd = sqrt(333.7)),
 test = "dunnettTest", alternative = "two.sided")
}
}
\seealso{
\code{\link[stats]{TDist}} \code{\link[mvtnorm]{qmvt}}
\code{\link{powerMCTests}}
}
\keyword{htest}
