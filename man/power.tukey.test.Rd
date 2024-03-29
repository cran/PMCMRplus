% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/power.tukey.test.R
\name{power.tukey.test}
\alias{power.tukey.test}
\title{Power Calculations for Balanced Tukey's
Multiple Comparison Test}
\source{
The Eqs. were taken from Lecture 5, \emph{Determining Sample Size},
Statistics 514, Fall 2015, Purdue University, IN, USA.
}
\usage{
power.tukey.test(n, groups, delta, within.var, sig.level = 0.05)
}
\arguments{
\item{n}{number of observations (per group)}

\item{groups}{number of groups}

\item{delta}{true difference in means}

\item{within.var}{within group variance}

\item{sig.level}{significance level (Type I error probability)}
}
\value{
Object of class \sQuote{\code{power.htest}},
a list of the arguments
(including the computed one) augmented with
\code{method} and \code{note} elements.
}
\description{
Compute average per-pair power of Tukey's test for
multiple comparison of means.
}
\details{
The function has implemented the following Eq.
to estimate average per-pair power for two-sided tests:

\deqn{
 1 - \beta = 1 - t(q_{\alpha v k}/\sqrt{2}, v, \mathrm{ncp}) +
  t(-q_{\alpha v k}/\sqrt{2}, v, \mathrm{ncp}),
}

with \eqn{q_{\alpha v k}} the upper \eqn{\alpha} quantile of
the studentised range distribution, with \eqn{v = k (n - 1)}
degree of freedom and \eqn{k} the number of groups;
and \eqn{t(. ~\mathrm{ncp})}
the probability function of the non-central student t-distribution
with non-centrality parameter

\deqn{
 \mathrm{ncp} = |\Delta| / \sqrt{s_{\mathrm{in}}^2 ~ 2 / n }.
}
}
\examples{
power.tukey.test(n = 11, groups = 5, delta = 30,
 within.var = 333.7)

## compare with t-test, Bonferroni-correction
power.t.test(n = 11, delta = 30, sd = sqrt(333.7),
sig.level = 0.05 / 10)

\dontrun{
powerMCTests(mu = c(rep(0,4), 30), n = 11,
 parms = list(mean = 0,sd = sqrt(333.7)),
 test = "tukeyTest")
}
}
\seealso{
\code{\link[stats]{TDist}} \code{\link[stats]{Tukey}}
\code{\link{powerMCTests}}
}
\keyword{htest}
