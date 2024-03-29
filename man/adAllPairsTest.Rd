% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adAllPairsTest.R
\name{adAllPairsTest}
\alias{adAllPairsTest}
\alias{adAllPairsTest.default}
\alias{adAllPairsTest.formula}
\title{Anderson-Darling All-Pairs Comparison Test}
\usage{
adAllPairsTest(x, ...)

\method{adAllPairsTest}{default}(x, g, p.adjust.method = p.adjust.methods, ...)

\method{adAllPairsTest}{formula}(
  formula,
  data,
  subset,
  na.action,
  p.adjust.method = p.adjust.methods,
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

\item{p.adjust.method}{method for adjusting p values (see \code{\link{p.adjust}}).}

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
Performs Anderson-Darling all-pairs comparison test.
}
\details{
For all-pairs comparisons in an one-factorial layout
with non-normally distributed residuals Anderson-Darling's
all-pairs comparison test can be used. A total of \eqn{m = k(k-1)/2}
hypotheses can be tested. The null hypothesis
H\eqn{_{ij}: F_i(x) = F_j(x)} is tested in the two-tailed test
against the alternative
A\eqn{_{ij}: F_i(x) \ne F_j(x), ~~ i \ne j}.

This function is a wrapper function that sequentially
calls \code{adKSampleTest} for each pair.
The calculated p-values for \code{Pr(>|T2N|)}
can be adjusted to account for Type I error multiplicity
using any method as implemented in \code{\link{p.adjust}}.
}
\examples{
adKSampleTest(count ~ spray, InsectSprays)

out <- adAllPairsTest(count ~ spray, InsectSprays, p.adjust="holm")
summary(out)
summaryGroup(out)

}
\references{
Scholz, F.W., Stephens, M.A. (1987) K-Sample Anderson-Darling Tests.
\emph{Journal of the American Statistical Association} \bold{82}, 918--924.
}
\seealso{
\code{\link{adKSampleTest}}, \code{\link{adManyOneTest}},
\code{\link[kSamples]{ad.pval}}.
}
\keyword{htest}
\keyword{nonparametric}
