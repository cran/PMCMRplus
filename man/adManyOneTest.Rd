% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adManyOneTest.R
\name{adManyOneTest}
\alias{adManyOneTest}
\alias{adManyOneTest.default}
\alias{adManyOneTest.formula}
\title{Anderson-Darling Many-To-One Comparison Test}
\usage{
adManyOneTest(x, ...)

\method{adManyOneTest}{default}(x, g, p.adjust.method = p.adjust.methods, ...)

\method{adManyOneTest}{formula}(
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

\item{p.adjust.method}{method for adjusting
p values (see \code{\link{p.adjust}}).}

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
Performs Anderson-Darling many-to-one comparison test.
}
\details{
For many-to-one comparisons (pairwise comparisons with one control)
in an one-factorial layout with non-normally distributed
residuals Anderson-Darling's non-parametric test can be performed.
Let there be \eqn{k} groups including the control,
then the number of treatment levels is \eqn{m = k - 1}.
Then \eqn{m} pairwise comparisons can be performed between
the \eqn{i}-th treatment level and the control.
H\eqn{_i: F_0 = F_i} is tested in the two-tailed case against
A\eqn{_i: F_0 \ne F_i, ~~ (1 \le i \le m)}.

This function is a wrapper function that sequentially
calls \code{adKSampleTest} for each pair.
The calculated p-values for \code{Pr(>|T2N|)}
can be adjusted to account for Type I error inflation
using any method as implemented in \code{\link{p.adjust}}.
}
\note{
Factor labels for \code{g} must be assigned in such a way,
that they can be increasingly ordered from zero-dose
control to the highest dose level, e.g. integers
\{0, 1, 2, ..., k\} or letters \{a, b, c, ...\}.
Otherwise the function may not select the correct values
for intended zero-dose control.

It is safer, to i) label the factor levels as given above,
and to ii) sort the data according to increasing dose-levels
prior to call the function (see \code{\link{order}}, \code{\link{factor}}).
}
\examples{
## Data set PlantGrowth
## Global test
adKSampleTest(weight ~ group, data = PlantGrowth)

##
ans <- adManyOneTest(weight ~ group,
                             data = PlantGrowth,
                             p.adjust.method = "holm")
summary(ans)
}
\references{
Scholz, F.W., Stephens, M.A. (1987) K-Sample Anderson-Darling Tests.
\emph{Journal of the American Statistical Association} \bold{82}, 918--924.
}
\seealso{
\code{\link{adKSampleTest}}, \code{\link{adAllPairsTest}},
\code{\link[kSamples]{ad.pval}}.
}
\keyword{htest}
\keyword{nonparametric}
