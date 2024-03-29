% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/steelTest.R
\name{steelTest}
\alias{steelTest}
\alias{steelTest.default}
\alias{steelTest.formula}
\title{Steel's Many-to-One Rank Test}
\usage{
steelTest(x, ...)

\method{steelTest}{default}(x, g, alternative = c("greater", "less"), ...)

\method{steelTest}{formula}(
  formula,
  data,
  subset,
  na.action,
  alternative = c("greater", "less"),
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

\item{alternative}{the alternative hypothesis. Defaults to \code{greater}}

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
A list with class \code{"osrt"} that contains the following components:

\describe{
\item{method}{a character string indicating what type of test was performed.}
\item{data.name}{a character string giving the name(s) of the data.}
\item{statistic}{the estimated statistic(s)}
\item{crit.value}{critical values for \eqn{\alpha = 0.05}.}
\item{alternative}{a character string describing the alternative hypothesis.}
\item{parameter}{the parameter(s) of the test distribution.}
\item{dist}{a string that denotes the test distribution.}
}
There are print and summary methods available.
}
\description{
Performs Steel's non-parametric many-to-one comparison
test for Wilcox-type ranked data.
}
\details{
For many-to-one comparisons (pairwise comparisons with one control)
in an one-factorial balanced layout with non-normally distributed
residuals Steels's non-parametric single-step test can be performed.
Let there be \eqn{k} treatment levels (excluding the control),
then \eqn{k} pairwise comparisons can be performed between
the \eqn{i}-th treatment level and the control.
H\eqn{_i: \theta_0 = \theta_i} is tested in the one-tailed case (less) against
A\eqn{_i: \theta_0 > \theta_i, ~~ (1 \le i \le k)}.

For each control - treatment level the data are ranked in increasing order.
The ranksum \eqn{R_i} for the \eqn{i}-th treatment level is compared
to a critical \eqn{R} value and is significantly(\eqn{p = 0.05}) less,
if \eqn{R_i \le R}. For the \code{alternative = "greater"} the sign is changed.

The function does not return p-values. Instead the critical \eqn{R}-values
as given in the tables of USEPA (2002) for \eqn{\alpha = 0.05} (one-sided, less)
are looked up according to the balanced sample sizes (\eqn{n}) and the order number of the
dose level (\eqn{i}).
}
\note{
Steel's Many-to-One Rank test is only applicable for balanced designs and
directional hypotheses. An error message will occur, if the design is unbalanced.
In the current implementation, only one-sided tests on
the level of \eqn{\alpha = 0.05} can be performed.
}
\section{Source}{

The critical rank sum values were taken from Table E.5 of USEPA (2002).

USEPA (2002) \emph{Short-term Methods for Estimating the
Chronic Toxicity of Effluents and Receiving
Waters to Freshwater Organisms}, 4th edition, EPA-821-R-02-013.
}

\examples{
## Example from Sachs (1997, p. 402)
x <- c(106, 114, 116, 127, 145,
110, 125, 143, 148, 151,
136, 139, 149, 160, 174)
g <- gl(3,5)
levels(g) <- c("0", "I", "II")

## Steel's Test
steelTest(x ~ g)


## Example from USEPA (2002):
## Reproduction data from a Ceriodaphnia dubia
## 7-day chronic test to several concentrations
## of effluent. Dose level 50\% is excluded.
x <- c(20, 26, 26, 23, 24, 27, 26, 23, 27, 24,
13, 15, 14, 13, 23, 26, 0, 25, 26, 27,
18, 22, 13, 13, 23, 22, 20, 22, 23, 22,
14, 22, 20, 23, 20, 23, 25, 24, 25, 21,
9, 0, 9, 7, 6, 10, 12, 14, 9, 13,
rep(0,10))
g <- gl(6, 10)
levels(g) <- c("Control", "3\%", "6\%", "12\%", "25\%", "50\%")

## NOEC at 3\%, LOEC at 6\%
steelTest(x ~ g, subset = g != "50\%", alternative = "less")



}
\references{
Steel, R. G. D. (1959) A multiple comparison rank sum test:
treatments versus control, \emph{Biometrics} \bold{15}, 560--572.
}
\seealso{
\code{\link[stats]{wilcox.test}}, \code{\link[stats]{pairwise.wilcox.test}},
\code{\link{manyOneUTest}}, \code{\link{flignerWolfeTest}},
\code{\link{shirleyWilliamsTest}}, \code{\link{kwManyOneDunnTest}},
\code{\link{kwManyOneNdwTest}}, \code{\link{kwManyOneConoverTest}},
\code{\link{print.osrt}}, \code{\link{summary.osrt}}
}
\keyword{htest}
