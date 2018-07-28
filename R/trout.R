#' @name trout
#' @docType data
#' @title
#' Data from a Dose-Response Experiment with Trouts
#' @description
#'  This data set contains results from a dose-response experiment with trouts.
#'  The experiment was conducted with five doses of 10, 25, 60, 150 and
#'  1000 ppm, respectively, plus a zero-dose control. The response is
#'  trout weight in mg.
#'
#' @format
#'  A data frame with 65 observations on the following 5 variables.
#'  \describe{
#'    \item{CONC}{a numeric vector of dose concentration in ppm}
#'    \item{DOSE}{a factor with levels \code{1} \code{2}
#'      \code{3} \code{4} \code{5} \code{6}}
#'    \item{REPA}{a factor with levels \code{1} \code{2}}
#'    \item{REPC}{a factor with levels \code{1} \code{2}}
#'    \item{Y}{a numeric vector of trout weight in mg}
#'  }
#'
#' @source
#'  ENV/JM/MONO(2006)18/ANN, page 113.
#'
#' @references
#'  OECD (ed. 2006) \emph{Current approaches in the statistical analysis
#' of ecotoxicity data: A guidance to application - Annexes}. OECD Series
#' on testing and assessment, No. 54, (ENV/JM/MONO(2006)18/ANN).
#'
#' @keywords datasets
NULL
