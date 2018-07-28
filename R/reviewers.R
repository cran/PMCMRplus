#' @name reviewers
#' @docType data
#' @title
#'  Reviewers
#'
#' @description
#'  9 reviewers (blocks) assigned ranks to 4 objects (groups).
#'
#' @format
#'  The format is a 9 x 4 Matrix with Friedman type rankings:
#' \describe{
#'  \item{rows}{reviewers, 1, 2, \ldots, 9}
#'  \item{columns}{groups, A, B, \ldots, D}
#' }
#'
#' @source
#'  Sachs (1997), p. 671 ff.
#'
#' @references
#'  Sachs, L. (1997) \emph{Angewandte Statistik}, New York: Springer.
#'
#' @examples
#' data(reviewers)
#' friedmanTest(reviewers)
#' pageTest(reviewers)
#' frdAllPairsExactTest(reviewers, p.adjust = "bonferroni")
#'
#' @keywords datasets
NULL
