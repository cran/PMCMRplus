#  Copyright (C) 2020 Thorsten Pohlert
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#
#' @title Summarize an osrt Object
#' @description
#' Summarize an object of class \emph{osrt}.
#' @method summary osrt
#' @aliases summary.osrt
#' @param object an object of class \code{"osrt"}.
#' @param \dots further arguments. Currenly ignored.
#' @keywords methods
#'
#' @seealso
#' \code{\link{print.osrt}}.
#' @export
summary.osrt <- function(object, ...)
{
  critVal <- as.numeric(object$crit.value)
  stat <- as.numeric(object$statistic)
  dist <- object$dist
  if (grepl(pattern = "Steel",
            x = object$method)) {
    dec <- ifelse(stat <= critVal , "reject", "accept")
  } else {
    dec <- ifelse(stat > critVal, "reject", "accept")
  }

  critDist <- paste0(dist, "-crit")

  if (!is.matrix(object$statistic)) {
    if (grepl(pattern = "Hayter's",
              x = object$method)) {
      H0 <- switch(
        object$alternative,
        greater = paste("Mean(xi) - Mean(xj) <= 0"),
        less = paste("Mean(xi) - Mean(xj) >= 0")
      )
    } else {
      H0 <- switch(
        object$alternative,
        greater = paste("Med(xi) - Med(xj) <= 0"),
        less = paste("Med(xi) - Med(xj) >= 0")
      )
    }
  } else {
    grp1 <- as.numeric(c(col(object$statistic)))
    grp2 <- as.numeric(c(row(object$statistic)))
    cnam <- colnames(object$statistic)
    rnam <- rownames(object$statistic)
    H0 <- switch(
      object$alternative,
      greater = paste(rnam[grp2], "-", cnam[grp1], "<=", "0"),
      less = paste(rnam[grp2], "-", cnam[grp1], ">=", "0")
    )
    ok <- !is.na(stat)
    stat <- stat[ok]
    H0 <- H0[ok]
    dec <- dec[ok]
  }
  dist <- paste0(dist, "-value")

  cat("\n\t", object$method, "\n\n")
  cat("data: ", object$data.name, "\n")
  if (!is.null(object$alternative)) {
    cat("alternative hypothesis: ", object$alternative, "\n")
  }

  paramName <- names(object$parameter)



  ## pretty names

  if (length(paramName) == 2) {
    suppressWarnings(
      expr =
        xdf <- data.frame(
          STATISTIC = round(stat, 3),
          PARAM1 = object$parameter[1],
          PARAM2 = object$parameter[2],
          CRITDIST = round(critVal, 3),
          DECISION = dec,
          ALPHA = 0.05
        )
    )
    names(xdf) <- c(dist,
                    paramName[1],
                    paramName[2],
                    critDist,
                    "decision",
                    "alpha")

  } else {
    suppressWarnings(
      expr =
        xdf <- data.frame(
          STATISTIC = round(stat, 3),
          PARAM1 = object$parameter[1],
          #  PARAM2 = object$parameter[2],
          CRITDIST = round(critVal, 3),
          DECISION = dec,
          ALPHA = 0.05
        )
    )
    names(xdf) <- c(dist,
                    paramName[1],
                    critDist,
                    "decision",
                    "alpha")
  }
  rownames(xdf) <- H0

  cat("\n")
  cat("H0\n")
  ##
  print(xdf)
  cat("---\n")
  invisible(object)
}

#' @title osrt Printing
#' @description
#' \code{print.osrt} is the \emph{osrt} method of the generic
#' \code{\link{print}} function which prints its argument
#' and returns it \emph{invisibly} (via \code{\link{invisible}(x)}).
#' @param x an object used to select a method.
#' @param \ldots further arguments. Currently ignored.
#' @aliases print.osrt
#' @method print osrt
#' @keywords print
#' @seealso summary.osrt
#' @export
print.osrt <-
  function(x, ...)
  {
    summary.osrt(x)
    invisible(x)
  }
