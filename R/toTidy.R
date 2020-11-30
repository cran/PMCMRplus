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
#' @name toTidy
#' @title Convert a PMCMR or osrt Object to a Data.Frame
#' @description
#' The functions converts a list object of class \code{"PMCMR"}
#' or \code{"osrt"} into a data.frame.
#'
#' @param mod an object of class \code{"PMCMR"}, \code{"trendPMCMR"} or \code{"osrt"}.
#' @param \ldots further arguments. Currently ignored.
#'
#' @return
#' A data.frame.
#'
#' @author
#' Indrajeet Patil (via email, 2020-1022),
#' modified by Thorsten Pohlert
#'
#' @examples
#' res <- tukeyTest(weight ~ Diet, data = ChickWeight, subset = Time == 21)
#' toTidy(res)
#'
#' @keywords utilities
#' @importFrom stats na.omit
#' @export
toTidy <- function(mod, ...) {
  matrix_to_tidy <- function(m, col_name = "value") {
    result <-
      data.frame(
        group1 = rep(rownames(m), each = ncol(m)),
        group2 = rep(colnames(m), times = nrow(m)),
        col3 = as.numeric(t(m)),
        stringsAsFactors = FALSE
      )

    names(result)[3] <- col_name
    na.omit(result)
  }

  # get rid of any \n and \t in method string
  regex1 <- function(x) {
    gsub(pattern = "(\\n|\\t)",
         replacement = "",
         x = x)
  }
  ## get rid of any multiple whitespaces
  regex2 <- function(x) {
    gsub(pattern = "( ){2,}",
         replacement = " ",
         x = x)
  }

  if (inherits(mod, "PMCMR") | inherits(mod, "trendPMCMR")) {
    # combine all components of the object in a single dataframe
    METH <- regex1(mod$method)
    METH <- regex2(METH)
    METH <- ifelse(length(METH) > 1,
                   paste(METH, collapse = ""),
                   METH)

    ans <- cbind(
      merge(
        x = matrix_to_tidy(mod$statistic, col_name = "statistic"),
        y = matrix_to_tidy(mod$p.value, col_name = "p.value"),
        by = c("group1", "group2")
      ),
      data.frame(alternative = ifelse(
        is.null(mod$alternative),
        "two.sided",
        mod$alternative
      )),
      data.frame(method = METH),
      data.frame(distribution = mod$dist),
      data.frame(p.adjust.method = mod$p.adjust.method)
    )
    return(ans)

  } else if (inherits(mod, "osrt")) {
    METH <- regex1(mod$method)
    METH <- regex2(METH)
    METH <- ifelse(length(METH) > 1,
                   paste(METH, collapse = ""),
                   METH)

    ans <- cbind(
      matrix_to_tidy(mod$statistic, col_name = "statistic"),
      data.frame(crit.value = mod$crit.value),
      data.frame(alternative = ifelse(
        is.null(mod$alternative),
        "two.sided",
        mod$alternative
      )),
      data.frame(method = METH),
      data.frame(distribution = mod$dist)
    )
    return(ans)
  } else {
    stop("Must be an object of class ", sQuote("PMCMR"), " or ", sQuote("osrt"),"!")
  }
}
