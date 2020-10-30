## adjust.hCrit.R
## Part of the R package: PMCMRplus
##
## Copyright (C) 2020 Thorsten Pohlert
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## A copy of the GNU General Public License is available at
## http://www.r-project.org/Licenses/

## aux function for all test that base on Hayter's h-distribution
adjust.hCrit <- function(hCrit, is.balanced = TRUE) {
  if (!is.balanced) {
    warning("Design is un-balanced. Using h/sqrt(2) approximation.")
    hCrit <- hCrit / sqrt(2)
  }
  return(hCrit)
}
