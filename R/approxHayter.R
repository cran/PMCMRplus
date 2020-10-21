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
approxHayter <- function(k, df) {

  ## use the CritTable values from Hayter
  ## for balanced design
  nrows <- nrow(TabCrit$hayter.h005)
  kk <- as.numeric(colnames(TabCrit$hayter.h005))
  dft <- as.numeric(rownames(TabCrit$hayter.h005))

  ## check for kk
  if (k > max(kk) | k < min(kk)) stop("No critical values for k = ", k)

  ## interpolate
  yh <- TabCrit$hayter.h005[, paste0(k)]
  hCrit <- approx(x = dft, y = yh, xout = df)$y

  hCrit
}
