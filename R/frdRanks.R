## internal function frdRanks
## part of package PMCMRplus
## Copyright (C) 2014-2019 Thorsten Pohlert
## License: GPL-3

## Purpose: Friedman type ranking
## taken from stats::friedman.test

frdRanks <- function(y, groups, blocks) {
  DNAME <- deparse(substitute(y))
  if (is.matrix(y)) {
    groups <- factor(c(col(y)))
    blocks <- factor(c(row(y)))

    GRPNAMES <- colnames(y)
    ROWNAMES <- rownames(y)
  }
  else {
    if (anyNA(groups) || anyNA(blocks))
      stop("NA's are not allowed in 'groups' or 'blocks'")
    if (any(diff(c(
      length(y), length(groups), length(blocks)
    )) != 0L))
      stop("'y', 'groups' and 'blocks' must have the same length")
    DNAME <- paste0(DNAME,
                    ", ",
                    deparse(substitute(groups)),
                    " and ",
                    deparse(substitute(blocks)))
    if (any(table(groups, blocks) != 1))
      stop("not an unreplicated complete block design")
    groups <- factor(groups)
    blocks <- factor(blocks)
    ## Need to ensure consistent order of observations within
    ## blocks.
    o <- order(groups, blocks)
    y <- y[o]
    groups <- groups[o]
    blocks <- blocks[o]

    GRPNAMES <- levels(groups)
    ROWNAMES <- levels(blocks)
  }

  k <- nlevels(groups)
  ## <FIXME split.matrix>
  y <- matrix(unlist(split(c(y), blocks)), ncol = k, byrow = TRUE)
  y <- y[complete.cases(y),]
  n <- nrow(y)
  r <- t(apply(y, 1L, rank))
  ## <FIXME split.matrix>

  colnames(r) <- GRPNAMES
  rownames(r) <- ROWNAMES

  inDF = data.frame(x = as.vector(y),
                    g = rep(groups,
                            length(as.vector(y)) /
                              length(groups)),
                    b = rep(blocks,
                            length(as.vector(y)) /
                              length(blocks)))
  levels(inDF$g) <- GRPNAMES
  levels(inDF$b) <- ROWNAMES
  ans <- list(r = r,
              inDF = inDF,
              DNAME = DNAME)
  return(ans)
}
