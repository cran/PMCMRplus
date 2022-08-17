## See Chen (1999, p. 1237)
SD1p <- function(p) {
  ## adjustment function
  ## Dunn-Sidak adjustment
  p.adj <- function(p, k) {
    pmin(1, 1 - (1 - p) ^ k)
  }

  k <- length(p)

  ## output vector
  padj.j <- numeric(k)

  for (i in seq_len(k)) {
    if (i == 1) {
      ## first step
      ki <- k
    } else {
      ki <- m - 1
    }

    if (ki < 1) {
      break
    }

    o <- order(p[1:ki], decreasing = TRUE)
    ii <- seq_along(o)

    m <- which(ii == o[ki])

    for (j in m:ki) {
      padj.j[j] <- p.adj(p[j], ki)
    }

  }
  padj.j
}
