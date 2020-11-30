## HStat
## computes Kruskal-Wallis statistic H
## r and g must be vectors of the same length
HStat <- function(r, g) {
  ni <- tapply(!is.na(r), g, length)
  N <- sum(ni)

  H <- (12 / (N * (N + 1))) *
    sum(tapply(r, g, "sum") ^ 2 / ni) - 3 * (N + 1)
  H
}
