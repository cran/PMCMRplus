## get ties for kruskalTest and
## Sachs proposed Nemenyi-Test for tie correction
gettiesKruskal <- function(x) {
  n <- length(x)
  t <- table(x)
  C <- 1 - sum(t^3 - t) / (n^3 - n)
  C <- min(1, C)
  return(C)
}
