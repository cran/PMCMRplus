# Taken from function 'pHayStonLSA'
# of package NSM3.
# (C) 2020 Grant Schneider, Eric Chicken, Rachel Becvarik
# GPL-2
#
# Function to compute the upper tail probability of the
# Hayter-Stone W (= h_{k v}, with v = \infty)
# asymptotic distribution for a given cutoff.
#
#' @importFrom stats dnorm pnorm
pHS <- function (h, k, delta = 0.001)
{
    our.grid <- seq(-8, 8, delta)
    len <- length(our.grid)
    init.grid <- pnorm(our.grid + h)
    if (k > 2) {
        for (i in 3:k) {
            new.grid <- cumsum(init.grid * dnorm(our.grid) *
                delta) + init.grid * (pnorm(our.grid + h) - pnorm(our.grid))
            init.grid <- new.grid
        }
    }
    out <- 1 - sum(dnorm(our.grid) * init.grid * delta)
    out <- min(1, out)
    return(out)
}
