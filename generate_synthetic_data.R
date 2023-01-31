# create increasing normal data
gen_inc_norm <- function(n) {
  x <- runif(n, 0, 10)
  y <- rnorm(n, x, x / 3 + 1)
  return(data.frame(x = x, y = y))
}

# create u-shaped normal data
gen_u_norm <- function(n) {
  x <- runif(n, 0, 10)
  y <- rnorm(n, abs(x - 5), (x - 5)^2 * 2 / 25 + 1)
  return(data.frame(x = x, y = y))
}

# create u-shaped normal data
gen_inv_u_norm <- function(n) {
  x <- runif(n, 0, 10)
  y <- rnorm(n, 5 - abs(x - 5), (x - 5)^2 * 2 / 25 + 1 * (x < 8))
  return(data.frame(x = x, y = y))
}

# create normal data with ties in x-values
gen_ties_norm <- function(n) {
  x <- runif(n, 0, 10)
  # on average 10 coinciding x-values
  x <- round(x, pmax(0, floor(log(n, base = 10)) - 2))
  y <- rnorm(n, sqrt(x), x / 3 + 1)
  return(data.frame(x = x, y = y))
}
