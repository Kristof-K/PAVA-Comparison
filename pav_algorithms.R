library(isotone)
library(Rcpp)
library(monotone)
library(dplyr)
library(data.table)

sourceCpp("pava_mean.cpp")

# Define PAVA implementations:
# Functions should receive a x and y vector and return the transformed vector

# wrapper of gpava function of package isotone
wrap_gpava <- function(x, y) {
  # ties secondary : average y values for ties in x
  res <- gpava(x, y, ties = "secondary")
  return(res$x[order(x)])
}

# wrapper for isoreg function
wrap_isoreg <- function(x, y) {
  res <- isoreg(x, y)
  return(res$yf)
}

wrap_monotone <- function(x, y) {
  xy <- data.table(x = x, y = y) %>% arrange(x, desc(y))
  return(monotone(xy$y))
}

# wrapper for PAV algorithm implemented in cpp
my_cpp_isoreg <- function(x, y) {
  res <- pav_mean_cpp(x, y)
  return(res[, 2])
}

# PAV algorithm implemented in  R
my_isoreg <- function(x, y) {
  stopifnot(length(x) == length(y), is.numeric(x), is.numeric(y))
  n <- length(x)
  ord <- order(x)
  x_sorted <- x[ord]
  g_count <- rep(0, n)
  g_val <- rep(0, n)

  pool_with_prev <- function(i) {
    new_count <- g_count[i] + g_count[i - 1]
    # update mean and count: <<- access bigger scope
    g_val[i-1] <<- (g_count[i] * g_val[i] + g_count[i-1] * g_val[i-1]) / new_count
    g_count[i-1] <<- new_count
  }

  # now walk through y-values and pool adjacent violators
  g_curr <- 1
  i <- 1
  while (i <= n) {
    # create new group
    g_val[g_curr] <- y[ord[i]]
    g_count[g_curr] <- 1
    while (i < n && x_sorted[i] == x_sorted[i + 1]) {    # coninciding x-values
      i <- i + 1
      g_val[g_curr] <- g_val[g_curr] + y[ord[i]]
      g_count[g_curr] <- g_count[g_curr] + 1
    }
    g_val[g_curr] <- g_val[g_curr] / g_count[g_curr]      # average

    if (g_curr > 1 && g_val[g_curr] < g_val[g_curr - 1]) {   # violations ?
      pool_with_prev(g_curr)   # pool with group before
      # new violations before? --> pool them!
      while ((g_curr > 2) && (g_val[g_curr - 1] < g_val[g_curr - 2])) {
	    g_curr <- g_curr - 1  	# decrease beforehand as we have checked for j-1 and j-2
        pool_with_prev(g_curr)
      }
    } else {
      g_curr <- g_curr + 1     # introduce next group
    }
    i <- i + 1
  }
  # build output together
  return(rep(g_val, g_count)[1:n])
}

# compile R implementation
my_r_isoreg <- compiler::cmpfun(my_isoreg)
