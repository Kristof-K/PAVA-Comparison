library(isotone)
library(Rcpp)

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

  # columns: groups, rows: (mean, count) of each group
  recal <- rbind(y[ord], rep(1, n))

  pool_groups <- function(i, recal, j = i) {
    I <- (i-1):j                                # merge group i:j into group i-1
    count <- sum(recal[2, I])                                  # new group count
    recal[1, i - 1] <- sum(recal[2, I] * recal[1, I]) / count  # new average
    recal[2, i - 1] <- count
    return(recal[, -(i:j), drop = F])   # delete merged group, keep dimensions
  }

  # find coinciding x-values
  no_change <- which(c(F, x[ord][2:n] == x[ord][1:(n-1)]))
  n_no_change <- length(no_change)
  if (n_no_change > 0) {
    # walk backwards, so that indices are valid despite pooling of recal
    j <- n_no_change
    for (i in n_no_change:1) {
      if (i > 1 && no_change[i-1] == no_change[i] - 1) {
        next      # decrement i until we are at the end of consecutive block
      } else {
        recal <- pool_groups(no_change[i], recal, no_change[j])
        j <- i - 1
      }
    }
  }

  i <- 2
  while (T) {
    if (recal[1, i] < recal[1, i-1]) {    # found violation --> pool
      recal <- pool_groups(i, recal)
      # pool with previous groups if now there are violations before
      for (j in (i-1):1) {
        if (j >= 2 && recal[1, j] < recal[1, j-1]) {
          recal <- pool_groups(j, recal)
          i <- i - 1      # we removed a group before i
        } else {
          break   # remainding groups are still increasing
        }
      }
    } else {     # go to next group
      i <- i + 1 # in "if" we deleted a group, so no need to increase i
    }
    if (i > ncol(recal))
      break   # iterated once through it all
  }
  return(rep(recal[1,], recal[2,]))
}

# compile R implementation
my_r_isoreg <- compiler::cmpfun(my_isoreg)
