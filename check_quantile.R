source("pav_algorithms.R")
source("generate_synthetic_data.R")


alpha <- c(0.25, 0.5, 0.75)

c <- 5
x <- list(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
           c(1, 2, 3, 4, 5, 5, 4, 3, 2, 1),
           c(1, 2, 3, 4, 5, 6, 7, 7, 9, 9),
           c(1, 2, 3, 4, 4, 5, 6, 7, 8, 8),
           c(1, 2, 3, 4, 5, 6))
y <- list(c(0, 0, 0, 1, 0, 0, 2, 0, 1, 3),
           c(5, 4, 3, 2, 1, 3, 4, 5, 6, 7),
           c(5, 4, 3, 2, 1, 3, 4, 5, 6, 7),
           c(1, 2, 3, 4, 5, 4, 3, 6, 7, 8),
           c(8, 2, 2, 3, 1, 1))

for (a in alpha) {
  for (i in 1:5) {
    x_rc1 <- my_cpp_isoreg_q(x[[i]], y[[i]], a)
    x_rc2 <- wrap_gpava_q(x[[i]], y[[i]], a)

    if (mean(abs(x_rc1 - x_rc2)) > 10^(-10)) {
      cat("!Error! Two different solutions: a =", a, "i =", i, "\n")
    }
  }
}

i <- 4
a <- 0.5
my_cpp_isoreg_q(x[[i]], y[[i]], a)
wrap_gpava_q(x[[i]], y[[i]], a)