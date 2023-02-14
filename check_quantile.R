source("pav_algorithms.R")
source("generate_synthetic_data.R")

library(lubridate)
library(ggplot2)

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

# check equality with synthetic data
n_list <- c(50, 100, 500, 1000)

l_dgp <- list(IncreasingNormal = gen_inc_norm, UShapedNormal = gen_u_norm,
              IncreasingNormalTies = gen_ties_norm,
              InvUShapedNormal = gen_inv_u_norm)

for (a in alpha) {
  for (n in n_list) {
    for (gen_data in names(l_dgp)) {
      data <- l_dgp[[gen_data]](n)
      x_rc1 <- my_cpp_isoreg_q(data$x, data$y, a)
      x_rc2 <- wrap_gpava_q(data$x, data$y, a)

      if (mean(abs(x_rc1 - x_rc2)) > 10^(-10)) {
        cat("!Error! Two different solutions: a =", a, "n =", n, "dgp =", gen_data,"\n")
      }
      if (any(x_rc2[-1] - x_rc2[-n] < 0)) {
        cat(" - GPava determined non-increasing solution\n")
      }
    }
  }
}

set.seed(1)
data <- gen_inv_u_norm(20)
x_rc1 <- my_cpp_isoreg_q(data$x, data$y, 0.75)
x_rc2 <- wrap_gpava_q(data$x, data$y, 0.75)

t <- cbind(x_rc1, x_rc2)
View(t)

# are there ties?
sum(table(data$x) >= 2)
ord <- order(data$x)

xy <- data.frame(x = data$x[ord], y = data$y[ord])
View(xy)

# check all types of quantiles
for (t in 1:9) {
  cat(quantile((data$y[ord])[5:20], probs = 0.75, type = t), " ")
}

# Check how fast solutions are determined
n_list <- 10^(3:6)
B <- 3
collect_res <- data.frame()

for (a in alpha) {
  for (gen_data in names(l_dgp)) {
    cat("\n", gen_data, ":")
    for (n in n_list) {
      cat(n, "")
      time <- 0
      data <- l_dgp[[gen_data]](n)

      # average runtimes
      for (b in 1:B) {
        start <- now()
        res <- my_cpp_isoreg_q(data$x, data$y, a)
        end <- now()
        time <- time + interval(start, end) / seconds(1)
      }
      collect_res <- rbind(collect_res, data.frame(D = gen_data, N = n, A = a, T = time / B))
    }
  }
}

ggplot(collect_res) +
  geom_line(aes(x = N, y = T, color = D, linetype = factor(A))) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_discrete(name = "Synthetic Data") +
  scale_linetype_discrete(name = expression(alpha)) +
  xlab("n") +
  ylab("runtime [s]") +
  ggtitle("PAVA Quantile Computation Time") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.box = "vertical")
ggsave("figures/quantile_runtime_5.pdf", width = 250, height = 180, unit = "mm")
