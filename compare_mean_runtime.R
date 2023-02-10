library(lubridate)
library(ggplot2)

source("pav_algorithms.R")
source("generate_synthetic_data.R")

compare_runtime <- function(pava_list, gen_data_list, n_list, plot_name,
                            my_colors, my_shapes, B = 1) {
  # - pava_list:
  # list of pava implementations following format used in pav_algorithms.R
  # - gen_data_list:
  # list of data generating functions following format used in generate_synthetic_data.R
  # - n_list:
  # vector of sizes of data
  # - plot_name:
  # name of the plot comparing runtimes
  # - my_colors:
  # named vector assigning algorithms in pava_list colors
  # - my_shapes:
  # named vector assigning algorithms in pava_list shape values
  # - B
  # number of times used to average runtimes
  collect_results <- data.frame()

  for (n in n_list) {
    cat("\n", n, ":", sep = "")
    for (dgp in names(gen_data_list)) {
      cat("\n-", dgp, ":", sep = "")
      cmp_res <- NULL
      data <- gen_data_list[[dgp]](n)
      for (a in names(pava_list)) {
        cat(" ", a, sep = "")
        time <- 0

        # average runtimes
        for (b in 1:B) {
          start <- now()
          res <- pava_list[[a]](data$x, data$y)
          end <- now()
          time <- time + as.numeric(seconds(end - start))
        }

        # check whether algorithms determine equal solutions
        if (!is.null(cmp_res) && mean(abs(res - cmp_res)) > 10e-10)
          cat("[Error] (", a, "): Varying solution ", mean(abs(res - cmp_res)),
              "\n", sep = "")

        collect_results <- rbind(
          collect_results,
          data.frame(n = n, T = time / B, Data = dgp, Algorithm = a)
        )
        cmp_res <- res    # store res to compare with next algorithm
      }
      cmp_res <- NULL     # forget res as we are moving on to new data
    }
  }
  ggplot(collect_results, aes(x = n, y = T, color = Algorithm)) +
    facet_wrap(~Data, scales = "free_y") +
    geom_point(aes(shape = Algorithm)) +
    geom_line(size = 0.5) +
    scale_color_manual(name = "Algorithm", values = my_colors,
                       breaks = names(pava_list)) +
    scale_shape_manual(name = "Algorithm", values = my_shapes,
                       breaks = names(pava_list)) +
    xlab("n") +
    ylab("runtime [s]") +
    ggtitle("PAV Implementations: Runtime vs. Input Size") +
    theme_bw() +
    theme(legend.justification = c(0, 1), legend.position = c(0.01, 0.99),
          legend.background = element_blank(),
          legend.box.background = element_rect(color = "gray"),
          panel.grid.major = element_line(size = 0.05),
          panel.grid.minor = element_line(size = 0.05),
          strip.background = element_blank())
  ggsave(plot_name, width = 250, height = 150, unit = "mm")
}


# data generating processes
l_dgp <- list(IncreasingNormal = gen_inc_norm, UShapedNormal = gen_u_norm,
              IncreasingNormalTies = gen_ties_norm,
              InvUShapedNormal = gen_inv_u_norm)
num_avg <- 3

# visualiize data
for (f in l_dgp) {
  with(f(1000), plot(x, y))
}

set.seed(999)

# We can use GPAVA only for small n values
l_algo <- list(GPAVA = wrap_gpava, ISOREG = wrap_isoreg, MY_R = my_r_isoreg,
               MY_CPP = my_cpp_isoreg, MONOTONE = wrap_monotone)
data_size <- c(100, 500, 1000, 5000, 10000)
cols <- setNames(RColorBrewer::brewer.pal(length(l_algo), name = "Dark2"),
                 names(l_algo))
shapes <- setNames(1:length(l_algo), names(l_algo))
name <- "figures/analysis1.pdf"
compare_runtime(l_algo, l_dgp, data_size, name, cols, shapes, B = num_avg)

# GPAVA ran into memory issues for n = 10^5

# If we skip GPAVA we can use a bit larger data, but then the R implementation
# slows down
l_algo <- list(ISOREG = wrap_isoreg, MY_R = my_r_isoreg, MY_CPP = my_cpp_isoreg)
data_size <- c(1000, 5000, 10000, 50000)
name <- "figures/analysis2_2.pdf"
compare_runtime(l_algo, l_dgp, data_size, name, cols, shapes, B = num_avg)

# For even larger data sets, we can only look at isoreg and the cpp implementation
l_algo <- list(ISOREG = wrap_isoreg, MY_CPP = my_cpp_isoreg,
               MONOTONE = wrap_monotone)
data_size <- c(10000, 50000, 100000, 500000, 10^6, 5*10^6, 10^7)
name <- "figures/analysis3.pdf"
compare_runtime(l_algo, l_dgp, data_size, name, cols, shapes, B = num_avg)

# For even larger data sets, we can only look at isoreg and the cpp implementation
l_algo <- list(MY_CPP = my_cpp_isoreg, MONOTONE = wrap_monotone)
data_size <- c(10000, 50000, 100000, 500000, 10^6, 5*10^6, 10^7)
name <- "figures/analysis4.pdf"
compare_runtime(l_algo, l_dgp, data_size, name, cols, shapes, B = num_avg)
